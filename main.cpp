#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <vector>
#include <utility>

#include "Constants.h"
#include "Types.h"
#include "Octant.h"
#include "Bhtree.h"


struct Body {
    body *ptr{nullptr};  // pointer to array of body
    size_t n{0};
    size_t cap{0};

    Body() = default;  // default constructor

    // keep simple: no accidental copies
    Body(const Body &o) = delete;
    Body &operator=(const Body &o) = delete;

    Body(Body &&o) noexcept // move-only type
        : ptr(o.ptr), n(o.n), cap(o.cap) {
        o.ptr = nullptr;
        o.n = o.cap = 0;
    }

    void reserve(size_t newCap) {
        if (newCap <= cap) return;
        body *np = new body[newCap];
        if (ptr) {
            // copy existing elements
            for (size_t i = 0; i < n; ++i) np[i] = ptr[i];
        }
        delete[] ptr;
        ptr = np;
        cap = newCap;
    }


    // This is working perfectly now! Don't mess with it!
    Body &operator=(Body &&o) noexcept {
        if (this != &o) {
            delete[] ptr;
            ptr = o.ptr;
            n = o.n;
            cap = o.cap;
            o.ptr = nullptr;
            o.n = o.cap = 0;
        }
        return *this;
    }

    // Everything else below also works perfectly!
    body *begin() { return ptr; }
    body *end() { return ptr + n; } // always ptr+n (works even if ptr==nullptr && n==0)
    const body *begin() const { return ptr; }
    const body *end() const { return ptr + n; }

    void resize(size_t newSize) {
        if (newSize > cap) reserve(newSize);
        n = newSize;
    }

    void clear() { n = 0; }

    bool empty() const { return n == 0; }
    size_t size() const { return n; }

    body &operator[](size_t i) { return ptr[i]; }
    const body &operator[](size_t i) const { return ptr[i]; }

    void push_back(const body &b) {
        if (n + 1 > cap) reserve(cap ? cap * 2 : 8);
        ptr[n++] = b;
    }
};

// ---------------------------------------------------------------------
#include "Collisions.h"

// ---------------- Utilities ---------------------
static vec3 barycentre(const Body &B) {
    long double mx = 0, my = 0, mz = 0, m = 0;
    for (const auto &b: B) {
        if (b.mass <= 0.0) continue;
        mx += (long double) b.mass * b.position.x;
        my += (long double) b.mass * b.position.y;
        mz += (long double) b.mass * b.position.z;
        m += (long double) b.mass;
    }
    if (m == 0) return {0, 0, 0};
    return {(double) (mx / m), (double) (my / m), (double) (mz / m)};
}

static inline int toPixel(double v, int sizePx, double view_half) {
    double t = (v / std::max(1e-9, view_half)) * 0.5 + 0.5;
    if (t < 0.0) t = 0.0;
    if (t > 1.0) t = 1.0; // clamp here
    return static_cast<int>(t * (sizePx - 1) + 0.5);
}


static void writePPM(const std::vector<uint8_t> &rgb, int W, int H, const std::string &file) {
    std::ofstream os(file, std::ios::binary);
    os << "P6\n" << W << " " << H << "\n255\n";
    os.write(reinterpret_cast<const char *>(rgb.data()), static_cast<std::streamsize>(rgb.size()));
}

static void plotPixel(std::vector<uint8_t> &im, int x, int y, uint8_t r, uint8_t g, uint8_t b) {
    if (x < 0 || y < 0 || x >= WIDTH || y >= HEIGHT) return;
    const int idx = 3 * (y * WIDTH + x);
    im[idx + 0] = (uint8_t) std::min(255, im[idx + 0] + (int) r);
    im[idx + 1] = (uint8_t) std::min(255, im[idx + 1] + (int) g);
    im[idx + 2] = (uint8_t) std::min(255, im[idx + 2] + (int) b);
}

static void renderSnapshot(const Body &B, const std::string &file, double zoom = 1.0) {
    std::vector<uint8_t> im(WIDTH * HEIGHT * 3, 0);

    // Camera centre: barycentre
    vec3 center = barycentre(B);

    // Zoom: view_half = window radius in AU
    const double view_half = VIEW_HALF_AU * zoom;

    for (size_t i = 0; i < B.size(); ++i) {
        if (B[i].mass <= 0.0) continue;
        double rx = B[i].position.x - center.x;
        double ry = B[i].position.y - center.y;

        if (std::abs(rx) > view_half || std::abs(ry) > view_half) continue;

        int px = toPixel(rx, WIDTH, view_half);
        int py = toPixel(ry, HEIGHT, view_half);

        if (B[i].mass > 0.5 * SOLAR_MASS) {
            // Draw star as a single yellow pixel (can change to 3x3 if preferred)
            plotPixel(im, px, py, 255, 220, 40);
        } else {
            plotPixel(im, px, py, 180, 200, 255);
        }
    }

    writePPM(im, WIDTH, HEIGHT, file);
}

// ---- Helpers for realistic g0 disks ------------------------------------
static inline double sampleExponentialRadius(double Rd, std::mt19937 &r) {
    std::uniform_real_distribution<double> U(0.0, 1.0);
    double u = std::max(1e-12, U(r) * U(r));
    return -Rd * std::log(u);
}

static inline double enclosedMassExponentialDisk(double R, double Md, double Rd) {
    const double x = R / std::max(Rd, 1e-12);
    return Md * (1.0 - std::exp(-x) * (1.0 + x));
}

static inline double randNormal(std::mt19937 &r, double sigma) {
    if (sigma <= 0.0) return 0.0;
    static thread_local std::normal_distribution<double> N(0.0, 1.0);
    return sigma * N(r);
}

static StarSpec makeMilkyWaySpec(vec3 pos, vec3 vel,
                                 int numBodies,
                                 double bulgeMass,
                                 double diskInnerAU,
                                 double diskOuterAU,
                                 double diskThicknessAU,
                                 double diskMassFrac) {
    StarSpec s{bulgeMass, pos, vel, diskInnerAU, diskOuterAU, diskThicknessAU, numBodies, diskMassFrac};
    return s;
}

static void initialiseStarsAndBodies_Exponential(Body &bb,
                                                 const std::vector<StarSpec> &s0,
                                                 const DiskKinematics &kin) {
    size_t total = s0.size();
    for (const auto &s: s0) total += (s.numBodies > 0 ? size_t(s.numBodies) : 0);
    bb.clear();
    bb.resize(total);

    std::mt19937 r(static_cast<uint32_t>(
        std::chrono::high_resolution_clock::now().time_since_epoch().count()));
    std::uniform_real_distribution<double> Uang(0.0, 2.0 * PI);
    std::uniform_real_distribution<double> Uz(-0.5, 0.5);

    size_t writeIdx = 0;
    for (const auto &s: s0) {
        body star{};
        star.position = s.position;
        star.velocity = s.velocity;
        star.accel = {0, 0, 0};
        star.mass = s.mass;
        bb[writeIdx++] = star;
    }

    // Disks
    for (const auto &s: s0) {
        if (s.numBodies <= 0 || s.mass <= 0.0) continue;
        const double Rin = std::max(1e-6, s.diskInnerAU);
        const double Rout = std::max(Rin + 1e-6, s.diskOuterAU);
        const double Rd = std::min(std::max(kin.RdAU, 0.1), Rout);
        const double Md = s.diskMassFrac * s.mass;
        const double mEach = Md / double(s.numBodies);

        for (int i = 0; i < s.numBodies; ++i) {
            double R = sampleExponentialRadius(Rd, r);
            int tries = 0;
            while ((R < Rin || R > Rout) && tries++ < 16) R = sampleExponentialRadius(Rd, r);
            R = std::clamp(R, Rin, Rout);

            const double ang = Uang(r);
            const double z = Uz(r) * s.diskThicknessAU;

            vec3 pos{s.position.x + R * std::cos(ang), s.position.y + R * std::sin(ang), s.position.z + z};

            const double MencBulge = s.mass;
            const double MencDisk = enclosedMassExponentialDisk(R, Md, Rd);
            const double MencTotal = std::max(0.0, MencBulge + MencDisk);

            const double vCirc = std::sqrt(G * MencTotal / std::max(R * AU, 1e-6));
            const double sinA = std::sin(ang), cosA = std::cos(ang);
            vec3 vtan{vCirc * sinA, -vCirc * cosA, 0.0};
            vtan.x += randNormal(r, kin.sigmaPhi);
            vtan.y += randNormal(r, kin.sigmaPhi);
            vec3 vrad{randNormal(r, kin.sigmaR) * cosA, randNormal(r, kin.sigmaR) * sinA, 0.0};
            vec3 vz{0.0, 0.0, randNormal(r, kin.sigmaZ)};

            body b{};
            b.position = pos;
            b.velocity = vtan + vrad + vz + s.velocity;
            b.accel = {0, 0, 0};
            b.mass = mEach;
            bb[writeIdx++] = b;
        }
    }

    // Remove net momentum once
    vec3 p{0, 0, 0};
    double msum = 0.0;
    for (const auto &b: bb) {
        if (b.mass > 0.0) {
            p = p + b.velocity * b.mass;
            msum += b.mass;
        }
    }
    if (msum > 0.0) {
        vec3 vcm = p / msum;
        for (auto &b: bb) if (b.mass > 0.0) b.velocity = b.velocity - vcm;
    }
}

// Dynamic r0 for t0
static Octant makeDynamicRootOctant(const Body &bb, double marginAU) {
    double minx = 1e9, miny = 1e9, minz = 1e9;
    double maxx = -1e9, maxy = -1e9, maxz = -1e9;
    for (const auto &b: bb) {
        if (b.mass <= 0.0) continue;
        minx = std::min(minx, b.position.x);
        miny = std::min(miny, b.position.y);
        minz = std::min(minz, b.position.z);
        maxx = std::max(maxx, b.position.x);
        maxy = std::max(maxy, b.position.y);
        maxz = std::max(maxz, b.position.z);
    }
    vec3 mid{0.5 * (minx + maxx), 0.5 * (miny + maxy), 0.5 * (minz + maxz)};
    double span = std::max({maxx - minx, maxy - miny, maxz - minz}) + marginAU;
    span = std::max(span, 2.0 * marginAU);
    return Octant{mid.x, mid.y, mid.z, span};
}

// One simulation stp with dynamic r0 and selectable drift fix
static void simulateStep(Body &bb) {
    // Build dynamic r0 box (AU)
    Octant r0 = makeDynamicRootOctant(bb, std::max(1.0, SYSTEM_SIZE_AU));
    Bhtree t0(std::move(r0));

    // Insert live bb
    for (auto &b: bb) {
        if (b.mass > 0.0 && t0.octant().contains(b.position)) {
            t0.insert(&b);
        }
    }

    // Reset accelerations and compute forces via t0
    for (auto &b: bb) {
        if (b.mass <= 0.0) continue;
        b.accel = {0, 0, 0};
        if (t0.octant().contains(b.position)) t0.interactInTree(&b);
    }

    // Note: Collisions are effectively off with tiny threshold
    mergeAllCollisions(bb, COLLISION_THRESHOLD_AU);

    // Semi-implicit Euler integration
    for (auto &b: bb) {
        if (b.mass <= 0.0) continue;
        b.velocity = b.velocity + b.accel * DT; // m/s
        b.position = b.position + (b.velocity * DT) / AU; // AU
    }

    // Drift fix selection
    switch (DRIFT_FIX) {
        case ANCHOR_FIRST_STAR:
            if (!bb.empty()) {
                bb[0].position = {0, 0, 0};
                bb[0].velocity = {0, 0, 0};
            }
            break;
        case REMOVE_BARYCENTRE_VELOCITY_EACH_STEP: {
            vec3 p{0, 0, 0};
            double msum = 0.0;
            for (const auto &b: bb)
                if (b.mass > 0.0) {
                    p = p + b.velocity * b.mass;
                    msum += b.mass;
                }
            if (msum > 0.0) {
                vec3 vcm = p / msum;
                for (auto &b: bb) if (b.mass > 0.0) b.velocity = b.velocity - vcm;
            }
            break;
        }
        case REMOVE_BARYCENTRE_VELOCITY:
        case NO_FIX:
        default:
            break;
    }
}

// TODO this file is getting too long to work on...
int main() {
    // this will execute the deleteImgs PowerShell script
    system("powershell -ExecutionPolicy Bypass -File \"../deleteImgs.ps1\"");
    // system("bash deleteImgs.bash");

    // Build a single Milky Way-like galaxy
    std::vector<StarSpec> s0;
    s0.push_back(makeMilkyWaySpec({0, 0, 0}, {0, 0, 0},
                                  NUM_BODIES,
                                  1.2 * SOLAR_MASS,
                                  0.3,
                                  8.0,
                                  0.2,
                                  0.25));

    DiskKinematics kin; // default dispersions and Rd
    Body bb;
    initialiseStarsAndBodies_Exponential(bb, s0, kin);

    std::cout << "Running Barnes-Hut n-body (" << bb.size()
            << " bb, " << STEPS << " steps)..." << std::endl;

    for (int s = 0; s < STEPS; ++s) {
        simulateStep(bb);
        std::ostringstream fname;
        fname << "images/frame_" << std::setw(3) << std::setfill('0') << s << ".ppm";
        renderSnapshot(bb, fname.str(), RENDER_ZOOM);
        if ((s + 1) % 50 == 0) {
            std::cout << " stp " << (s + 1) << "/" << STEPS << "\r" << std::flush;
        }
    }

    std::cout << "\nDone. Wrote " << STEPS << " frames as PPM images." << std::endl;

    // this will execute the createVideo PowerShell script
    system("powershell -ExecutionPolicy Bypass -File \"../createVideo.ps1\"");
    // system("bash createVideo.bash");
    return 0;
}
