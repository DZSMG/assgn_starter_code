// Task 4: Created Simulation.cpp to optimize main.cpp

#include "Simulation.h"
#include "Collisions.h"
#include <chrono>
#include <cmath>
#include <algorithm>

// Moved the implementations of the helper functions from main.cpp
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

// Moved the implementations of the declared functions
StarSpec makeMilkyWaySpec(vec3 pos, vec3 vel,
                                 int numBodies,
                                 double bulgeMass,
                                 double diskInnerAU,
                                 double diskOuterAU,
                                 double diskThicknessAU,
                                 double diskMassFrac) {
    StarSpec s{bulgeMass, pos, vel, diskInnerAU, diskOuterAU, diskThicknessAU, numBodies, diskMassFrac};
    return s;
}

void initialiseStarsAndBodies_Exponential(Body &bb,
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
Octant makeDynamicRootOctant(const Body &bb, double marginAU) {
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
void simulateStep(Body &bb) {
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