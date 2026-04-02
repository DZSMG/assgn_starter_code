// Task 4: Created Renderer.cpp to optimize main.cpp

#include "Renderer.h"
#include "Constants.h"
#include <fstream>
#include <cmath>
#include <algorithm>

// Moved from main.cpp. Keep the internal to this file using static
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

// Moved implementation of public function from main.cpp
void renderSnapshot(const Body &B, const std::string &file, double zoom) {
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

        // Refactor to a coloured rendering
        if (B[i].mass > 0.5 * SOLAR_MASS) {
            // Make the central star pop by drawing a 3x3 pixel square instead of 1 pixel
            for (int dx = -1; dx <= 1; dx++) {
                for (int dy = -1; dy <= 1; dy++) {
                    plotPixel(im, px + dx, py + dy, 255, 220, 40); // Bright Yellow
                }
            }
        } else {
            // --- DYNAMIC COLOR GRADIENT ---
            // 1. Calculate the 2D distance from the galactic center
            double distance = std::sqrt(rx * rx + ry * ry);

            // 2. Normalize the distance (0.0 at center, 1.0 at screen edge)
            double normalizedDist = std::min(1.0, distance / view_half);

            // 3. Map to RGB: Inner = Cyan/White, Outer = Deep Pink/Purple
            uint8_t r = static_cast<uint8_t>(255 * normalizedDist);         // Red increases outward
            uint8_t g = static_cast<uint8_t>(200 * (1.0 - normalizedDist)); // Green is high in center
            uint8_t b = 255;                                                // Constant high blue

            plotPixel(im, px, py, r, g, b);
        }
    }

    writePPM(im, WIDTH, HEIGHT, file);
}
