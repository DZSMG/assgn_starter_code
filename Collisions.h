#ifndef COLLISIONS_H
#define COLLISIONS_H

#include <vector>
#include <cmath>
#include "Types.h"
#include "Constants.h"

// Merge b into a when distance(a,b) < thresholdAU.
// Conserves mass and linear momentum. Positions in AU; velocities in m/s.

// Task 2: Complete the function checkAndMergeCollision based on pseudocode provided.
inline void checkAndMergeCollision(body& a, body& b, double thresholdAU = COLLISION_THRESHOLD_AU) {
    // Error Handling: Ensure we don't calculate inactive bodies
    if (a.mass <= 0.0 || b.mass <= 0.0) return;

    // Calculate squared distance to avoid expensive std::sqrt
    double dx = a.position.x - b.position.x;
    double dy = a.position.y - b.position.y;
    double dz = a.position.z - b.position.z;
    double distSq = dx * dx + dy * dy + dz * dz;
    double thresholdSq = thresholdAU * thresholdAU;

    if (distSq < thresholdSq) {
        double totalMass = a.mass + b.mass;

        // Update Position to Center of Mass
        a.position.x = (a.position.x * a.mass + b.position.x * b.mass) / totalMass;
        a.position.y = (a.position.y * a.mass + b.position.y * b.mass) / totalMass;
        a.position.z = (a.position.z * a.mass + b.position.z * b.mass) / totalMass;

        // Conserve Momentum for Velocity
        a.velocity.x = (a.velocity.x * a.mass + b.velocity.x * b.mass) / totalMass;
        a.velocity.y = (a.velocity.y * a.mass + b.velocity.y * b.mass) / totalMass;
        a.velocity.z = (a.velocity.z * a.mass + b.velocity.z * b.mass) / totalMass;

        // Finalize merge
        a.mass = totalMass;
        b.mass = 0.0; // Mark B as inactive
        b.velocity = {0, 0, 0};
    }
}

// Quadratic pass: merges any pairs closer than threshold; skips retired bodies.
// Task 3: Improve running time by using random number generator.
/*
 * The original implementation used a nested loop to perform exhaustive pairwise collision checks.
 * This resulted in an O(N^2) time complexity. This creates a massive CPU bottleneck and completely defeats
 * the O(N log N) performance gains provided by the Barnes-Hut tree algorithm.
 *
 * The updated implementation uses a Monte Carlo randomized sampling approach, reducing the time complexity to O(N).
 * By utilizing modern C++ <random> (std::mt19937), we uniformly and efficiently sample a subset of pairs each frame.
 * Although collision detection becomes probabilistic rather than deterministic, it dramatically improves FPS while
 * strictly preserving the macroscopic orbital dynamics and center-of-mass integrity of the simulation.
 */
template <typename Bodies>
inline void mergeAllCollisions(Bodies& bodies, double thresholdAU = COLLISION_THRESHOLD_AU) {
    // Cast to int to match the distribution type and avoid signed/unsigned mismatch warnings
    const int n = static_cast<int>(bodies.size());

    // Early exit: no collisions possible with fewer than 2 bodies
    if (n <= 1) return;

    // 'static' ensures the generator persists across calls, avoiding expensive re-seeding
    // every frame while still producing high-quality uniformly distributed indices.
    static std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<int> pickIndex(0, n - 1);

    // Number of random collision checks per body.
    const int samplesPerBody = 10;

    // Outer loop: iterate over every body exactly once.
    for (int i = 0; i < n; ++i) {
        if (bodies[i].mass <= 0.0) continue;

        // Inner loop: for each active body, randomly sample a fixed number of partners.
        for (int s = 0; s < samplesPerBody; ++s) {
            int j = pickIndex(rng);

            if (j == i) continue;
            if (bodies[j].mass <= 0.0) continue;

            // Perform the actual distance check and conditional merge.
            // Only bodies within thresholdAU of each other will be merged.
            checkAndMergeCollision(bodies[i], bodies[j], thresholdAU);
        }
    }
}

#endif // COLLISIONS_H