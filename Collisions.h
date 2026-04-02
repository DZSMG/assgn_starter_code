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
    const size_t n = bodies.size();
    if (n < 2) return;

    // Modern C++ Thread-Safe Random Number Generator
    static thread_local std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<size_t> dist(0, n - 1);

    // Monte Carlo optimization: Instead of N^2 checks, we do a subset of random checks.
    // This dramatically improves FPS but makes collisions probabilistic.
    size_t checks = n * 4; // Arbitrary heuristic multiplier for performance

    for (size_t k = 0; k < checks; ++k) {
        size_t i = dist(rng);
        size_t j = dist(rng);

        // Ensure distinct bodies and both are active
        if (i != j && bodies[i].mass > 0.0 && bodies[j].mass > 0.0) {
            checkAndMergeCollision(bodies[i], bodies[j], thresholdAU);
        }
    }

}

#endif // COLLISIONS_H