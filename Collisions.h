#ifndef COLLISIONS_H
#define COLLISIONS_H

#include <vector>
#include <cmath>
#include "Types.h"
#include "Constants.h"

// Merge b into a when distance(a,b) < thresholdAU.
// Conserves mass and linear momentum. Positions in AU; velocities in m/s.
// TODO Complete this function ...
inline void checkAndMergeCollision(body& a, body& b, double thresholdAU = COLLISION_THRESHOLD_AU) {



}

// Quadratic pass: merges any pairs closer than threshold; skips retired bodies.
template <typename Bodies>
inline void mergeAllCollisions(Bodies& bodies, double thresholdAU = COLLISION_THRESHOLD_AU) {
    const int n = static_cast<int>(bodies.size());
   /*
    for (int i = 0; i < n; ++i) {
        if (bodies[i].mass <= 0.0) continue;
        for (int j = i + 1; j < n; ++j) {
            if (bodies[j].mass <= 0.0) continue;
            checkAndMergeCollision(bodies[i], bodies[j], thresholdAU);
        }
    }
    */


    // TODO make this function faster by using random number generator ...


}

#endif // COLLISIONS_H