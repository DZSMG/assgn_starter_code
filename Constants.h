#ifndef CONSTANTS_H
#define CONSTANTS_H

// ---------------- Image output ----------------
#define WIDTH 800
#define HEIGHT 800
// ---------------- Simulation setup -------------
#define NUM_BODIES  (10*1024)   // default number of bodies in disk
#define SYSTEM_SIZE_AU  10.0    // baseline view radius in AU
#define DISK_THICKNESS_AU 0.08    // default disk thickness in AU
#define INNER_BOUND_AU 0.3     // default disk inner radius in AU

// ---------------- Physics -------------------
#define PI 3.14159265358979323846
#define AU 1.496e11           // meters per AU
#define G 6.67408e-11        // m^3 kg^-1 s^-2
#define SOLAR_MASS 2.0e30             // kg
#define DISK_MASS_FRAC 0.2             // default disk fraction
#define SOFTENING_M (0.02 * AU)          // gravitational softening (m)
#define BH_THETA 0.5                // Barnes–Hut opening angle

// ---------------- Integrator -----------------
#define DT (1000.0 * 1024.0) // seconds per step (smaller for stability)
#define STEPS 100 // number of frames, lower = shorter simulation

// ---------------- Rendering ------------------
#define VIEW_HALF_AU SYSTEM_SIZE_AU // +/- AU span mapped to image
#define RENDER_ZOOM 1.0 // 1.0=default; <1 zoom in, >1 out

// ---------------- Collisions -----------------
#define COLLISION_THRESHOLD_AU 0.00001 // 1.0e-8 // very small values = collisions are off

// ---------------- Drift / centering ----------
#define NO_FIX 0
#define ANCHOR_FIRST_STAR 1
#define REMOVE_BARYCENTRE_VELOCITY 2
#define REMOVE_BARYCENTRE_VELOCITY_EACH_STEP 3
#define DRIFT_FIX REMOVE_BARYCENTRE_VELOCITY_EACH_STEP

#endif
