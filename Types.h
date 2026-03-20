#ifndef TYPES_H
#define TYPES_H

#include <cmath>

struct vec3 {
    double x{0}, y{0}, z{0};
};

inline vec3 operator+(const vec3& a, const vec3& b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
inline vec3 operator-(const vec3& a, const vec3& b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
inline vec3 operator*(const vec3& a, double s)     { return {a.x * s, a.y * s, a.z * s}; }
inline vec3 operator/(const vec3& a, double s)     { return {a.x / s, a.y / s, a.z / s}; }
inline double dot(const vec3& a, const vec3& b)    { return a.x*b.x + a.y*b.y + a.z*b.z; }
inline double norm(const vec3& a)                  { return std::sqrt(dot(a, a)); }

struct body {
    vec3 position{0,0,0}; // AU
    vec3 velocity{0,0,0}; // m/s
    vec3 accel{0,0,0};    // m/s^2
    double mass{0.0}; // kg
};

// Flexible per-galaxy specification
struct StarSpec {
    double mass{};                       // central bulge/black-hole surrogate (kg)
    vec3   position{0,0,0};      // AU
    vec3   velocity{0,0,0};      // m/s
    double diskInnerAU{};               // AU
    double diskOuterAU{};               // AU
    double diskThicknessAU{};           // AU
    int    numBodies{};                 // number of disk bodies
    double diskMassFrac{};              // fraction of central mass in the disk
};

// Kinematic parameters for disk sampling
struct DiskKinematics {
    double RdAU = 2.0;      // exponential scale length (AU)
    double sigmaR = 60.0;   // m/s
    double sigmaPhi = 60.0; // m/s
    double sigmaZ = 30.0;   // m/s
    bool   disableCollisions = true; // galaxies are effectively collision-less
};
#endif