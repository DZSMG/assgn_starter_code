#ifndef OCTANT_H
#define OCTANT_H

#include "Types.h"

class Octant {
    vec3   mid;    // center (AU)
    double length; // full box width (AU)
public:
    Octant(double x, double y, double z, double l) : mid{x, y, z}, length{l} {}
    [[nodiscard]] Octant mUNW() const; // +y +z -x
    [[nodiscard]] Octant mUNE() const; // +y +z +x
    [[nodiscard]] Octant mUSW() const; // -y +z -x
    [[nodiscard]] Octant mUSE() const; // -y +z +x
    [[nodiscard]] Octant mDNW() const; // +y -z -x
    [[nodiscard]] Octant mDNE() const; // +y -z +x
    [[nodiscard]] Octant mDSW() const; // -y -z -x
    [[nodiscard]] Octant mDSE() const; // -y -z +x
    [[nodiscard]] bool   contains(const vec3& p) const;
    [[nodiscard]] double getLength() const { return length; }
    [[nodiscard]] const vec3& center() const { return mid; }
};

#endif