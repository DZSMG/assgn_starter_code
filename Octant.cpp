#include "Octant.h"

static Octant makeChild(const vec3& mid, double length, double sx, double sy, double sz) {
    const double quarter = length / 4.0;
    return { mid.x + sx * quarter, mid.y + sy * quarter, mid.z + sz * quarter, length / 2.0 };
}

bool Octant::contains(const vec3& p) const {
    const double h = length / 2.0;
    return (p.x <= mid.x + h && p.x >= mid.x - h) &&
           (p.y <= mid.y + h && p.y >= mid.y - h) &&
           (p.z <= mid.z + h && p.z >= mid.z - h);
}
Octant Octant::mUNW() const { return makeChild(mid, length, -1, +1, +1); }
Octant Octant::mUNE() const { return makeChild(mid, length, +1, +1, +1); }
Octant Octant::mUSW() const { return makeChild(mid, length, -1, -1, +1); }
Octant Octant::mUSE() const { return makeChild(mid, length, +1, -1, +1); }
Octant Octant::mDNW() const { return makeChild(mid, length, -1, +1, -1); }
Octant Octant::mDNE() const { return makeChild(mid, length, +1, +1, -1); }
Octant Octant::mDSW() const { return makeChild(mid, length, -1, -1, -1); }
Octant Octant::mDSE() const { return makeChild(mid, length, +1, -1, -1); }
