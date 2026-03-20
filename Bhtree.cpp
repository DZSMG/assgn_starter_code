#include "Bhtree.h"
#include <cmath>
#include <utility>

Bhtree::Bhtree(Octant &&oct) : octy(std::move(oct)) {
}

Bhtree::~Bhtree() {
    delete UNW;
    delete UNE;
    delete USW;
    delete USE;
    delete DNW;
    delete DNE;
    delete DSW;
    delete DSE;
}

bool Bhtree::isExternal() const {
    return !UNW && !UNE && !USW && !USE && !DNW && !DNE && !DSW && !DSE;
}

void Bhtree::insert(body* insertBod) {
    // Discard bodies outside this octant
    if (!octy.contains(insertBod->position)) return;

    // Empty node: take ownership of the first body
    if (myBod.mass == 0.0) {
        myBod = *insertBod;
        return;
    }

    // Leaf node (external) with one body: subdivide and insert both the
    // existing occupant and the incoming body into children, then update
    // the aggregate mass/centre of mass once.
    if (isExternal()) {
        // Preserve the existing occupant's properties before they are overwritten
        body old = myBod;
        // Compute new aggregate mass and centre of mass (in AU)
        const double totalMass = old.mass + insertBod->mass;
        myBod.position.x = (old.position.x * old.mass + insertBod->position.x * insertBod->mass) / totalMass;
        myBod.position.y = (old.position.y * old.mass + insertBod->position.y * insertBod->mass) / totalMass;
        myBod.position.z = (old.position.z * old.mass + insertBod->position.z * insertBod->mass) / totalMass;
        myBod.mass = totalMass;
        // Insert the saved occupant into the appropriate child
        insertIntoChild(&old);
        // Insert the incoming body into the appropriate child
        insertIntoChild(insertBod);
        return;
    }

    // Internal node: update aggregate mass and centre of mass, then pass
    // the body down to the appropriate child
    const double totalMass = myBod.mass + insertBod->mass;
    myBod.position.x = (myBod.position.x * myBod.mass + insertBod->position.x * insertBod->mass) / totalMass;
    myBod.position.y = (myBod.position.y * myBod.mass + insertBod->position.y * insertBod->mass) / totalMass;
    myBod.position.z = (myBod.position.z * myBod.mass + insertBod->position.z * insertBod->mass) / totalMass;
    myBod.mass = totalMass;
    // Forward insertion to a child octant
    insertIntoChild(insertBod);
}

// Choose the appropriate child octant for the given body and insert it
void Bhtree::insertIntoChild(body* route) {
    const vec3& routePos = route->position;
    // Determine which of the eight sub-octants contains the body.
    // Note: Octant::contains uses inclusive bounds, so on boundaries the first matching branch will be selected.
    if (Octant unw = octy.mUNW(); unw.contains(routePos)) {
        if (!UNW) UNW = new Bhtree(std::move(unw));
        UNW->insert(route);
    } else if (Octant une = octy.mUNE(); une.contains(routePos)) {
        if (!UNE) UNE = new Bhtree(std::move(une));
        UNE->insert(route);
    } else if (Octant usw = octy.mUSW(); usw.contains(routePos)) {
        if (!USW) USW = new Bhtree(std::move(usw));
        USW->insert(route);
    } else if (Octant use = octy.mUSE(); use.contains(routePos)) {
        if (!USE) USE = new Bhtree(std::move(use));
        USE->insert(route);
    } else if (Octant dnw = octy.mDNW(); dnw.contains(routePos)) {
        if (!DNW) DNW = new Bhtree(std::move(dnw));
        DNW->insert(route);
    } else if (Octant dne = octy.mDNE(); dne.contains(routePos)) {
        if (!DNE) DNE = new Bhtree(std::move(dne));
        DNE->insert(route);
    } else if (Octant dsw = octy.mDSW(); dsw.contains(routePos)) {
        if (!DSW) DSW = new Bhtree(std::move(dsw));
        DSW->insert(route);
    } else {
        Octant dse = octy.mDSE();
        if (!DSE) DSE = new Bhtree(std::move(dse));
        DSE->insert(route);
    }
}

void Bhtree::singleInteract(body *target, const body *other, bool singlePart) {
    // r = target - other (meters)
    vec3 rM = (target->position - other->position) * AU;
    const double r2 = rM.x * rM.x + rM.y * rM.y + rM.z * rM.z;
    if (singlePart && r2 == 0.0) return; // skip self-force on leaf
    const double denom = std::sqrt(r2 + SOFTENING_M * SOFTENING_M);
    const double invDenom3 = 1.0 / (denom * denom * denom);
    const double factor = -G * other->mass * invDenom3; // s^-2 per meter
    target->accel.x += factor * rM.x;
    target->accel.y += factor * rM.y;
    target->accel.z += factor * rM.z;
}

void Bhtree::interactInTree(body *bod) const {
    if (myBod.mass == 0.0) return;
    if (isExternal()) {
        singleInteract(bod, &myBod, true);
        return;
    }

    // Distance to node aggregate COM (AU, for BH criterion)
    const double dx = myBod.position.x - bod->position.x;
    const double dy = myBod.position.y - bod->position.y;
    const double dz = myBod.position.z - bod->position.z;
    const double d = std::sqrt(dx * dx + dy * dy + dz * dz);

    if (d <= 1e-12) {
        // Avoid degeneracy: if very close, open the node
        if (UNW) UNW->interactInTree(bod);
        if (UNE) UNE->interactInTree(bod);
        if (USW) USW->interactInTree(bod);
        if (USE) USE->interactInTree(bod);
        if (DNW) DNW->interactInTree(bod);
        if (DNE) DNE->interactInTree(bod);
        if (DSW) DSW->interactInTree(bod);
        if (DSE) DSE->interactInTree(bod);
        return;
    }

    const double s = octy.getLength(); // box width (AU)
    if (s / d < BH_THETA) {
        // Treat as pseudo-body at COM
        singleInteract(bod, &myBod, false);
    } else {
        if (UNW) UNW->interactInTree(bod);
        if (UNE) UNE->interactInTree(bod);
        if (USW) USW->interactInTree(bod);
        if (USE) USE->interactInTree(bod);
        if (DNW) DNW->interactInTree(bod);
        if (DNE) DNE->interactInTree(bod);
        if (DSW) DSW->interactInTree(bod);
        if (DSE) DSE->interactInTree(bod);
    }
}
