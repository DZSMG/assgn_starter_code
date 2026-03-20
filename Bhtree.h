#ifndef BHTREE_H
#define BHTREE_H

#include "Types.h"
#include "Octant.h"
#include "Constants.h"

// Barnes–Hut tree node for 3D n-body simulation
class Bhtree {
    body   myBod{}; // aggregate centre of mass and total mass for this node
    Octant octy;
    Bhtree* UNW{nullptr}; Bhtree* UNE{nullptr}; Bhtree* USW{nullptr}; Bhtree* USE{nullptr};
    Bhtree* DNW{nullptr}; Bhtree* DNE{nullptr}; Bhtree* DSW{nullptr}; Bhtree* DSE{nullptr};
    // Insert a body into the appropriate child octant.  Assumes this node's
    // octant already contains the body's position.  This helper will
    // allocate the child node on demand and forward the insertion.
    void insertIntoChild(body* route);
public:
    explicit Bhtree(Octant&& oct);
    ~Bhtree();
    [[nodiscard]] const Octant& octant() const { return octy; }
    [[nodiscard]] bool isExternal() const;
    // Insert a body, updating aggregate mass and centre of mass
    void insert(body* insertBod);
    // Compute force on a body by traversing the tree
    static void singleInteract(body* target, const body* other, bool singlePart);
    void interactInTree(body* bod) const;
};

#endif // BHTREE_H