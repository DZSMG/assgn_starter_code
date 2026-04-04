// Task 4: Created Simulation.h to optimize main.cpp

#ifndef ASSGN_SIMULATION_H
#define ASSGN_SIMULATION_H
#pragma once
#include <vector>
#include <random>
#include "Types.h"
#include "Octant.h"
#include "Bhtree.h"
#include "Constants.h"
#include "Parser.h"

// Moved the entire Body struct from main.cpp
struct Body {
    body *ptr{nullptr};  // pointer to array of body
    size_t n{0};
    size_t cap{0};

    Body() = default;  // default constructor

    // keep simple: no accidental copies
    Body(const Body &o) = delete;
    Body &operator=(const Body &o) = delete;

    Body(Body &&o) noexcept // move-only type
        : ptr(o.ptr), n(o.n), cap(o.cap) {
        o.ptr = nullptr;
        o.n = o.cap = 0;
    }

    void reserve(size_t newCap) {
        if (newCap <= cap) return;
        body *np = new body[newCap];
        if (ptr) {
            // copy existing elements
            for (size_t i = 0; i < n; ++i) np[i] = ptr[i];
        }
        delete[] ptr;
        ptr = np;
        cap = newCap;
    }


    // This is working perfectly now! Don't mess with it!
    Body &operator=(Body &&o) noexcept {
        if (this != &o) {
            delete[] ptr;
            ptr = o.ptr;
            n = o.n;
            cap = o.cap;
            o.ptr = nullptr;
            o.n = o.cap = 0;
        }
        return *this;
    }

    // Everything else below also works perfectly!
    body *begin() { return ptr; }
    body *end() { return ptr + n; } // always ptr+n (works even if ptr==nullptr && n==0)
    const body *begin() const { return ptr; }
    const body *end() const { return ptr + n; }

    void resize(size_t newSize) {
        if (newSize > cap) reserve(newSize);
        n = newSize;
    }

    void clear() { n = 0; }

    bool empty() const { return n == 0; }
    size_t size() const { return n; }

    body &operator[](size_t i) { return ptr[i]; }
    const body &operator[](size_t i) const { return ptr[i]; }

    void push_back(const body &b) {
        if (n + 1 > cap) reserve(cap ? cap * 2 : 8);
        ptr[n++] = b;
    }
};

// Declare the Physics and Initialization functions
StarSpec makeMilkyWaySpec(vec3 pos, vec3 vel, int numBodies, double bulgeMass, double diskInnerAU, double diskOuterAU, double diskThicknessAU, double diskMassFrac);

void initialiseStarsAndBodies_Exponential(Body &bb, const std::vector<StarSpec> &s0, const DiskKinematics &kin);

Octant makeDynamicRootOctant(const Body &bb, double marginAU);

void simulateStep(Body &bb);

#endif //ASSGN_SIMULATION_H
