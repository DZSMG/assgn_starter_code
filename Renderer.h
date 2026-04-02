// Task 4: Created Renderer.h to optimize main.cpp

#ifndef ASSGN_RENDERER_H
#define ASSGN_RENDERER_H
#pragma once

#include <string>
#include <vector>
#include <cstdint>
#include "Simulation.h"

void renderSnapshot(const Body &B, const std::string &file, double zoom = 1.0);


#endif //ASSGN_RENDERER_H
