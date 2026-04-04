// Created for command line parser

#ifndef ASSGN_PARSER_H
#define ASSGN_PARSER_H

#pragma once
#include <cstddef>
#include "Constants.h"

// The runtime configuration object (defaults pulled from Constants.h)
struct SimConfig {
    size_t numBodies   = NUM_BODIES;
    double systemSize  = SYSTEM_SIZE_AU;
    double viewHalf    = VIEW_HALF_AU;
    double dt          = DT;
    int    steps       = STEPS;
    double renderZoom  = RENDER_ZOOM;
};

// Parses command-line arguments and updates the config object.
// Returns false if the program should exit (e.g., error occurred or --help requested).
bool parseArguments(int argc, char* argv[], SimConfig& config);

#endif //ASSGN_PARSER_H
