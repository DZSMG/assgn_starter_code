// Created for command line parser

#include "Parser.h"
#include <iostream>
#include <string>
#include <stdexcept>

static void printHelp() {
    std::cout << "========================================================\n";
    std::cout << "          N-Body Simulation Help Menu                   \n";
    std::cout << "========================================================\n";
    std::cout << "Usage: ./assgn.exe [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  --help, -h          Show this help message\n";
    std::cout << "  --bodies <num>      Set total number of bodies (integer)\n";
    std::cout << "  --size <val>        Set system size in AU (double)\n";
    std::cout << "  --steps <num>       Set total simulation steps (integer)\n";
    std::cout << "========================================================\n";
}

bool parseArguments(int argc, char* argv[], SimConfig& config) {
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "--help" || arg == "-h") {
            printHelp();
            return false;
        }

        try {
            if (arg == "--size" && i + 1 < argc) {
                config.systemSize = std::stod(argv[++i]);
                config.viewHalf = config.systemSize; // Sync rendering camera
                if (config.systemSize <= 0) throw std::invalid_argument("Size must be strictly positive.");
            }
            else if (arg == "--bodies" && i + 1 < argc) {
                config.numBodies = std::stoull(argv[++i]);
                if (config.numBodies == 0) throw std::invalid_argument("Bodies cannot be zero.");
            }
            else if (arg == "--steps" && i + 1 < argc) {
                config.steps = std::stoi(argv[++i]);
                if (config.steps <= 0) throw std::invalid_argument("Steps must be greater than zero.");
            }
            else {
                std::cerr << "[Error] Unknown argument: " << arg << "\n";
                return false;
            }
        }
        catch (const std::exception& e) {
            std::cerr << "[Error] Invalid value for " << arg << ": " << e.what() << "\n";
            std::cerr << "Simulation safely aborted to prevent memory corruption.\n";
            return false;
        }
    }
    return true;
}