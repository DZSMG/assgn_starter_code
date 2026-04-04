#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <vector>

#include "Constants.h"
#include "Types.h"
#include "Simulation.h"
#include "Renderer.h"
#include "Parser.h"


/*
 * Task 4: Cleaned up the program to strictly act as the controller.
 * Implement function encapsulation to prevent misuses and namespace pollution.
 */
int main(int argc, char* argv[]) {
    SimConfig config;

    // Run Parser - exit if invalid parameters or help menu called
    if (!parseArguments(argc, argv, config)) {
        std::cout.flush();
        return 0;
    }

    // this will execute the deleteImgs PowerShell script
    system("powershell -ExecutionPolicy Bypass -File \"../deleteImgs.ps1\"");
    // system("bash deleteImgs.bash");

    // Build a single Milky Way-like galaxy
    std::vector<StarSpec> s0;
    s0.push_back(makeMilkyWaySpec({0, 0, 0}, {0, 0, 0},
                                  config.numBodies,
                                  1.2 * SOLAR_MASS,
                                  0.3,
                                  config.systemSize,
                                  0.2,
                                  0.25));

    DiskKinematics kin; // default dispersions and Rd
    Body bb;
    initialiseStarsAndBodies_Exponential(bb, s0, kin);

    std::cout << "Running Barnes-Hut n-body (" << bb.size()
            << " bb, " << STEPS << " steps)..." << std::endl;

    for (int s = 0; s < STEPS; ++s) {
        simulateStep(bb);
        std::ostringstream fname;
        fname << "images/frame_" << std::setw(3) << std::setfill('0') << s << ".ppm";
        renderSnapshot(bb, fname.str(), RENDER_ZOOM);
        if ((s + 1) % 50 == 0) {
            std::cout << " stp " << (s + 1) << "/" << STEPS << "\r" << std::flush;
        }
    }

    std::cout << "\nDone. Wrote " << STEPS << " frames as PPM images." << std::endl;

    // this will execute the createVideo PowerShell script
    system("powershell -ExecutionPolicy Bypass -File \"../createVideo.ps1\"");
    // system("bash createVideo.bash");
    return 0;
}
