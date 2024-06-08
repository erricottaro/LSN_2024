#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <iomanip>
#include "random.h"
#include "system.h"

using namespace std;

#ifndef __annealing__
#define __annealing__

class annealing
{
private:
    Random _rnd;        // Random number generator instance
    System _SYS;        // Instance of the System class (quantum system to be optimized)
    double _delta_mu;   // Step size for mu parameter adjustments
    double _delta_sigma;// Step size for sigma parameter adjustments
    int _nsteps;        // Number of simulation steps for the annealing process
    energy_meas _old_energy; // Energy measurement of the system with old parameters
    energy_meas _new_energy; // Energy measurement of the system with new parameters
    double _measurement;// Current measurement of energy
    double _average;    // Average energy measurement
    int _naccepted;     // Number of accepted parameter changes
    int _nattempts;     // Number of attempted parameter changes
    double _temp;       // Current temperature for the annealing process
    double _beta;       // Inverse temperature (1/temp)

public:
    annealing();        // Constructor
    ~annealing();       // Destructor

    // Getters for number of steps, beta (inverse temperature), mu, and sigma
    int get_nsteps();
    double get_beta();
    double get_mu();
    double get_sigma();

    void step(); // Function to perform a single annealing step
    bool metro(); // Metropolis algorithm to accept or reject new parameters

    void finalize(); // Function to finalize the annealing process
    energy_meas get_energy(); // Getter for the current energy measurement
    void measure(); // Function to perform measurements of energy
    double error(double acc, double acc2, int blk); // Function to compute the statistical error
};

#endif