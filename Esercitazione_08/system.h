#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <iomanip>
#include <vector>
#include "random.h"

using namespace std;

#ifndef __energy_measure__
#define __energy_measure__

struct energy_meas{
    double _energy;
    double _error;
};

#endif

#ifndef __System__
#define __System__

class System
{
private:
    double _position;     // Current position in the simulation
    double _mu, _sigma;   // Parameters of the trial wave function
    double _delta;        // Step size for Metropolis moves
    int _nblocks;         // Number of blocks for block averaging
    int _nsteps;          // Number of simulation steps in each block
    int _nattempts;       // Number of attempted Metropolis moves
    int _naccepted;       // Number of accepted Metropolis moves
    double _block_av;     // Block average for energy
    double _global_av;    // Global average of energy
    double _global_av2;   // Global average of the square of the energy
    double _average;      // Current average of measurements
    double _measurement;  // Current measurement of energy
    int _nbins;           // Number of bins for histogram of sampled psi^2
    double* _sampled_psi2;// Array to store histogram of sampled psi^2 values
                          // First 100 bins are for values between (-2, 0) and last 100 for (0, 2)
    double _binsize;      // Size of each bin in the histogram
    Random _rnd;          // Random number generator instance
public:
    System();             // Constructor
    ~System();            // Destructor

    // Setters for mu, sigma, and delta
    void set_mu(double);
    void set_sigma(double);
    void set_delta(double);

    // Getters for mu, sigma, number of blocks, and number of steps per block
    double get_mu();
    double get_sigma();
    int get_nbl();
    int get_nsteps();

    double wavefunction(double); // Function to evaluate the trial wave function at a given position    
    double wavefunc2der(double); // Function to evaluate the second derivative of the trial wave function at a given position
    double potential(double); // Function to evaluate the potential energy at a given position

    void step(); // Function to perform a single Metropolis step    
    bool metro(double); // Function to perform the Metropolis algorithm to sample the distribution
    energy_meas finalize(); // Function to finalize the simulation and return the last progressive average of the Hamiltonian with error
    void block_reset(int blk); // Function to reset block averages at the beginning of a new block
    void measure(); // Function to perform measurements of energy
    void average_wo_saving(); // Function to compute the average without saving intermediate values
    void average(int blk); // Function to compute and save the average for the current block
    double error(double acc, double acc2, int blk); // Function to compute the statistical error
};

#endif
