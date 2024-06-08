#include <iostream>  // Input-output stream library
#include <fstream>   // File stream library
#include <cmath>     // Math functions
#include <cstdlib>   // Standard library definitions
#include "random.h"  // Custom random number generator class

using namespace std;

// Function prototype for calculating the standard deviation
double sigma(double, double, int);

int main(){
    // Instance of the class that produces the random numbers
    Random rnd;
    // Set up the generator
    int seed[4];
    int p1, p2;
    
    // Read prime numbers from file "Primes" for random seed initialization
    ifstream Primes("Primes");
    if (Primes.is_open()){
       Primes >> p1 >> p2 ;
    } else cerr << "PROBLEM: Unable to open Primes" << endl;
    Primes.close();
    
    // Read initial random seed from file "seed.in"
    ifstream input("seed.in");
    string property;
    if (input.is_open()){
       while ( !input.eof() ){
          input >> property;
          if( property == "RANDOMSEED" ){
             input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
             rnd.SetRandom(seed,p1,p2);
          }
       }
       input.close();
    } else cerr << "PROBLEM: Unable to open seed.in" << endl;

    // Total number of throws
    int M = 100000;
    // Number of blocks
    int N = 100;
    // Save parameters of the simulation to an external file
    ofstream param;
    param.open("parameters.dat");
    param << M << " " << N << endl;
    param.close();
    // Number of throws per block
    int L = M/N;
    // Accumulator variables for the expected value
    double* average = new double[N];
    double* average2 = new double[N];
    // Accumulator variables for the variance
    double* variance = new double[N];
    double* variance2 = new double[N];

    // Cycle over the number of blocks
    for (int i = 0; i < N; i++){
        // Accumulator for the expected value
        double sum_block = 0;
        // Accumulator for the variance
        double var_block = 0;
        for (int j = 0; j < L; j++){
            // Generate a random number and add it to the total
            double rand_value = rnd.Rannyu();
            sum_block += rand_value;
            var_block += pow(rand_value - 0.5, 2);
        }
        average[i] = sum_block / L;
        average2[i] = pow(average[i], 2);

        variance[i] = var_block / L;
        variance2[i] = pow(variance[i], 2);
    }
    // <r> related variables
    ofstream output1;
    output1.open("average.out");
    // sigma^2 related variables
    ofstream output2;
    output2.open("variance.out");
    
    for (int i = 0; i < N; i++){
        // Initialize accumulators
        // <r>
        double cumul_average = 0;
        double cumul_average2 = 0;
        double cumul_error;
        // sigma^2
        double cumul_variance = 0;
        double cumul_variance2 = 0;
        double cumul_errorvar;
        for (int j = 0; j < i + 1; j++){  // All the contributions up to the i-th one are summed
            cumul_average += average[j];
            cumul_average2 += average2[j];

            cumul_variance += variance[j];
            cumul_variance2 += variance2[j];
        }
        // <r>
        cumul_average /= (i + 1);
        cumul_average2 /= (i + 1);
        cumul_error = sigma(cumul_average, cumul_average2, i);
        output1 << cumul_average << " " << cumul_error << endl;
        // sigma^2
        cumul_variance /= (i + 1);
        cumul_variance2 /= (i + 1);
        cumul_errorvar = sigma(cumul_variance, cumul_variance2, i);
        output2 << cumul_variance << " " << cumul_errorvar << endl;
    }

    output1.close();
    delete [] average;
    delete [] average2;
    output2.close();
    delete [] variance;
    delete [] variance2;
    
    return 0;
}

// Function for calculating the standard deviation
double sigma(double av, double av2, int dim){
    if(dim == 0){
        return 0;
    } else {
        return sqrt((av2 - pow(av, 2)) / dim);
    }
}
