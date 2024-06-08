#include <iostream>  // Input-output stream library
#include <fstream>   // File stream library
#include <cmath>     // Math functions
#include <cstdlib>   // Standard library definitions
#include "random.h"  // Custom random number generator class

using namespace std;

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

    // Parameters of the simulation
    int N = 10000;      // Number of throws
    int M = 100;        // Number of subintervals (bins)
    int K = 100;        // Number of chi-squared values to compute
    double probability = 1 / double(M);  // Probability of each bin
    double average = N / double(M);      // Average expected count per bin

    // Open an output file for storing chi-squared values
    ofstream chi;
    chi.open("chi.out");

    // Loop over K iterations to compute chi-squared values
    for (int l = 0; l < K; l++){
        // Array of counters to store the number of elements in each bin
        int* count = new int[M];
        for (int i = 0; i < M; i++) {count[i] = 0;} // Initialize the counters to zero

        // Generate random numbers and count them into bins
        for (int i = 0; i < N; i++){
            double random = rnd.Rannyu();  // Generate a pseudo-random number

            // Iterate over the number of intervals to assign the random number to a bin
            for (int j = 0; j < M; j++){
                if (random >= j * probability && random < (j + 1) * probability){
                    count[j]++;  // Increment the count for the corresponding bin
                    break;       // No need to check further as intervals do not overlap
                }
            }
        }

        // Compute the chi-squared value
        double chi_squared = 0;
        for (int i = 0; i < M; i++){
            chi_squared += (pow(count[i] - average, 2) / average);
        }
        chi << chi_squared << endl;  // Write the chi-squared value to the output file

        delete [] count;  // Free memory allocated for the count array
    }

    chi.close();  // Close the output file
    return 0;
}
