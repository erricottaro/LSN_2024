#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include "random.h"

using namespace std;

// Function to compute statistical error
double sigma(double, double, int);
// Function to compute squared norm of a vector                      
double norm2(double*, int); 
// Function to perform a random step on a sphere                                 
double random_step_sphere(double*, double, double, double);  

int main() {
    // Initialization
    Random rnd;
    int M = 100000;  // Total number of walkers
    int N = 100;     // Number of blocks
    int m = M / N;   // Number of walkers per block
    
    double a = 1.;   // Step magnitude
    
    // Open output file
    ofstream output;
    output.open("distances_sphere.out");

    // Iterate over different lengths of the random walks
    for (int length = 0; length < 100; length++) {
        // Arrays to store distances and their squares for each block
        double* average_distance = new double[N];
        double* average_distance2 = new double[N];
        
        // Initialize arrays
        for (int i = 0; i < N; i++) {
            average_distance[i] = 0;
            average_distance2[i] = 0;
            
            // Iterate over walkers in the block
            for (int j = 0; j < m; j++) {
                // Initialize position
                double distance;
                double x[] = {0, 0, 0};  // Initial position at origin
                
                // Perform random walk
                for (int l = 0; l < length + 1; l++) {
                    // Generate random solid angles
                    double theta = rnd.Sine();
                    double phi = rnd.Rannyu(0, 2 * M_PI);
                    
                    // Perform step and compute distance
                    distance = random_step_sphere(x, a, theta, phi); 
                }
                
                // Accumulate distance
                average_distance[i] += distance; 
            }
            
            // Compute average distance and its square
            average_distance[i] /= m;
            average_distance[i] = sqrt(average_distance[i]);
            average_distance2[i] = average_distance[i] * average_distance[i];
        }
        
        // Compute cumulative average with error
        double cum_average = 0;
        double cum_average2 = 0;
        double cum_error = 0;
        for (int i = 0; i < N; i++) {
            cum_average += average_distance[i];
            cum_average2 += average_distance2[i];
        }
        cum_average /= N;
        cum_average2 /= N;
        cum_error = sigma(cum_average, cum_average2, N - 1);
        
        // Write results to output file
        output << length + 1 << " " << cum_average << " " << cum_error << " " << endl;

        // Clean up memory
        delete[] average_distance;
        delete[] average_distance2; 
    }
    
    // Close output file
    output.close();

    return 0;
}

// Function to compute statistical error
double sigma(double av, double av2, int dim) {
    if ((dim) == 0) {
        return 0;
    } else {
        return sqrt((av2 - pow(av, 2)) / (dim));
    }
}

// Function to compute squared norm of a vector
double norm2(double* x, int dim) {
    double distance = 0;
    for (int i = 0; i < dim; i++) {
        distance += x[i] * x[i];
    }
    return distance;
}

// Function to perform a random step on a sphere
double random_step_sphere(double* position, double step, double theta, double phi) {
    // Update position
    position[0] += step * sin(theta) * cos(phi);
    position[1] += step * sin(theta) * sin(phi);
    position[2] += step * cos(theta);
    
    // Compute squared distance from the origin
    double distance = norm2(position, 3);
    return distance;
}
