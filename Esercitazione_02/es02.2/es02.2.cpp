#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include "random.h"

using namespace std;

// Function to compute the statistical error
double sigma(double, double, int);
// Function to compute the squared norm of a vector
double norm2(double*, int);
// Function to perform a 3D random walk and compute distance
double random_step3D(double*, double, int, int);

int main(){
    // Initialization of random number generator
    Random rnd;
    int M = 100000; // Number of walkers
    int N = 100;    // Number of blocks
    int m = M / N;  // Number of walkers per block

    // Save parameters to a file
    ofstream param;
    param.open("parameters.dat");
    param << M << " " << N << endl;
    param.close();

    double a = 1.; // Step magnitude

    // Open stream to output file
    ofstream output;
    output.open("distances.out");

    // Cycle on the array of lengths of the walk
    for (int length = 0; length < 100; length++){
        // Arrays to store average distances and their squares
        double* average_distance = new double[N];
        double* average_distance2 = new double[N];
        for (int i = 0; i < N; i++){
            average_distance[i] = 0;        
            average_distance2[i] = 0;

            for (int j = 0; j < m; j++){
                // Performs a 3D random walk
                double distance;
                double x[] = {0, 0, 0}; // Initial position at the origin
                for (int l = 0; l < length + 1; l++){
                    double directions[] = {0, 1, 2}; // Array of directions
                    double prob_dir[] = {1/3, 1/3, 1/3}; // Array of probabilities of directions (equally likely)
                    double choices[] = {-1, 1};
                    double prob_choice[] = {0.5, 0.5};
                    
                    // Generate a random step in a random direction
                    double dir = rnd.Sample(directions, prob_dir, 3);
                    double choice = rnd.Sample(choices, prob_choice, 2);
                    
                    // Register the distance travelled (useful only for the last step)
                    distance = random_step3D(x, a, dir, choice);
                }
                average_distance[i] += distance; 
            }
            
            // Take the average of |r|^2
            average_distance[i] /= m;
            // Square root of |r|^2
            average_distance[i] = sqrt(average_distance[i]);
            average_distance2[i] = average_distance[i] * average_distance[i];
        }
        cout << length << "-th chunk of RW computed" << endl << endl;
        
        // Compute last cumulative average with error
        double cum_average = 0;
        double cum_average2 = 0;
        double cum_error = 0;
        for (int i = 0; i < N; i++){
            cum_average += average_distance[i];
            cum_average2 += average_distance2[i];
        }
        cum_average /= N;
        cum_average2 /= N;
        cum_error = sigma(cum_average, cum_average2, N - 1);
        output << length + 1 << " " << cum_average << " " << cum_error << " " << endl;

        delete [] average_distance;
        delete [] average_distance2; 
    }
    output.close();
    return 0;
}

// Function to compute the statistical error
double sigma(double av, double av2, int dim){
    if((dim) == 0){
        return 0;
    }
    else{
        return sqrt((av2 - pow(av, 2)) / (dim));
    }
}

// Function to compute the squared norm of a vector
double norm2(double* x, int dim){
    double distance = 0;
    for (int i = 0; i < dim; i++){
        distance += x[i] * x[i];
    }
    return distance;
}

// Function to perform a 3D random walk and compute distance
double random_step3D(double* position, double step, int direction, int choice){
    // Update position
    position[direction] += choice * step;
    // Compute the distance from the origin
    double distance = norm2(position, 3);
    return distance;
}
