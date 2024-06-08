#include <iostream>   // Input-output stream library
#include <fstream>    // File stream library
#include <string>     // String manipulation library
#include "random.h"   // Custom random number generator class

using namespace std;

int main(){
    // Instance of the random number generator class
    Random rnd;
    // Number of realizations of each mean variable
    int realizations = 10000;
    // Array with the dimensions of the samples to take the average of
    int dims[] = {1, 2, 10, 100};
    // Dimension of the previous array
    int n_iterations = sizeof(dims) / sizeof(dims[0]);

    // Open a file to store parameters
    ofstream param;
    param.open("parameters.dat");
    for (int i = 0; i < n_iterations; i++){
        param << dims[i] << endl;  // Write sample dimensions to the file
    }

    // Iteration over the dimensions of samples
    for (int i = 0; i < n_iterations; i++){
        // Declare and open streams to output files for different distributions
        ofstream out_unif;
        ofstream out_exp;
        ofstream out_lorentz;
        out_unif.open("unif" + to_string(dims[i]) + ".out");
        out_exp.open("exp" + to_string(dims[i]) + ".out");
        out_lorentz.open("lorentz" + to_string(dims[i]) + ".out");

        // Iteration over realizations
        for (int k = 0; k < realizations; k++){
            // Initialize variables to store averages for different distributions
            double average_unif = 0;
            double average_exp = 0;
            double average_lorentz = 0;

            // Iteration over sample size
            for (int j = 0; j < dims[i]; j++){
                // Generate random numbers from different distributions and accumulate the sum
                // Uniformly distributed
                double random = rnd.Rannyu();
                average_unif += random;
                // Exponentially distributed
                random = rnd.Exponential(1.);
                average_exp += random;
                // Lorentzian distribution
                random = rnd.Lorentz(0., 1.);
                average_lorentz += random;
            }
            
            // Compute the average for each distribution
            average_unif /= dims[i];
            average_exp /= dims[i];
            average_lorentz /= dims[i];

            // Write the averages to the respective output files
            out_unif << average_unif << endl;
            out_exp << average_exp << endl;
            out_lorentz << average_lorentz << endl;
        }
        // Close the output files for each distribution
        out_unif.close();
        out_exp.close();
        out_lorentz.close();
    }
    return 0;
}
