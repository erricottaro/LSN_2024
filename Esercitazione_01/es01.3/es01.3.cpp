#include <iostream>   // Input-output stream library
#include <fstream>    // File stream library
#include <string>     // String manipulation library
#include <cstdlib>    // Standard library
#include <cmath>      // Math library
#include "random.h"   // Custom random number generator class

using namespace std;

// Function declaration for computing statistical error
double sigma(double, double, int);

int main(){
    // Parameters of the simulation: grid displacement d (unitary), needle length L (L < d)
    double d = 1.;
    double L = 0.5;  // Arbitrary choice (not too small)
    int M = 1000000; // Number of throws
    int N = 100;     // Number of blocks
    int m = M/N;     // Number of throws per block
    double scaling = 2 * L / d; // Scaling constant

    // Save parameters to a file
    ofstream param;
    param.open("parameters.dat");
    param << d << " " << L << " " << M << " " << N << endl;
    param.close();

    Random rnd; // Instance of the random number generator class

    // Arrays to store average frequency per block and its square for computing error
    // Frequency = N_hits/N
    double* pi_estim = new double[N];
    double* pi_estim2 = new double[N];

    // Cycle over the number of blocks
    for (int i = 0; i < N; i++){
        pi_estim[i] = 0;
        // Generate the needle
        for (int j = 0; j < m; j++){
            // x coordinate of the midpoint of the needle
            double x_mid = rnd.Rannyu(0, d);
            // Slope of the needle
            // Generate two random points on the plane to get a random slope
            double x_1 = rnd.Rannyu();
            double y_1 = rnd.Rannyu();
            double x_2 = rnd.Rannyu();
            double y_2 = rnd.Rannyu();
            double slope = (y_2 - y_1) / (x_2 - x_1);
            // Endpoints
            double x_left = x_mid - 0.5 * L / sqrt(1 + slope * slope);
            double x_right = x_mid + 0.5 * L / sqrt(1 + slope * slope);
            // Check if needle crosses the boundaries of the stripe of length d
            if (x_left < 0 || x_right > d){
                pi_estim[i]++;
            }
        }
        // Compute estimation of π from frequency
        pi_estim[i] /= m;
        pi_estim[i] = scaling / pi_estim[i]; // π depends on 1/frequency
        pi_estim2[i] = pow(pi_estim[i], 2);
    }

    // Computation of π and error
    ofstream output;
    output.open("pi_estim.out");
    for (int i = 0; i < N; i++){
        double cumul_av = 0; 
        double cumul_av2 = 0;
        double cumul_error = 0;
        for (int j = 0; j < i + 1; j++){
            cumul_av += pi_estim[j];
            cumul_av2 += pi_estim2[j];
        }
        cumul_av /= (i + 1);
        cumul_av2 /= (i + 1);
        cumul_error = sigma(cumul_av, cumul_av2, i + 1);
        output << cumul_av << " " << cumul_error << endl;  
    }

    delete [] pi_estim;
    delete [] pi_estim2;
    output.close();

    return 0;
}

// Function for computing statistical error
double sigma(double av, double av2, int dim){
    if ((dim - 1) == 0){
        return 0;
    }
    else{
        return sqrt((av2 - pow(av, 2)) / (dim - 1));
    }
}
