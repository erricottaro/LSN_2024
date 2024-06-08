#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include "random.h" // Include the custom random number generator header file

using namespace std;

// Function to define the integrand
double integrand(double);
// Function to define the distribution for importance sampling
double distro(double);
// Function to compute the statistical error
double sigma(double, double, int);

int main(){
    // Initialization of the random number generator
    Random rnd;
    int M = 100000; // Total number of throws
    int N = 100;    // Number of blocks
    int m = M / N;  // Number of throws per block

    // Save parameters to a file
    ofstream param;
    param.open("parameters.dat");
    param << M << " " << N << endl;
    param.close();

    // Integral computed with uniform distribution
    // Arrays to store integral values and their squares
    double* integral = new double[N];
    double* integral2 = new double[N];
    for (int i = 0; i < N; i++){
        integral[i] = 0;
        integral2[i] = 0;

        // Monte Carlo integration
        for (int j = 0; j < m; j++){
            // Evaluate the integrand at a random point (uniform distribution)
            integral[i] += integrand(rnd.Rannyu());
        }
        // Compute the average integral value and its square
        integral[i] /= m;
        integral2[i] = pow(integral[i], 2);
    }

    // Computation of the integral and error
    ofstream output;
    output.open("integ_unif.out");
    for (int i = 0; i < N; i++){
        double cumul_av = 0;
        double cumul_av2 = 0;
        double cumul_error = 0;
        for (int j = 0; j < i + 1; j++){
            cumul_av += integral[j];
            cumul_av2 += integral2[j];
        }
        cumul_av /= (i + 1);
        cumul_av2 /= (i + 1);
        cumul_error = sigma(cumul_av, cumul_av2, i);
        output << cumul_av << " " << cumul_error << endl;
    }
    output.close();
    delete [] integral;
    delete [] integral2;

    // Integral computed with importance sampling
    // Arrays to store integral values and their squares
    double* integral_imp = new double[N];
    double* integral_imp2 = new double[N];
    for (int i = 0; i < N; i++){
        integral_imp[i] = 0;
        integral_imp2[i] = 0;

        // Monte Carlo integration with importance sampling
        for (int j = 0; j < m; j++){
            double x = rnd.Linear();
            // Evaluate the integrand at a random point (importance sampling)
            integral_imp[i] += integrand(x) / distro(x);
        }
        // Compute the average integral value and its square
        integral_imp[i] /= m;
        integral_imp2[i] = pow(integral_imp[i], 2);
    }

    // Computation of the integral and error with importance sampling
    ofstream output_imp;
    output_imp.open("integ_importance.out");
    for (int i = 0; i < N; i++){
        double cumul_av = 0;
        double cumul_av2 = 0;
        double cumul_error = 0;
        for (int j = 0; j < i + 1; j++){
            cumul_av += integral_imp[j];
            cumul_av2 += integral_imp2[j];
        }
        cumul_av /= (i + 1);
        cumul_av2 /= (i + 1);
        cumul_error = sigma(cumul_av, cumul_av2, i);
        output_imp << cumul_av << " " << cumul_error << endl;
    }
    output_imp.close();

    delete [] integral_imp;
    delete [] integral_imp2;

    return 0;
}

// Definition of the integrand function
double integrand(double x){
    return M_PI * 0.5 * cos(M_PI * 0.5 * x);
}

// Definition of the distribution for importance sampling
double distro(double x){
    return 2 * (1 - x);
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
