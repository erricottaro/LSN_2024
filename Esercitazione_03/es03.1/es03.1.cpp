#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include "random.h"

using namespace std;

// Function to compute statistical error
double sigma(double, double, int);
// Function to compute call option price
double call_price(double, double, double, double);
// Function to compute put option price  
double put_price(double, double, double, double);

int main() {
    Random rnd;

    // Parameters for the option pricing
    double S_0 = 100.;    // Initial asset price
    double T = 1;         // Time of delivery
    double K = 100;       // Strike price
    double r = 0.1;       // Risk-free interest rate
    double volat = 0.25;  // Volatility

    // Save parameters to a file
    ofstream param;
    param.open("economic_parameters.dat");
    param << S_0 << " " << T << " " << K << " " << r << " " << volat << endl;
    param.close();

    // Simulation parameters for one-step method
    int M_onestep = 10000000;  // Number of throws for one step method
    int N = 100;                // Number of blocks
    int m_onestep = M_onestep / N; // Number of throws per block

    // Save simulation parameters to a file
    param.open("sim_onestep_parameters.dat");
    param << M_onestep << " " << N << endl;
    param.close();

    // Arrays to store prices and their squares for call and put options with one-step method
    double* prices_call_onestep = new double[N];
    double* prices_call_onestep2 = new double[N];
    double* prices_put_onestep = new double[N];
    double* prices_put_onestep2 = new double[N];

    // Simulation for call and put options using one-step method
    for (int i = 0; i < N; i++) {
        prices_call_onestep[i] = 0;
        prices_call_onestep2[i] = 0;
        prices_put_onestep[i] = 0;
        prices_put_onestep2[i] = 0;

        for (int j = 0; j < m_onestep; j++) {
            // Generate a random fluctuation for the asset price for call option
            double Z_T_call = rnd.Gauss(0., 1.);
            double S_T_call = S_0 * exp((r - 0.5 * volat * volat) * T + volat * Z_T_call * sqrt(T));
            prices_call_onestep[i] += call_price(S_T_call, r, K, T);

            // Generate a random fluctuation for the asset price for put option
            double Z_T_put = rnd.Gauss(0., 1.);
            double S_T_put = S_0 * exp((r - 0.5 * volat * volat) * T + volat * Z_T_put * sqrt(T));
            prices_put_onestep[i] += put_price(S_T_put, r, K, T);
        }
        prices_call_onestep[i] /= m_onestep;
        prices_call_onestep2[i] = pow(prices_call_onestep[i], 2);
        prices_put_onestep[i] /= m_onestep;
        prices_put_onestep2[i] = pow(prices_put_onestep[i], 2);
    }

    // Computation of averages and error for one-step method for call options
    ofstream output_call_onestep;
    output_call_onestep.open("call_onestep.out");
    for (int i = 0; i < N; i++) {
        double cumul_av_call_onestep = 0;
        double cumul_av2_call_onestep = 0;
        double cumul_error_call_onestep = 0;

        for (int j = 0; j < i + 1; j++) {
            cumul_av_call_onestep += prices_call_onestep[j];
            cumul_av2_call_onestep += prices_call_onestep2[j];
        }
        cumul_av_call_onestep /= (i + 1);
        cumul_av2_call_onestep /= (i + 1);
        cumul_error_call_onestep = sigma(cumul_av_call_onestep, cumul_av2_call_onestep, i);
        output_call_onestep << cumul_av_call_onestep << " " << cumul_error_call_onestep << endl;
    }
    output_call_onestep.close();

    // Computation of averages and error for one-step method for put options
    ofstream output_put_onestep;
    output_put_onestep.open("put_onestep.out");
    for (int i = 0; i < N; i++) {
        double cumul_av_put_onestep = 0;
        double cumul_av2_put_onestep = 0;
        double cumul_error_put_onestep = 0;

        for (int j = 0; j < i + 1; j++) {
            cumul_av_put_onestep += prices_put_onestep[j];
            cumul_av2_put_onestep += prices_put_onestep2[j];
        }
        cumul_av_put_onestep /= (i + 1);
        cumul_av2_put_onestep /= (i + 1);
        cumul_error_put_onestep = sigma(cumul_av_put_onestep, cumul_av2_put_onestep, i);
        output_put_onestep << cumul_av_put_onestep << " " << cumul_error_put_onestep << endl;
    }
    output_put_onestep.close();

    // Cleanup memory
    delete[] prices_call_onestep;
    delete[] prices_call_onestep2;
    delete[] prices_put_onestep;
    delete[] prices_put_onestep2;

    // Simulation parameters for discrete walk method
    int M_walk = 1000000;
    int m_walk = M_walk / N;
    int walk_length = 100;
    double time_int = T / walk_length;

    // Save simulation parameters for discrete walk to a file
    param.open("sim_walk_parameters.dat");
    param << M_walk << " " << N << " " << walk_length << endl;
    param.close();

    // Arrays to store prices and their squares for call and put options with discrete walk
    double* prices_call_walk = new double[N];
    double* prices_call_walk2 = new double[N];
    double* prices_put_walk = new double[N];
    double* prices_put_walk2 = new double[N];

    // Simulation for call and put options using discrete walk
    for (int i = 0; i < N; i++) {
        prices_call_walk[i] = 0;
        prices_call_walk2[i] = 0;
        prices_put_walk[i] = 0;
        prices_put_walk2[i] = 0;

        for (int j = 0; j < m_walk; j++) {
            // Recursive method for call option
            double Z_T_call_i = rnd.Gauss(0, 1);
            double S_call = S_0 * exp((r - 0.5 * volat * volat) * time_int + volat * Z_T_call_i * sqrt(time_int));
            for (int t = 0; t < walk_length - 1; t++) {
                Z_T_call_i = rnd.Gauss(0, 1);
                S_call *= exp((r - 0.5 * volat * volat) * time_int + volat * Z_T_call_i * sqrt(time_int));
            }
            prices_call_walk[i] += call_price(S_call, r, K, T);

            // Recursive method for put option
            double Z_T_put_i = rnd.Gauss(0, 1);
            double S_put = S_0 * exp((r - 0.5 * volat * volat) * time_int + volat * Z_T_put_i * sqrt(time_int));
            for (int t = 0; t < walk_length - 1; t++) {
                Z_T_put_i = rnd.Gauss(0, 1);
                S_put *= exp((r - 0.5 * volat * volat) * time_int + volat * Z_T_put_i * sqrt(time_int));
            }
            prices_put_walk[i] += put_price(S_call, r, K, T);
        }
        prices_call_walk[i] /= m_walk;
        prices_call_walk2[i] = pow(prices_call_walk[i], 2);
        prices_put_walk[i] /= m_walk;
        prices_put_walk2[i] = pow(prices_put_walk[i], 2);
    }

    // Computation of averages and error for discrete walk method for call options
    ofstream output_call_walk;
    output_call_walk.open("call_walk.out");
    for (int i = 0; i < N; i++) {
        double cumul_av_call_walk = 0;
        double cumul_av2_call_walk = 0;
        double cumul_error_call_walk = 0;

        for (int j = 0; j < i + 1; j++) {
            cumul_av_call_walk += prices_call_walk[j];
            cumul_av2_call_walk += prices_call_walk2[j];
        }
        cumul_av_call_walk /= (i + 1);
        cumul_av2_call_walk /= (i + 1);
        cumul_error_call_walk = sigma(cumul_av_call_walk, cumul_av2_call_walk, i);
        output_call_walk << cumul_av_call_walk << " " << cumul_error_call_walk << endl;
    }
    output_call_walk.close();

    // Computation of averages and error for discrete walk method for put options
    ofstream output_put_walk;
    output_put_walk.open("put_walk.out");
    for (int i = 0; i < N; i++) {
        double cumul_av_put_walk = 0;
        double cumul_av2_put_walk = 0;
        double cumul_error_put_walk = 0;

        for (int j = 0; j < i + 1; j++) {
            cumul_av_put_walk += prices_put_walk[j];
            cumul_av2_put_walk += prices_put_walk2[j];
        }
        cumul_av_put_walk /= (i + 1);
        cumul_av2_put_walk /= (i + 1);
        cumul_error_put_walk = sigma(cumul_av_put_walk, cumul_av2_put_walk, i);
        output_put_walk << cumul_av_put_walk << " " << cumul_error_put_walk << endl;
    }
    output_put_walk.close();

    // Cleanup memory
    delete[] prices_call_walk;
    delete[] prices_call_walk2;
    delete[] prices_put_walk;
    delete[] prices_put_walk2;

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

// Function to compute call option price
double call_price(double S, double r, double K, double T) {
    if (S > K) {
        return exp(-r * T) * (S - K);
    } else {
        return 0;
    }
}

// Function to compute put option price
double put_price(double S, double r, double K, double T) {
    if (S < K) {
        return exp(-r * T) * (K - S);
    } else {
        return 0;
    }
}
