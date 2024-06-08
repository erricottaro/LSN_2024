#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>

using namespace std;

// Function to compute autocorrelation
double autocorrelation(int t, vector<double> data);

// Function to compute variance
double variance(vector<double> data);

int main(){
    ifstream energies;
    ofstream autocorrs;
    string phases[3] = {"Solid", "Liquid", "Gas"};
    
    // Loop through different phases
    for (int i = 0; i < 3; i++)
    {  
        vector<double> data;
        string inputfile = phases[i]+"/potential_energy.dat";
        cout << inputfile << endl;
        
        // Open input file
        energies.open(inputfile);
        
        // Read and discard the first line
        string firstline;
        getline(energies, firstline);
        
        double block, appo_energy, average, error;
        
        // Read data from file
        for(int j=0; j< 500000; j++)
        {
            energies >> block >> appo_energy >> average >> error;
            data.push_back(appo_energy);
        }
        
        // Close input file
        energies.close();
        cout << data.size() << endl;
        
        // Output file for autocorrelation
        string outputfile = phases[i]+"/autocorrelation.dat";
        autocorrs.open(outputfile);
        autocorrs << "      T:          AUTOCORRELATION:" << endl;
        autocorrs << setprecision(10);
        double autoc;
        
        // Compute variance of the data
        double var = variance(data);
        
        // Compute autocorrelation
        for (int j = 0; j < data.size()/1000; j++)
        {
            autoc = autocorrelation(j, data);
            autoc /= var;
            autocorrs << setw(10) << j
                      << " " << setw(10) << autoc << endl; 
        }
        
        // Close output file
        autocorrs.close(); 
    }
    
    return 0;
}

// Function to compute autocorrelation
double autocorrelation(int t, vector<double> data) { 
    int n = data.size(); // Dimension of the dataset
    int t_max = n-1; // Max value of the index

    // Return 0 if t_max is equal to t
    if (t_max==t) { return 0; }

    int time_difference = t_max-t;
    double sum_product  = 0.;
    double sum_tauplus  = 0.;
    double sum_tau      = 0.;

    // Compute sum of products and sums
    for (int i = 0; i <= time_difference; i++)
    {
        sum_product     += data[i]*data[i+t];
        sum_tau += data[i];
        sum_tauplus     += data[i+t];
    }
    double numerator = (sum_product - sum_tau*sum_tauplus/double(time_difference))/double(time_difference);

    return numerator;
}

// Function to compute variance
double variance(vector<double> data) {
    int n = data.size();
    
    // Return 0 if there's only one element
    if (n==1) { return 0; }
    
    double sum_sq = 0;
    double sum    = 0;
    
    // Compute sum of squares and sum of elements
    for (int i = 0; i < n; i++)
    {
        sum    += data[i];
        sum_sq += data[i]*data[i];
    }
    
    // Compute and return variance
    sum_sq /= double(n-1);
    sum    /= double(n-1);
    return sum_sq - sum*sum;
}
