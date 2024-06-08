#include <iostream>   // Input-output stream library
#include <fstream>    // File stream library
#include <string>     // String manipulation library
#include <cstdlib>    // Standard library
#include <cmath>      // Math library
#include "random.h"   // Custom random number generator class

using namespace std;

int main(){
    Random rnd; // Instance of the random number generator class
    int N = 100000; // Number of variables to generate

    // Output stream to write random slopes to a file
    ofstream output;
    output.open("test.out");

    for (int i = 0; i < N; i++){
        // Generate two random points on the plane to get a random slope
        double x_1 = rnd.Rannyu();
        double y_1 = rnd.Rannyu();
        double x_2 = rnd.Rannyu();
        double y_2 = rnd.Rannyu();
        double slope = (y_2 - y_1) / (x_2 - x_1); // Calculate the slope
        output << slope << endl; // Write the slope to the file
    }
    output.close(); // Close the file

    return 0;
}
