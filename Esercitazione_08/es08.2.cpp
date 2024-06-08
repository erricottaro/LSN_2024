#include <iostream>
#include <cmath>
#include "annealing.h"

using namespace std;

int main(){ //change of strategy. cycle over temperatures in a python script

    annealing simulator;

    //thermalization of annealing metropolis (this one is needed)
    for (int i = 0; i < 100; i++){
        simulator.step();
    }
    for (int i = 0; i < simulator.get_nsteps(); i++){ //loop over total number of steps
        simulator.step();
        simulator.measure();
    }
    simulator.finalize();

    ofstream final_params;
    final_params.open("final_parameters.dat");

    final_params << "MU: " << simulator.get_mu() << endl;
    final_params << "SIGMA: " << simulator.get_sigma() << endl;

    final_params.close();

    return 0;
}
