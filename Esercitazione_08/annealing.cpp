#include "annealing.h"

using namespace std;

annealing :: annealing() : _SYS() { //constructor that initializes all parameters. It also calls _SYS constructor
    int p1, p2; // Read from Primes a pair of numbers to be used to initialize the RNG
    ifstream Primes("INPUT/Primes");
    Primes >> p1 >> p2 ;
    Primes.close();
    int seed[4]; // Read the seed of the RNG
    ifstream Seed("INPUT/seed.in");
    string property;
    Seed >> property;
    if( property == "RANDOMSEED" ){
            Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
         }
    _rnd.SetRandom(seed,p1,p2);
    Seed.close();


    ofstream coutf;
    coutf.open("OUTPUT/annealing_output.dat");
    //reading input file
    ifstream input("INPUT/annealing_input.dat");
    while (!input.eof())
    {
        input >> property;
        if (property=="TEMP"){
            input >> _temp;
            _beta = 1./_temp;
            coutf << "TEMP: " << _temp << endl;
        }
        else if (property=="DELTAMU"){
            input >> _delta_mu;
            coutf << "DELTAMU: " << _delta_mu << endl;
        }
        else if (property=="DELTASIGMA"){
            input >> _delta_sigma;
            coutf << "DELTASIGMA: " << _delta_sigma << endl;
        }
        else if (property=="NSTEPS"){
            input >> _nsteps;
            coutf << "NSTEPS: " << _nsteps << endl;
        }
        else if (property=="ENDINPUT"){
            coutf << "Parameters correctly read" << endl;
        }
        else cerr << "PROBLEM: unknown input" << endl;
    }
    coutf << "System initialized" << endl;
    coutf.close();
    coutf.open("OUTPUT/annealing_energy.dat");
    coutf << "     MU:      SIGMA:     E_MEASURED:    ERROR:" << endl;
    coutf.close();
    coutf.open("OUTPUT/annealing_accep.dat");
    coutf << "     ACCEPTANCE:" << endl;
    coutf.close();

    _average = 0.;
    _naccepted = 0.;
    _nattempts = 0.;
    //first computation of the energy
    _new_energy = this->get_energy();
}

annealing :: ~annealing(){}

int annealing :: get_nsteps(){
    return _nsteps;
}

double annealing :: get_beta(){
    return _beta;
}

double annealing :: get_mu(){
    return _SYS.get_mu();
}

double annealing :: get_sigma(){
    return _SYS.get_sigma();
}

void annealing :: step() {
    //old energy is previous new energy
    _old_energy = _new_energy;
    //generate move in parameter space
    double shift_mu = _rnd.Rannyu(-1.0, 1.0)*_delta_mu;
    double shift_sigma = _rnd.Rannyu(-1.0, 1.0)*_delta_sigma;
    //set new parameters
    double new_mu = _SYS.get_mu() + shift_mu;
    double new_sigma = _SYS.get_sigma() + shift_sigma;
    _SYS.set_mu(new_mu);
    _SYS.set_sigma(new_sigma);
    _new_energy = this->get_energy(); //only time that compute energy is called 
    if (this->metro()){
        //if move accepted, everything has already been done
        _naccepted++;
    }
    else {
        //restore previous parameters
        _SYS.set_mu(new_mu - shift_mu);
        _SYS.set_sigma(new_sigma - shift_sigma);
        //the new energy must be retrieved
        _new_energy = _old_energy;
    }
    _nattempts++;
    return;
}

bool annealing :: metro() {
    //Compute boltzmann factors ratio between energy of new configuration and old one
    double acceptance =exp(-_beta*(_new_energy._energy - _old_energy._energy));
    bool decision = false;
    if (_rnd.Rannyu() < acceptance){ decision=true; }

    return decision;
}

void annealing :: finalize() {
    _rnd.SaveSeed();
    ofstream coutf;
    coutf.open("OUTPUT/annealing_output.dat", ios::app);
    coutf << "Simulation completed!" << endl;
    coutf.close();
    return;
}

energy_meas annealing :: get_energy() { //very time consuming method (set params of _SYS such that error on this measure is in % order)
    //here we compute the total energy by a complete MC simulation
    _SYS.block_reset(0);
    
    for (int i = 0; i < _SYS.get_nbl(); i++){
        for (int j = 0; j < _SYS.get_nsteps(); j++){
            _SYS.step();
            _SYS.measure();
        }
        _SYS.average_wo_saving();
        _SYS.block_reset(i+1);
    }
    energy_meas energy = _SYS.finalize();

    return energy;
}

void annealing :: measure() {
    //energy has already been computed during Metropolis move
    ofstream coutf;
    coutf.open("OUTPUT/annealing_energy.dat", ios::app);

    coutf << setprecision(5)
          << setw(12) << _SYS.get_mu()
          << setw(12) << _SYS.get_sigma() 
          << setw(12) << _new_energy._energy
          << setw(12) << _new_energy._error << endl;
    coutf.close();

    double fraction;
    coutf.open("OUTPUT/annealing_accep.dat", ios::app);
    if(_nattempts > 0) fraction = double(_naccepted)/double(_nattempts);
    else fraction = 0.0; 
    coutf << " " << setw(12) << fraction << endl;
    coutf.close();
    return;
}

double annealing :: error(double acc, double acc2, int blk) {
    if(blk <= 1) return 0.0;
    else return sqrt( fabs(acc2/double(blk) - pow( acc/double(blk) ,2) )/double(blk) );
}
