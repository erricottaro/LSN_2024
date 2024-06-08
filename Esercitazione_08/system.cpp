#include "system.h"

System :: System() {
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

    ifstream params("INPUT/system_input.dat");
    ofstream coutf;
    coutf.open("OUTPUT/output.dat");
    while (!params.eof())
    {
        params >> property;
        if (property=="MU")
        {
            params >> _mu;
            coutf << "MU = " << _mu << endl;
        }
        else if (property=="SIGMA")
        {
            params >> _sigma;
            coutf << "SIGMA = " << _sigma << endl;
        }
        else if (property=="DELTA")
        {
            params >> _delta;
            coutf << "DELTA = " << _delta << endl;
        }
        else if (property=="NBLOCKS")
        {
            params >> _nblocks;
            coutf << "NBLOCKS = " << _nblocks << endl;
        }
        else if (property=="NSTEPS")
        {
            params >> _nsteps;
            coutf << "NSTEPS = " << _nsteps << endl;
        }
        else if (property == "ENDINPUT")
        {
            coutf << "Parameters correctly read" << endl;
        } 
        else cerr << "PROBLEM: unknown input" << endl; 
    }
    coutf << "System initialized" << endl;
    coutf.close();
    coutf.open("OUTPUT/energy.dat");
    coutf << "#     BLOCK:   ACTUAL_E:    E_AVE:      ERROR:" << endl;
    coutf.close();
    coutf.open("OUTPUT/acceptance.dat");
    coutf << "#   N_BLOCK:  ACCEPTANCE:" << endl;
    coutf.close();
    coutf.open("OUTPUT/distribution.dat");
    coutf << "   POSITION:     NORM_DISTRIB:" << endl;
    coutf.close();

    _position = 0.;
    _nbins = 200;
    _sampled_psi2 = new double[_nbins];
    for (int i = 0; i < _nbins; i++)
    {
        _sampled_psi2[i]=0;
    }
    //dimension of a bin
    _binsize = 4./double(_nbins);
    
    _average = 0.;
    _naccepted = 0;
    _nattempts = 0;
    _global_av = 0.;
    _global_av2 = 0.;
}

System :: ~System(){}

void System :: set_mu(double mu){
    _mu = mu;

    return;
}

void System :: set_sigma(double sigma){
    _sigma = sigma;

    return;
}

void System :: set_delta(double delta){
    _delta = delta;
    return;
}

double System :: get_mu(){
    return _mu;
}

double System :: get_sigma(){
    return _sigma;
}

int System :: get_nbl(){
    return _nblocks;
}

int System :: get_nsteps(){
    return _nsteps;
}

double System :: wavefunction(double x){
    double argument_minus, argument_plus;
    argument_minus = 0.5*(x - _mu)*(x - _mu)/(_sigma*_sigma);
    argument_plus  = 0.5*(x + _mu)*(x + _mu)/(_sigma*_sigma);

    return exp(-argument_minus)+exp(-argument_plus);
}

double System :: wavefunc2der(double x){
    double argument_minus, argument_plus;
    double exp_minus, exp_plus;
    double psi_minus, psi_plus;
    argument_minus = 0.5*(x - _mu)*(x - _mu)/(_sigma*_sigma);
    argument_plus  = 0.5*(x + _mu)*(x + _mu)/(_sigma*_sigma);

    exp_minus = exp(-argument_minus);
    exp_plus  = exp(-argument_plus);

    psi_minus = exp_minus/(_sigma*_sigma)*(argument_minus - 0.5);
    psi_plus  = exp_plus/(_sigma*_sigma)*(argument_plus - 0.5);

    return psi_minus + psi_plus;
}

double System :: potential(double x){
    return pow(x, 4) - 2.5*x*x;
}

void System :: step(){
    double shift;
    shift = _rnd.Rannyu(-1.0, 1.0)*_delta;
    // cout << "shift: " << shift << endl;
    if (this->metro(shift)){
        _position += shift;
        _naccepted++;
    }
    _nattempts++;
    return;
}

bool System :: metro(double shift){
    double p_new = pow(this->wavefunction(_position+shift), 2);
    double p_old = pow(this->wavefunction(_position), 2);
    double acceptance = p_new/p_old;
    bool decision = false;
    if (_rnd.Rannyu() < acceptance){ decision=true; }

    return decision;
}

energy_meas System :: finalize(){
    _rnd.SaveSeed();
    ofstream coutf;
    coutf.open("OUTPUT/output.dat", ios::app);
    coutf << "Simulation completed!" << endl;
    coutf.close();
    energy_meas last_ave;
    last_ave._energy = _global_av/double(_nblocks);
    last_ave._error = this->error(_global_av, _global_av2, _nblocks);
    _global_av = 0.;
    _global_av2 = 0.;
    _position = 0.;
    return last_ave;
}

void System :: measure(){
    double potential_energy = this->potential(_position);
    double kinetic_energy = - this->wavefunc2der(_position)/this->wavefunction(_position);

    double total_energy = kinetic_energy + potential_energy;

    _measurement = total_energy;
    _block_av += _measurement;

    //measurement of sampled distribution
    if (abs(_position)<2.){
        //determine distance from the origin
        int bin = floor(_position/_binsize);
        int bin_index = bin + _nbins/2;
        //increment counter
        _sampled_psi2[bin_index]++;
    }

    return;
}

void System :: block_reset(int blk){
    ofstream coutf;
    if (blk>0){
        coutf.open("OUTPUT/output.dat",ios::app);
        coutf << "Block completed: " << blk << endl;
        coutf.close();
    }
    _block_av=0.;
    return;
}

//to be called inside annealing to compute energies without opening and closing streams
void System :: average_wo_saving(){
    _average = _block_av / double(_nsteps);
    _global_av += _average;
    _global_av2 += _average*_average;
}

//to be called in main to collect all data
void System :: average(int blk){
    ofstream coutf;
    _average = _block_av / double(_nsteps);
    _global_av += _average;
    _global_av2 += _average*_average;

    coutf.open("OUTPUT/energy.dat", ios::app);
    double average, sum_average, sum_ave2;
    average = _average;
    sum_average = _global_av;
    sum_ave2 = _global_av2;
    coutf << setprecision(5)
        << setw(12) << blk 
        << " " << setw(12) << average
        << " " << setw(12) << sum_average/double(blk)
        << " " << setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close();

    double fraction;
    coutf.open("OUTPUT/acceptance.dat", ios::app);
    if(_nattempts > 0) fraction = double(_naccepted)/double(_nattempts);
    else fraction = 0.0; 
    coutf << setw(12) << blk << " " << setw(12) << fraction << endl;
    coutf.close();

    if (blk == _nblocks)
    {
        coutf.open("OUTPUT/distribution.dat", ios::app);
        double sum=0;
        for (int i = 0; i < _nbins; i++){
            sum+=_sampled_psi2[i];
        }
        sum*=_binsize;
        for (int i = 0; i < _nbins; i++){
            coutf << setprecision(5)
                  << setw(12) << double(i*_binsize)-2.
                  << " " << setw(12) << double(_sampled_psi2[i]/sum ) << endl;
        }
        //coutf << "TOTAL: " << sum << endl;
        coutf.close();
    }
    

    return;
}

double System :: error(double acc, double acc2, int blk){
  if(blk <= 1) return 0.0;
  else return sqrt( fabs(acc2/double(blk) - pow( acc/double(blk) ,2) )/double(blk) );
}