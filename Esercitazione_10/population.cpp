#include "population.h"

using namespace std;

population::population() {}

population::population(vector<int> &input_cities, int pop_size, int rank): _generation(pop_size) {
    int p1, p2; // Read from Primes a pair of numbers to be used to initialize the RNG
    ifstream Primes("INPUT/Primes");
    for (int j = 0; j < rank; j++) {
        string line;
        getline(Primes, line);
    } 
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

    _pop_size = _generation.size();
    //create a random population
    vector<int> dummy = input_cities;
    for (int i = 0; i < _generation.size(); i++){
        for (int j = 0; j < dummy.size()/2; j++){
            this->scramble(dummy);
        }
        
        _generation[i].set_genome(dummy);
        /*
        cout << "Individual # " << i+1 << endl;
        cout << "Check: " << this->check(_generation[i]) << endl;
        for (int j = 0; j < _generation[i].get_ncities(); j++){
          cout << "City index: " << _generation[i].get_gene(j) << endl;
        } 
        */
    }
    
}

population::~population(){}

int population::get_pop_size(){
    return _pop_size;
}

individual population::get_individual(int i){
    return _generation[i];
}

void population::set_individual(individual new_indiv, int i){
    vector<int> new_genome = new_indiv.get_genome();
    _generation[i].set_genome(new_genome);
    return;
}

void population::order_population(vector<double> &fitness){
    // Create a vector of pairs (individual, key)
    vector<pair<individual, double>> paired;
    for (size_t i = 0; i < _generation.size(); ++i) {
        paired.emplace_back(_generation[i], fitness[i]);
    }

    // Sort the vector of pairs based on the key
    sort(paired.begin(), paired.end(), [](const pair<individual, double>& a, const pair<individual, double>& b) {return a.second < b.second;});

    // Extract the sorted _generation based on the sorted pairs
    for (size_t i = 0; i < _generation.size(); ++i) {
        _generation[i] = paired[i].first;
    }
    return;
}

void population::selection(double expo){
    vector<individual> new_generation(_generation.size());
    for (int i = 0; i < _generation.size(); i++) {
        int index = int(pow(_rnd.Rannyu(), expo));
        new_generation[i] = _generation[index];
    }
    this->replace_population(new_generation);
    return;
}

void population::scramble(vector<int> &genes){
    int first_index = int(_rnd.Rannyu(1, genes.size()));
    int second_index = int(_rnd.Rannyu(1, genes.size()));
    swap(genes[first_index], genes[second_index]);
    return;
}

void population::replace_population(vector<individual> &new_gen){
    for (int i = 0; i < _generation.size(); i++){
        _generation[i] = new_gen[i];
    }
    return;
}

bool population::check(individual ind){
    bool check = true;
    vector<int> genome = ind.get_genome();
    if (genome[0]!=0){
        check = false;
    }
    vector<int>::iterator it = unique(genome.begin(), genome.end());
    if (it!=genome.end()){
        check = false;
    }
    
    return check;
}

//mutations on individuals
void population::swap_mutation(int index){
    vector<int> dummy = _generation[index].get_genome();
    this->scramble(dummy);
    _generation[index].set_genome(dummy);
    // cout << this->check(_generation[index]) << endl;
    return;
}

void population::shift_mutation(int index){
    vector<int> genome = _generation[index].get_genome();
    int string_index = _rnd.Rannyu(1, genome.size()-1);
    int string_length = _rnd.Rannyu(1, genome.size()-string_index);
    //number of steps to shift the string to the right
    int shift = _rnd.Rannyu(1, genome.size()-string_index-string_length);
    //shift the elements
    rotate(genome.begin()+string_index, genome.begin()+string_index+string_length, genome.end());
    rotate(genome.begin()+string_index+string_length, genome.end()-shift, genome.end());
    _generation[index].set_genome(genome);
    // cout << this->check(_generation[index]) << endl;
}

void population::swap_group_mutation(int index){
    vector<int> genome = _generation[index].get_genome();
    //create index in the first half (first gene must remain fixed)
    int string_index;
    double r = _rnd.Rannyu();
    if (r < 0.5) {
        //half of the times, the string is made to start at the first available gene, and is swapped to the right
        string_index = 1;
        //create length of the string of genes to swap
        int string_length = int(_rnd.Rannyu(1, genome.size()/2));
        //starting index of the second string of genes to swap
        int string_index_two = string_index+string_length;
        swap_ranges(genome.begin()+string_index, genome.begin()+string_index+string_length, genome.begin()+string_index_two);
        
    } 
    else {
        //half of the times, the string is made to start at the first available gene, and is swapped to the left
        string_index = genome.size()/2;
        // cout << string_index << endl << endl;
        //create length of the string of genes to swap
        int string_length = int(_rnd.Rannyu(1, genome.size()/2));
        // cout << string_length << endl << endl;
        //starting index of the second string of genes to swap
        int string_index_two = string_index - string_length;
        swap_ranges(genome.begin()+string_index, genome.begin()+string_index+string_length, genome.begin()+string_index_two);
    }
    /*
    for (int i = 0; i < genome.size(); i++){
        cout << genome[i] << endl;
    } 
    cout << endl; */
    _generation[index].set_genome(genome);

    // cout << this->check(_generation[index]) << endl;
    return;
}

void population::inversion_mutation(int index){
    vector<int> genome = _generation[index].get_genome();
    int string_index = int(_rnd.Rannyu(1, genome.size()/2));
    int string_length = int(_rnd.Rannyu(1, genome.size()-string_index));
    reverse(genome.begin()+string_index, genome.begin()+string_index+string_length);
    _generation[index].set_genome(genome);
    // cout << this->check(_generation[index]) << endl;
    return;
}

//mutations on all population
void population::swap_mutation(double prob){
    for (int i = 0; i < _pop_size; i++){
        double r = _rnd.Rannyu();
        if (r < prob){ this->swap_mutation(i);}
    }
    return;
}

void population::shift_mutation(double prob){
    for (int i = 0; i < _pop_size; i++){
        double r = _rnd.Rannyu();
        if (r < prob) { this->shift_mutation(i); }    
    }
    return;
}

void population::swap_group_mutation(double prob){
    for (int i = 0; i < _pop_size; i++){
        double r = _rnd.Rannyu();
        if (r < prob){ this->swap_group_mutation(i);}
    }
    return;
}

void population::inversion_mutation(double prob){
    for (int i = 0; i < _pop_size; i++){
        double r = _rnd.Rannyu();
        if (r < prob){ this->inversion_mutation(i);}
    }
    return;
}

//crossing
void population::crossing(int ind1, int ind2){
    vector<int> genome1 = _generation[ind1].get_genome();
    vector<int> genome2 = _generation[ind2].get_genome();
    int genome_size = genome1.size();
    //generate random index in correspondence of which to slice
    int slicing_index = _rnd.Rannyu(1, genome_size);
    //generate slices
    vector<int> slice1(genome1.begin()+slicing_index, genome1.end());
    vector<int> slice2(genome2.begin()+slicing_index, genome2.end());
    //arrays of ranks
    vector<int> rank1=rank_order(slice1);
    vector<int> rank2=rank_order(slice2);
    //reorder the slices according to the partner's ranks
    vector<int> new_slice1(slice1.size());    
    vector<int> new_slice2(slice2.size());
    for (int i = 0; i < rank1.size(); i++) {
        //find in the array of ranks the index of the element which is equal to the current element of the other array of ranks
        auto it1 = find(rank1.begin(), rank1.end(), rank2[i]);
        auto it2 = find(rank2.begin(), rank2.end(), rank1[i]);
        int index1 = distance(rank1.begin(), it1);
        int index2 = distance(rank2.begin(), it2);
        //now fill the new slices
        new_slice1[i] = slice1[index1];
        new_slice2[i] = slice2[index2];
    }
    //substitute genomes of the two individuals
    copy(new_slice1.begin(), new_slice1.end(), genome1.begin()+slicing_index);
    copy(new_slice2.begin(), new_slice2.end(), genome2.begin()+slicing_index);
    //update individuals
    _generation[ind1].set_genome(genome1);
    _generation[ind2].set_genome(genome2);
    // cout << this->check(_generation[ind1]) << endl;
    // cout << this->check(_generation[ind2]) << endl;
    return;
}

void population::crossing(double prob){
    for (int i = 0; i < _pop_size; i+=2){
        double r = _rnd.Rannyu();
        if (r < prob) { this->crossing(i, i+1); }
    }
    return;
}


//used in crossing
vector<int> rank_order(const vector<int>& vec) {
    // Create a vector of pairs (value, original_index)
    vector<pair<int, int>> value_index_pairs;
    value_index_pairs.reserve(vec.size());

    for (int i = 0; i < vec.size(); ++i) {
        value_index_pairs.emplace_back(vec[i], i);
    }

    // Sort the pairs based on the values
    sort(value_index_pairs.begin(), value_index_pairs.end());

    // Create a vector to store the rank of each element
    vector<int> ranks(vec.size());

    for (int rank = 0; rank < value_index_pairs.size(); ++rank) {
        ranks[value_index_pairs[rank].second] = rank;
    }

    return ranks;
}