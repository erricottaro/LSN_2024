#include "individual.h"

using namespace std;

individual::individual() {}

individual::~individual(){}

int individual::get_gene(int i){
    return _genes[i];
}

vector<int> individual::get_genome(){
    return _genes;
}

int individual::get_genome_size(){
    return _genome_size;
}

void individual::set_gene(int new_city, int i){
    _genes[i] = new_city;
}

void individual::set_genome(vector<int> &input_cities){
    _genes = input_cities;
    _genome_size = _genes.size();
}
