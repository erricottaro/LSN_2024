#include "individual.h"
#include "random.h"
#include <vector>
#include <algorithm>
#include <fstream>

using namespace std;

#ifndef __population__
#define __population__

class population {
private:
    int _pop_size;                        // Size of the population
    vector<individual> _generation;       // Vector containing the current generation of individuals
    Random _rnd;                          // Random number generator

public:
    // Constructor generates a random population given a vector of cities
    population(vector<int> &input_cities, int pop_size);
    ~population();                        // Destructor

    int get_pop_size();                   // Returns the size of the population
    individual get_individual(int i);     // Returns the individual at index i

    // Orders the population according to fitness and applies selection operator
    void order_population(vector<double> &fitness);
    void selection(double expo);          // Applies selection based on a given exponent parameter

    // Scrambles the genes leaving the first one unchanged (the starting city remains the same)
    void scramble(vector<int> &genes);

    // Replaces the current population with a new generation
    void replace_population(vector<individual> &new_gen);

    // Checks the validity of an individual
    bool check(individual ind);

    // Mutation that acts on a single individual of the population
    void swap_mutation(int index);        // Swaps two genes in the individual at index
    void shift_mutation(int index);       // Shifts genes in the individual at index
    void swap_group_mutation(int index);  // Swaps groups of genes in the individual at index
    void inversion_mutation(int index);   // Inverts a sequence of genes in the individual at index

    // Mutation that acts on the entire population with a given probability
    void swap_mutation(double prob);      // Applies swap mutation to the population with a given probability
    void shift_mutation(double prob);     // Applies shift mutation to the population with a given probability
    void swap_group_mutation(double prob);// Applies group swap mutation to the population with a given probability
    void inversion_mutation(double prob); // Applies inversion mutation to the population with a given probability

    // Crossing (crossover) operator
    void crossing(int ind1, int ind2);    // Applies crossover between two individuals
    void crossing(double prob);           // Applies crossover to the population with a given probability
};


#endif

vector<int> rank_order(const vector<int>& vec); //auxiliary function used in the crossing operator