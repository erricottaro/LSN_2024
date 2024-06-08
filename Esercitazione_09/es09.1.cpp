#include "individual.h"
#include "population.h"
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>

using namespace std;

double compute_fitness(individual, vector<city> &); //compute fitness of a given individual (length of the route)
double distancesq(city, city); //compute distance squared between two cities
int pbc(int i, int dim); //periodic boundary condition over indeces

int main(){
    Random rnd;
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
    rnd.SetRandom(seed,p1,p2);
    Seed.close();

    //generate array of cities on a circle of radius one
    int n_cities = 34;
    vector<city> cities(n_cities);
    vector<int>  labels(n_cities);
    
    for (int i = 0; i < cities.size(); i++){
        double theta = rnd.Rannyu(0. , 2.*M_PI);
        cities[i].set_coord(cos(theta), 0);
        cities[i].set_coord(sin(theta), 1);
        labels[i]=i;
        
    }
    
    //save them on a file
    ofstream ocities("OUTPUT/cities_circle.dat");
    ocities << "     INDEX" << "          x" << "      y" << endl;
    for (int i = 0; i < cities.size(); i++){
        ocities << setw(8) << i << setw(15) << cities[i].get_coord(0) << setw(15) << cities[i].get_coord(1) << endl;
    }
    ocities.close();

    int n_generations = 1000;              // Number of generations for the simulation
    int pop_size = 500;                    // Population size
    double mut_prob_circle = 0.05;         // Probability for mutation
    double cross_prob_circle = 0.95;       // Probability for crossover
    double exponent_circle = 2.;           // Exponent used in the selection process

    population pop(labels, pop_size);      // Initialize the population with given labels and population size
    vector<double> fitness(pop_size);      // Vector to store fitness values of individuals

    ofstream ofitness("OUTPUT/fitness_circle.dat");     // Output file for average fitness per generation
    ofitness << "       GEN" << "           L2" << endl; // Write header to the fitness file

    ofstream obestfitness("OUTPUT/best_fitness_circle.dat"); // Output file for best fitness per generation
    obestfitness << "       GEN" << "           L2" << endl; // Write header to the best fitness file

    for (int i = 0; i < n_generations; i++) {
        // Compute fitness of each individual
        for (int j = 0; j < pop_size; j++) {
            individual ind = pop.get_individual(j);       // Get the individual at index j
            fitness[j] = compute_fitness(ind, cities);   // Compute the fitness of the individual
        }
        // Order population according to its fitness
        pop.order_population(fitness);
        // Apply selection operator
        pop.selection(exponent_circle);
        // Apply crossing on the population
        pop.crossing(cross_prob_circle);
        // Apply mutations on the population
        pop.swap_mutation(mut_prob_circle);
        pop.shift_mutation(mut_prob_circle);
        pop.swap_group_mutation(mut_prob_circle);
        pop.inversion_mutation(mut_prob_circle);

        // Sort fitness values to find the best and calculate the average of the first half
        sort(fitness.begin(), fitness.end());
        // Print best fitness
        obestfitness << setw(8) << i << setw(15) << fitness[0] << endl;
        double accu = 0.;
        for (int j = 0; j < pop_size / 2; j++) {
            accu += fitness[j];   // Accumulate the fitness values of the first half
        }
        accu /= double(pop_size) / 2;  // Compute average fitness of the first half
        // Print average fitness to the file
        ofitness << setw(8) << i << setw(15) << accu << endl;
    }

    ofitness.close();         // Close the average fitness output file
    obestfitness.close();     // Close the best fitness output file

    // Now find the best individual among the last generation
    ofstream ogenes("OUTPUT/genes_circle.dat");
    ogenes << "      GENE" << endl; // Write header to the genes output file

    // Compute fitness for all individuals again
    for (int j = 0; j < pop_size; j++) {
        individual ind = pop.get_individual(j);       // Get the individual at index j
        fitness[j] = compute_fitness(ind, cities);   // Compute the fitness of the individual
    }
    pop.order_population(fitness);  // Order population by fitness
    individual best = pop.get_individual(0);  // Get the best individual

    // Write the genes of the best individual to the file
    for (int i = 0; i < best.get_genome_size(); i++) {
        ogenes << setw(8) << best.get_gene(i) << endl;
    }
    ogenes.close();  // Close the genes output file
    //*********************************************************//
    //*********************************************************//
    //repeat for cities in a square [0, 1]x[0, 1]
    for (int i = 0; i < cities.size(); i++){
    double x = rnd.Rannyu();
    double y = rnd.Rannyu();
    cities[i].set_coord(x, 0);
    cities[i].set_coord(y, 1);
    labels[i]=i;
        
    }
    //save them on file
    ocities.open("OUTPUT/cities_square.dat");
    ocities << "     INDEX" << "          x" << "      y" << endl;
    for (int i = 0; i < cities.size(); i++){
        ocities << setw(8) << i << setw(15) << cities[i].get_coord(0) << setw(15) << cities[i].get_coord(1) << endl;
    }
    ocities.close();

    double mut_prob_square = 0.11;
    double cross_prob_square = 0.95;
    double exponent_square = 2.;

    //new population
    population pop_2(labels, pop_size);
    ofitness.open("OUTPUT/fitness_square.dat");
    ofitness << "       GEN" << "           L2" << endl;
    obestfitness.open("OUTPUT/best_fitness_square.dat");
    obestfitness << "       GEN" << "           L2" << endl;

    for (int i = 0; i < n_generations; i++) {
        //cout << i << "th generation running..." << endl;
        //compute fitness of each individual
        for (int j = 0; j < pop_size; j++){
            individual ind = pop_2.get_individual(j);
            fitness[j] = compute_fitness(ind, cities);
        }
        //order population according to its fitness
        pop_2.order_population(fitness);
        pop_2.selection(exponent_square);
        //apply crossing on the population
        pop_2.crossing(cross_prob_square);
        //Apply mutations on the population
        pop_2.swap_mutation(mut_prob_square);
        pop_2.shift_mutation(mut_prob_square);
        pop_2.swap_group_mutation(mut_prob_square);
        pop_2.inversion_mutation(mut_prob_square); 

        //sort fitness and compute average of the first half
        sort(fitness.begin(), fitness.end());
        //print best fitness
        obestfitness << setw(8) << i << setw(15) << fitness[0] << endl;
        double accu = 0.;
        for (int j = 0; j < pop_size/2; j++){
            accu += fitness[j];
        }
        accu/= double(pop_size)/2;
        //print it on file
        ofitness << setw(8) << i << setw(15) << accu << endl;
        
    }
    ofitness.close();
    obestfitness.close();

    //now find the best among all
    ogenes.open("OUTPUT/genes_square.dat");
    ogenes << "      GENE" << endl;

    for (int j = 0; j < pop_size; j++){
        individual ind = pop_2.get_individual(j);
        fitness[j] = compute_fitness(ind, cities);
    }
    pop_2.order_population(fitness);
    best = pop_2.get_individual(0);
    for (int i = 0; i < best.get_genome_size(); i++){
        ogenes << setw(8) << best.get_gene(i) << endl;
    }
    ogenes.close();

    return 0;
}

double compute_fitness(individual indiv, vector<city> &cities){
    double sum2 = 0.;
    if (indiv.get_genome_size()!=cities.size()){
        cerr << "ERROR: #of cities not corresponding" << endl;
        exit(-1);
    }
    for (int i = 0; i < cities.size(); i++){
        sum2 += distancesq(cities[indiv.get_gene(pbc(i, cities.size()))], cities[indiv.get_gene(pbc(i+1, cities.size()))]);
    }
    return sum2;
}

double distancesq(city city1, city city2){
    double sum2 = 0.;
    for (int i = 0; i < 2; i++)
    {
        sum2 += pow(city1.get_coord(i) - city2.get_coord(i), 2);
    }
    return sum2;
}

int pbc(int i, int dim){
    return i - int(i/dim)*dim;
}