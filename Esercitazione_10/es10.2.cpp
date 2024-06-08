#include "individual.h"
#include "population.h"
#include "mpi.h"
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>

using namespace std;

double compute_fitness(individual, vector<city> &);   // Function to compute fitness of an individual
double distancesq(city, city);                        // Function to compute squared distance between two cities
void permutation(vector<int>, Random);                // Function to create a permutation of integers
int pbc(int i, int dim);                              // Function to apply periodic boundary conditions

int main(int argc, char* argv[]) {
    // Initialize MPI environment
    int size, rank;
    MPI_Init(&argc, &argv);
    // Get the total number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    // Get the rank (ID) of this process
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);            

    // Vector to store ranks of all processes
    vector<int> vec_ranks(size);                      
    for (int i = 0; i < size; i++) {
        // Initialize vec_ranks with rank IDs
        vec_ranks[i] = i;                             
    }

    // Random number generator instance
    Random rnd;                                       
    int p1, p2;                                      
    ifstream Primes("INPUT/Primes");                

    // Read lines up to the rank of this process
    for (int j = 0; j < rank; j++) {
        string line;
        getline(Primes, line);                        
    }

    // Read the primes for this process
    Primes >> p1 >> p2;                               
    Primes.close();                                   
    int seed[4];                                      
    ifstream Seed("INPUT/seed.in");                   
    string property;
    Seed >> property;                                 
    if (property == "RANDOMSEED") {
        Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    }
    // Initialize the RNG with seed and primes
    rnd.SetRandom(seed, p1, p2);                      
    Seed.close();                                    

    // Load cities from file
    ifstream cap_prov("cap_prov_ita.dat");            
    vector<city> cities;                              
    vector<int> labels;                              
    int i = 0;

    while (!cap_prov.eof()) {
        // Read coordinates of a city
        double x, y;
        cap_prov >> x >> y;                           
        city dummy(x, y);                             
        cities.push_back(dummy);                      
        labels.push_back(i);                          
        i++;
    }
    // Remove the last city which is duplicated
    cities.erase(cities.end());
    int n_cities = cities.size();                     
    labels.erase(labels.end() - 1);                   

    // Open the input file for simulation parameters
    ifstream input("INPUT/input.dat");                
    int n_generations, pop_size, n_migrations;
    bool migrate; //if true, processes exchange individuals, otherwise run independent simulations
    double mut_prob, cross_prob;
    while (!input.eof()) {
        input >> property;
        if (property == "GENERATIONS") {
            // Read the number of generations
            input >> n_generations;                   
        }
        else if (property == "POPULATION") {
            // Read the population size
            input >> pop_size;                        
        }
        else if (property == "MIGRATE") {
            // Read the migrate flag
            input >> migrate;                         
        }
        else if (property == "NMIGRATIONS") {
            // Read the number of migrations
            input >> n_migrations;                    
        }
        else if (property == "MUTPROB") {
            // Read the mutation probability
            input >> mut_prob;                        
        }
        else if (property == "CROSSPROB") {
            // Read the crossover probability
            input >> cross_prob;                      
        }
    }
    input.close();                                   
    int migr_step = n_generations / n_migrations;     // Calculate migration step interval
    double exponent = 2.;                             // Set the exponent for selection

    // Initialize the population with labels, size, and rank
    population pop(labels, pop_size, rank);           
    vector<double> fitness(pop_size);                 
    ofstream ofitness("OUTPUT" + to_string(migrate) + "/fitness" + to_string(rank) + ".dat");
    ofitness << "       GEN" << "           L2" << endl;  // Write header to the fitness output file
    ofstream obestfitness("OUTPUT" + to_string(migrate) + "/best_fitness" + to_string(rank) + ".dat");
    obestfitness << "       GEN" << "           L2" << endl; // Write header to the best fitness output file

    for (int i = 0; i < n_generations; i++) {
        // Compute fitness of each individual in the population
        for (int j = 0; j < pop_size; j++) {
            individual ind = pop.get_individual(j);
            fitness[j] = compute_fitness(ind, cities);  // Calculate fitness for individual j
        }
        // Order the population according to fitness
        pop.order_population(fitness);
        // Sort fitness values (now the best fitness is in the first position)
        sort(fitness.begin(), fitness.end());
        // Print the best fitness value of the current generation
        obestfitness << setw(8) << i << setw(15) << fitness[0] << endl;
        // Compute the average fitness of the first half of the population
        double accu = 0.;
        for (int j = 0; j < pop_size / 2; j++) {
            accu += fitness[j];
        }
        accu /= double(pop_size) / 2;
    
        // Print the average fitness to the output file
        ofitness << setw(8) << i << setw(15) << accu << endl;

        // Check if it's time to migrate (exchange individuals between nodes)
        if (i % migr_step == 0 && size > 1 && migrate) {
            // Apply permutation on the vector of ranks
            permutation(vec_ranks, rnd);
            auto it_rank = find(vec_ranks.begin(), vec_ranks.end(), rank);
            int index_rank = distance(vec_ranks.begin(), it_rank);
        
            // Determine the ranks of nodes to send and receive data
            int idest = vec_ranks[pbc(index_rank + 1, vec_ranks.size())];
            int isource = vec_ranks[pbc(index_rank - 1, vec_ranks.size())];
            int itag_send = idest;
            int itag_receive = rank;
        
            // Unpack the information to exchange (genome of the best individual)
            individual best = pop.get_individual(0);
            vector<int> best_genome = best.get_genome();
            int genome_size = best_genome.size();
            
            // Prepare the data to send
            int* imesg = new int[genome_size];
            for (int j = 0; j < genome_size; j++) {
                imesg[j] = best_genome[j];
            }  
            // Prepare the data to receive
            int* imesg2 = new int[genome_size]; 
            
            MPI_Status status;
            // To avoid deadlocks, if rank is odd, first send then receive, and vice versa
            if (rank % 2) {
                MPI_Send(&imesg[0], genome_size, MPI_INTEGER, idest, itag_send, MPI_COMM_WORLD);
                MPI_Recv(&imesg2[0], genome_size, MPI_INTEGER, isource, itag_receive, MPI_COMM_WORLD, &status);
            } else {
                MPI_Recv(&imesg2[0], genome_size, MPI_INTEGER, isource, itag_receive, MPI_COMM_WORLD, &status);
                MPI_Send(&imesg[0], genome_size, MPI_INTEGER, idest, itag_send, MPI_COMM_WORLD);
            }  
            // Copy the received genome to the best_genome
            for (int j = 0; j < genome_size; j++) {
                best_genome[j] = imesg2[j];
            }
            best.set_genome(best_genome); 
            // Recompute the fitness of the first individual
            fitness[0] = compute_fitness(best, cities);
            // Set the new individual as the best one
            pop.set_individual(best, 0);
            // Reorder the population (since the immigrant could be worse than the current best individual)
            pop.order_population(fitness);
        }
        // Apply selection process
        pop.selection(exponent);
        // Apply crossover operations on the population
        pop.crossing(cross_prob);
        // Apply various mutation operations on the population
        pop.swap_mutation(mut_prob);
        pop.shift_mutation(mut_prob);
        pop.swap_group_mutation(mut_prob);
        pop.inversion_mutation(mut_prob); 
    }
    // Close the output files
    ofitness.close();
    obestfitness.close();


    // Now find the best individual among all
    ofstream ogenes("OUTPUT" + to_string(migrate) + "/genes" + to_string(rank) + ".dat");
    ogenes << "      GENE" << endl;  // Write header for the genes output file

    // Compute fitness for each individual in the population
    for (int j = 0; j < pop_size; j++) {
        individual ind = pop.get_individual(j);
        fitness[j] = compute_fitness(ind, cities);  // Calculate fitness for individual j
    }

    // Order the population according to fitness
    pop.order_population(fitness);
    // Get the best individual (the one with the highest fitness, i.e. shortest path)
    individual best = pop.get_individual(0);
    // Write the genes of the best individual to the output file
    for (int i = 0; i < best.get_genome_size(); i++) {
        ogenes << setw(8) << best.get_gene(i) << endl;
    }
    ogenes.close();  // Close the genes output file

    // Finalize the MPI environment
    MPI_Finalize();

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

void permutation(vector<int> ranks, Random rnd) {
    for (int i = 0; i < ranks.size(); i++) {
        int ind1 = int(rnd.Rannyu(0, ranks.size()));
        int ind2 = int(rnd.Rannyu(0, ranks.size()));
        swap(ranks[ind1], ranks[ind2]);
    }
}

int pbc(int i, int dim){
    if (i >= dim) { i = i - dim;} 
    else if (i < 0) { i = i + dim; }
    return i;
}