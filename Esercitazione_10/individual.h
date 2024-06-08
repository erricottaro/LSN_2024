#include "city.h"
#include "random.h"
#include <vector>
#include <algorithm>
#include <fstream>

using namespace std;

#ifndef __individual__
#define __individual__

class individual {
private:
    int _genome_size;
    vector<int> _genes;
public:
    individual(); // Default constructor
    ~individual(); // Destructor
        
    int get_gene(int i); // Method to get the gene (city) at position i
    vector<int> get_genome(); // Method to get the entire genome (sequence of cities)
    int get_genome_size(); // Method to get the size of the genome

    void set_gene(int new_city, int i);  // Method to set the gene (city) at position i with a new city
    void set_genome(vector<int> &input_cities);  // Method to set the entire genome with a given sequence of cities
};

#endif