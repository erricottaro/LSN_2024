#include <iostream>
#include <cmath>

using namespace std;

#ifndef __city__
#define __city__

class city {
private:
    double _x;
    double _y;
public:
    city();
    city(double x, double y);
    ~city();

    double get_coord(int i);

    void set_coord(double new_coord, int i);
};

#endif