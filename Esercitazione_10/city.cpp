#include "city.h"

using namespace std;

city::city(){}

city::city(double x, double y){
    _x = x;
    _y = y;
}

city::~city(){}

double city::get_coord(int i){
    if (i==0){
        return _x;
    }
    else if (i==1){
        return _y;
    }
    else {
        cerr << "ERROR: invalid index" << endl; 
        exit(-1);
    }
}

void city::set_coord(double new_coord, int i){
    if (i==0){
        _x=new_coord;
    }
    else if (i==1){
        _y=new_coord;
    }
    else {cerr << "ERROR: invalid index" << endl; }
}