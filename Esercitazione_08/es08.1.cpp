#include <iostream>
#include <cmath>
#include "system.h"

using namespace std;

int main (int argc, char *argv[]){

  System SYS;
  SYS.block_reset(0);

  for(int i=0; i < SYS.get_nbl(); i++){ //loop over blocks
    for(int j=0; j < SYS.get_nsteps(); j++){ //loop over steps in a block
      SYS.step();
      SYS.measure();
    }
    SYS.average(i+1);
    SYS.block_reset(i+1);
  }
  SYS.finalize();

  return 0;
}