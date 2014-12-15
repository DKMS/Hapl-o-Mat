#include <iostream>

#include "Locus.h"

void PhasedLocus::resolve(){

  for(auto it : phasedLocus){
    for(auto it2 : it){
      std::cout << it2 << std::endl;
    }
    std::cout << std::endl;
  }
}

void UnphasedLocus::resolve(){

  for(auto it : unphasedLocus){
    for(auto it2 : it){
      std::cout << it2 << std::endl;
    }
    std::cout << std::endl;
  }

}
