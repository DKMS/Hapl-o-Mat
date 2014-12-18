#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream> 
#include <algorithm>

#include "Haplotype.h"
#include "Phenotype.h"
#include "Utility.h"


/*
void EMAlgorithm(PhenotypeList & phenotypes, HaplotypeList & haplotypes){

  std::fstream file;
  file.open(parameters.getEpsilonFileName(), std::ifstream::out);
  if(!file.is_open()) {
    std::cout << "Could not open file: "
              << parameters.getEpsilonFileName()
              << std::endl;
  }
  file.precision(14);
  file << std::scientific;

  bool stop;
  size_t counter = 0;
  do{
    counter ++;
    stop = false;

    phenotypes.expectationStep(haplotypes);
    double largestEpsilon = 0.;
    haplotypes.maximizationStep(phenotypes, largestEpsilon);
    
    if(largestEpsilon - parameters.getEpsilon() < ZERO ){
      stop = true;
    }
    file << largestEpsilon << std::endl;
  } while(!stop);
  std::cout << "Used " << counter <<" steps" << std::endl;
  file.close();
}
*/

void HaplotypeList::initialiseFrequencies(const PhenotypeList & phenotypes){
  
  switch(initType){
  case Parameters::initialisationHaplotypeFrequencies::numberOccurence:
    {
      initialiseNumberOccurence(phenotypes);
      break;
    }
  case Parameters::initialisationHaplotypeFrequencies::perturbation:
    {
      initialiseNumberOccurence(phenotypes);
      initialisePerturbation();
      break;
    }
  case Parameters::random:
    {
      double normalisation = 0.;
      for(auto it = listBegin();
	  it != listEnd();
	  it ++){
	double number = rng() + 1.;
	normalisation += number;
	it->second.setFrequency(number);
      }
      for(auto it = listBegin();
	  it != listEnd();
	  it ++){
	double newFreq = it->second.getFrequency()/normalisation;
	it->second.setFrequency(newFreq);
      }
      break;
    }
  }//switch
}

void HaplotypeList::initialiseNumberOccurence(const PhenotypeList & phenotypes){

  double factor = static_cast<double>(pow(2, numberLoci));
  factor *= static_cast<double>(numberDonors);

  auto itPhenoEnd = phenotypes.c_listEnd();
  for(auto itPheno = phenotypes.c_listBegin();
      itPheno != itPhenoEnd;
      itPheno++)
    {
      double numberPhenotype = itPheno->second.getNumInDonors() / factor;
      auto itDiploEnd = itPheno->second.c_diplotypeListEnd();
      for(auto itDiplo = itPheno->second.c_diplotypeListBegin();
          itDiplo != itDiploEnd;
          itDiplo ++)
        {
	  auto itHaplo1 = hashList.find(itDiplo->id1);
	  auto itHaplo2 = hashList.find(itDiplo->id2);
	  itHaplo1->second.addFrequency(numberPhenotype);
	  itHaplo2->second.addFrequency(numberPhenotype);
	}
    }
}

void HaplotypeList::initialisePerturbation(){

  std::uniform_real_distribution<> dis(0, 1);
  double sum = 0.;
  for(auto it = listBegin();
      it != listEnd();
      it ++){
    double oldFreq = it->second.getFrequency();
    double newFreq = oldFreq + dis(rng) * oldFreq;

    it->second.addFrequency(newFreq);
    sum += newFreq;
  }

  double normalisation = 1./sum;
  for(auto it = listBegin();
      it != listEnd();
      it ++){
    it->second.multiplyFrequency(normalisation);
  }
}

void HaplotypeList::writeFrequenciesToFile() const{

  std::ifstream inFile;
  openFileToRead(haplotypesFileName ,inFile);
  std::ofstream outFile;
  openFileToWrite(haplotypeFrequenciesFileName ,outFile);
  outFile.precision(14);
  outFile << std::fixed;

  std::string code;
  while(inFile >> code){
    size_t hashValue = string_hash(code);
    auto pos = hashList.find(hashValue);
    if(pos != hashList.end()){
      double freq = pos->second.getFrequency();

      //      if(parameters.getMonteCarlo()){
      //	fileOut << code
      //		<< "\t";
      //	fileOut << freq
      //		<< "\n";
      //}//if Monte Carlo
    //      else{
      if(freq > epsilon){
	outFile << code
		<< "\t";
	outFile << freq
		<< "\n";
      }
      //      }//else Montecarlo
    }
  }
  inFile.close();
  outFile.close();
}

/*
void HaplotypeList::maximizationStep(const PhenotypeList & phenotypes, double & largestEpsilon){

    //save old frequencies, initialize haplotype frequencies to zero
    std::vector<double> oldFrequencies;
    oldFrequencies.reserve(hashList.size());
    auto itEnd = hashList.end();
    for(auto it = hashList.begin();
	it != itEnd;
	it++){
      oldFrequencies.push_back(it->second.getFrequency());
      it->second.setFrequency(0.);
    }

    maximization(phenotypes);

    auto itOld = oldFrequencies.cbegin();
    for(auto it = hashList.begin();
	it != hashList.end();
	it++, itOld ++){
      it->second.multiplyFrequency(.5 / static_cast<double>(parameters.getNumberDonors()));
      
      double possibleEpsilon = fabs(it->second.getFrequency() - *itOld);
      if(possibleEpsilon > largestEpsilon){
	largestEpsilon = possibleEpsilon; 
      }
    }

    for(auto it = hashList.begin();
	it != hashList.end();
	it++){
      if(it->second.getFrequency() < ZERO){
	hashList.erase(it);
      }
    }
}

void HaplotypeList::maximization(const PhenotypeList & phenotypes){

  auto itPhenoEnd = phenotypes.c_listEnd();
  for(auto itPheno = phenotypes.c_listBegin();
      itPheno != itPhenoEnd;
      ++ itPheno){

    double factor = itPheno->second.getNumInDonors() / itPheno->second.computeSummedFrequencyDiplotypes();
		
    auto itDiploEnd = itPheno->second.c_diplotypeListEnd();
    for(auto itDiplo = itPheno->second.c_diplotypeListBegin();
	itDiplo != itDiploEnd;
	++ itDiplo)
      { 
	double freq = factor * itDiplo->frequency;
	auto itHaplo1 = hashList.find(itDiplo->id1);
	if(itHaplo1 != hashList.end()){
	  itHaplo1->second.addFrequency(freq);
	}
	if(itDiplo->id1 == itDiplo->id2){
	  itHaplo1->second.addFrequency(freq);
	}
	else{
	  auto itHaplo2 = hashList.find(itDiplo->id2);
	  if(itHaplo2 != hashList.end()){
	    itHaplo2->second.addFrequency(freq);
	  }
	}
      }//for diplo
  }
}

*/
