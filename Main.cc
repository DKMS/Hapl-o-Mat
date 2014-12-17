#include <string>
#include <iostream>
#include <memory>

#include "DataProcessing.h"
#include "Allele.h"
#include "Locus.h"
#include "Typedefs.h"
#include "Phenotype.h"
#include "Haplotype.h"
#include "Parameters.h"

int main(int argc, char *argv[]){

  std::string format;
  if (argc == 2)
    format = argv[1];
  else{
    std::cerr << "Specify a file format (DKMS, GL)" << std::endl;
    exit(EXIT_FAILURE);
  }

  std::unique_ptr<Parameters> pParameters;
  if(format == "DKMS"){
    std::unique_ptr<Parameters> pParametersTmp(new ParametersDKMS());
    pParameters = move(pParametersTmp);
  }
  else if(format == "GL"){
    std::unique_ptr<Parameters> pParametersTmp(new ParametersGL());
    pParameters = move(pParametersTmp);
  }
  else{
    std::cerr << "Specify one of the file formats (DKMS, GL)" << std::endl;
    exit(EXIT_FAILURE);
  }
  
  strVec_t lociToDo;
  lociToDo.push_back("A");
  lociToDo.push_back("C");
  lociToDo.push_back("None");
  lociToDo.push_back("DRB1");
  lociToDo.push_back("DQB1");

  //  std::unique_ptr<DataProcessing> pDataProcessing (new GLDataProcessing("reports.pull", "results/haplotypes.dat", "results/phenotypes.dat", "reports.glid", lociToDo,   Allele::codePrecision::fourDigit));
  std::unique_ptr<DataProcessing> pDataProcessing (new DKMSDataProcessing("reports.txt",
									  "results/haplotypes.dat",
									  "results/phenotypes.dat",
									  Allele::codePrecision::g,
									  .0001));
  PhenotypeList pList;
  HaplotypeList hList;
  pDataProcessing->dataProcessing(pList, hList);

  std::cout << "\t Removed reports: " << pDataProcessing->getNumberRemovedDonors() << std::endl;
  std::cout << "\t Leftover Reports: " << pDataProcessing->getNumberDonors() << std::endl;
  std::cout << "\t Phenotypes: " << pList.getSize() << std::endl;
  std::cout << "\t Haplotypes: " << hList.getSize() << std::endl;
  std::cout << std::endl;
  std::cout << "EM-algorithm" << std::endl;

  /*
  hList.initialiseFrequencies(pList);
  EMAlgorithm(pList, hList);
  hList.writeFrequenciesToFile();
  */

}
