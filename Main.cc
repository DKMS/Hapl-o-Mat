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
  std::unique_ptr<DataProcessing> pDataProcessing;
  if(format == "DKMS"){
    std::unique_ptr<ParametersDKMS> pParametersTmp(new ParametersDKMS());
    std::unique_ptr<DataProcessing> pDataProcessingTmp(new DKMSDataProcessing(*pParametersTmp));
    pParameters = std::move(pParametersTmp);
    pDataProcessing = std::move(pDataProcessingTmp);
  }
  else if(format == "GL"){
    std::unique_ptr<ParametersGL>pParametersTmp(new ParametersGL());
    std::unique_ptr<DataProcessing> pDataProcessingTmp(new GLDataProcessing(*pParametersTmp));
    pParameters = std::move(pParametersTmp);
    pDataProcessing = std::move(pDataProcessingTmp);
  }
  else{
    std::cerr << "Specify one of the file formats (DKMS, GL)" << std::endl;
    exit(EXIT_FAILURE);
  }

  std::cout << "#########Data-preprocessing" << std::endl;
  PhenotypeList pList;
  HaplotypeList hList(*pParameters);
  pDataProcessing->dataProcessing(pList, hList);

  std::cout << "\t Number loci: " << pDataProcessing->getNumberLoci() << std::endl;
  std::cout << "\t Removed reports: " << pDataProcessing->getNumberRemovedDonors() << std::endl;
  std::cout << "\t Leftover Reports: " << pDataProcessing->getNumberDonors() << std::endl;
  std::cout << "\t Phenotypes: " << pList.getSize() << std::endl;
  std::cout << "\t Haplotypes: " << hList.getSize() << std::endl;
  std::cout << std::endl;

  std::cout << "#########EM-algorithm" << std::endl;

  hList.initialiseFrequencies(pList);
  hList.EMAlgorithm(pList);
  hList.writeFrequenciesToFile();
}
