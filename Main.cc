#include <string>
#include <iostream>
#include <memory>
#include <chrono>

#include "DataProcessing.h"
#include "Allele.h"
#include "Locus.h"
#include "Typedefs.h"
#include "Phenotype.h"
#include "Haplotype.h"
#include "Parameters.h"
#include "Utility.h"
#include "Report.h"

int main(int argc, char *argv[]){

  std::string format;
  if (argc == 2)
    format = argv[1];
  else{
    std::cerr << "Specify a file format (DKMS, GL)" << std::endl;
    exit(EXIT_FAILURE);
  }

  std::cout << "#########Initialisation" << std::endl;
  timePoint t1;
  timePoint t2;
  double timeDataPreProcessing = 0.;
  double timeEMAlgorithm = 0.;
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
  t1 = getTime();
  PhenotypeList pList;
  HaplotypeList hList(*pParameters);
  pDataProcessing->dataProcessing(pList, hList);
  std::cout << "\t Number loci: " << pDataProcessing->getNumberLoci() << std::endl;
  std::cout << "\t Removed reports: " << pDataProcessing->getNumberRemovedDonors() << std::endl;
  std::cout << "\t Leftover Reports: " << pDataProcessing->getNumberDonors() << std::endl;
  std::cout << "\t H0 loci/reports: " <<  Report::getNumberH0Loci() << "/" << Report::getNumberH0Reports() << std::endl;
  std::cout << "\t H1 loci/reports: " <<  Report::getNumberH1Loci() << "/" << Report::getNumberH1Reports() << std::endl;
  std::cout << "\t H2 loci/reports: " <<  Report::getNumberH2Loci() << "/" << Report::getNumberH2Reports() << std::endl;
  std::cout << "\t H2M loci/reports: " <<  Report::getNumberH2MLoci() << "/" << Report::getNumberH2MReports() << std::endl;
  std::cout << "\t I loci/reports: " <<  Report::getNumberILoci() <<  "/" << Report::getNumberIReports() <<std::endl;
  std::cout << "\t Phenotypes: " << pList.getSize() << std::endl;
  std::cout << "\t Haplotypes: " << hList.getSize() << std::endl;
  std::cout << std::endl;
  t2 = getTime();
  timeDataPreProcessing = getTimeDifference(t1, t2);

  std::cout << "#########EM-algorithm" << std::endl;
  t1 = getTime();
  hList.initialiseFrequencies(pList);
  hList.EMAlgorithm(pList);
  hList.writeFrequenciesToFile();
  t2 = getTime();
  timeEMAlgorithm = getTimeDifference(t1, t2);

  std::cout << "#########Time" << std::endl;
  std::cout << "\t Data pre-processing time: " << timeDataPreProcessing << " mus" << std::endl;
  std::cout << "\t EM-algorithm time: " << timeEMAlgorithm << " mus" << std::endl;
}
