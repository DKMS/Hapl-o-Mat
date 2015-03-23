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
#include "Utility.h"
#include "Report.h"
#include "FisherInformation.h"

int main(int argc, char *argv[]){

  std::string format;
  if (argc == 2)
    format = argv[1];
  else{
    std::cerr << "Specify a file format (DKMS, GL, READ)" << std::endl;
    exit(EXIT_FAILURE);
  }

  std::cout << "#########Initialisation" << std::endl;
  timePoint t1;
  timePoint t2;
  double timeDataPreProcessing = 0.;
  double timeEMAlgorithm = 0.;
  double timeVariance = 0.;
  std::unique_ptr<Parameters> pParameters;
  std::unique_ptr<Data> pData;
  if(format == "DKMS"){
    std::unique_ptr<ParametersDKMS> pParametersTmp(new ParametersDKMS());
    std::unique_ptr<Data> pDataProcessingTmp(new DKMSDataProcessing(*pParametersTmp));
    pParameters = std::move(pParametersTmp);
    pData = std::move(pDataProcessingTmp);
  }
  else if(format == "GL"){
    std::unique_ptr<ParametersGL>pParametersTmp(new ParametersGL());
    std::unique_ptr<Data> pDataProcessingTmp(new GLDataProcessing(*pParametersTmp));
    pParameters = std::move(pParametersTmp);
    pData = std::move(pDataProcessingTmp);
  }
  else if(format == "READ"){
    std::unique_ptr<ParametersReadin>pParametersTmp(new ParametersReadin());
    std::unique_ptr<Data> pDataProcessingTmp(new DataReadin(*pParametersTmp));
    pParameters = std::move(pParametersTmp);
    pData = std::move(pDataProcessingTmp);
  }
  else{
    std::cerr << "Specify one of the file formats (DKMS, GL)" << std::endl;
    exit(EXIT_FAILURE);
  }

  std::cout << "#########Data-preprocessing" << std::endl;
  t1 = getTime();
  PhenotypeList pList;
  HaplotypeList hList(*pParameters);
  pData->dataProcessing(pList, hList);
  pData->printStatistics();
  std::cout << "Memory requirement haplotypes: " << hList.computeSizeInBytes() << " bytes" << std::endl;
  std::cout << "Memory requirement phenoypes: " << pList.computeSizeInBytes() << " bytes" << std::endl;
  t2 = getTime();
  timeDataPreProcessing = getTimeDifference(t1, t2);

  std::cout << "#########EM-algorithm" << std::endl;
  t1 = getTime();
  hList.initialiseFrequencies(pList);
  hList.EMAlgorithm(pList);
  t2 = getTime();
  timeEMAlgorithm = getTimeDifference(t1, t2);

  t1 = getTime();
  if(pParameters->getDoVariance()){
    std::cout << "#########Variance" << std::endl;
    pList.expectationAndRemoveStep(hList);
    fisherInformationParallel(hList, pList);
  }
  else{
    hList.writeFrequenciesToFile();
  }
  t2 = getTime();
  timeVariance = getTimeDifference(t1, t2);

  std::cout << "#########Time" << std::endl;
  std::cout << "\t Data pre-processing time: " << timeDataPreProcessing << " mus" << std::endl;
  std::cout << "\t EM-algorithm time: " << timeEMAlgorithm << " mus" << std::endl;
  std::cout << "\t Variance time: " << timeVariance << " mus" << std::endl;
}
