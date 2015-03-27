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

  std::string inputFileFormat;
  if (argc == 2)
    inputFileFormat = argv[1];
  else{
    std::cerr << "Specify a input file format (DKMS, GL, READ)" << std::endl;
    exit(EXIT_FAILURE);
  }

  std::cout << "#########Initialisation" << std::endl;
  timePoint startTime;
  timePoint endTime;
  double timeTakenForDataPreProcessing = 0.;
  double timeTakenForEMAlgorithm = 0.;
  double timeTakenForVariance = 0.;
  std::unique_ptr<Parameters> pParameters;
  std::unique_ptr<InputFile> pInputFile;
  if(inputFileFormat == "DKMS"){
    std::unique_ptr<ParametersDKMS> pParametersTmp(new ParametersDKMS());
    std::unique_ptr<InputFile> pInputFileTmp(new DKMS(*pParametersTmp));
    pParameters = std::move(pParametersTmp);
    pInputFile = std::move(pInputFileTmp);
  }
  else if(inputFileFormat == "GL"){
    std::unique_ptr<ParametersGL>pParametersTmp(new ParametersGL());
    std::unique_ptr<InputFile> pInputFileTmp(new GL(*pParametersTmp));
    pParameters = std::move(pParametersTmp);
    pInputFile = std::move(pInputFileTmp);
  }
  else if(inputFileFormat == "READ"){
    std::unique_ptr<ParametersReadin>pParametersTmp(new ParametersReadin());
    std::unique_ptr<InputFile> pInputFileTmp(new InputFileToRead(*pParametersTmp));
    pParameters = std::move(pParametersTmp);
    pInputFile = std::move(pInputFileTmp);
  }
  else{
    std::cerr << "Specify one of the input file formats (DKMS, GL, READ)" << std::endl;
    exit(EXIT_FAILURE);
  }

  std::cout << "#########Data-preprocessing" << std::endl;
  startTime = getTime();
  PhenotypeList phenotypes;
  HaplotypeList haplotypes(*pParameters);
  pInputFile->dataProcessing(phenotypes, haplotypes);
  pInputFile->printStatistics();
  std::cout << "\t Memory requirement haplotypes: " << haplotypes.computeSizeInBytes() << " bytes" << std::endl;
  std::cout << "\t Memory requirement phenoypes: " << phenotypes.computeSizeInBytes() << " bytes" << std::endl;
  endTime = getTime();
  timeTakenForDataPreProcessing = getTimeDifference(startTime, endTime);

  std::cout << "#########EM-algorithm" << std::endl;
  startTime = getTime();
  haplotypes.initialiseFrequencies(phenotypes);
  haplotypes.EMAlgorithm(phenotypes);
  endTime = getTime();
  timeTakenForEMAlgorithm = getTimeDifference(startTime, endTime);

  startTime = getTime();
  if(pParameters->getDoVariance()){
    std::cout << "#########Variance" << std::endl;
    double memory = haplotypes.getSize() * haplotypes.getSize() * 8. / 1024. / 1024.;
    if(memory >= MAX_MEMORY){
      std::cout << "\t Size of information matrix exceeds memory." << std::endl;
    }
    else{
      phenotypes.expectationAndRemoveStep(haplotypes);
      fisherInformation(haplotypes, phenotypes);
    }
  }
  else{
    haplotypes.writeFrequenciesToFile();
  }
  endTime = getTime();
  timeTakenForVariance = getTimeDifference(startTime, endTime);

  std::cout << "#########Time" << std::endl;
  std::cout << "\t Data pre-processing time: " << timeTakenForDataPreProcessing << " mus" << std::endl;
  std::cout << "\t EM-algorithm time: " << timeTakenForEMAlgorithm << " mus" << std::endl;
  std::cout << "\t Variance time: " << timeTakenForVariance << " mus" << std::endl;
}
