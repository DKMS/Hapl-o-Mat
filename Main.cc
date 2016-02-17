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
  double timeTakenForWriting = 0.;
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
  Phenotypes phenotypes;
  Haplotypes haplotypes(*pParameters);
  pInputFile->dataProcessing(phenotypes, haplotypes);
  pInputFile->printStatistics();
  std::cout << "\t Memory requirement haplotypes: " << haplotypes.computeSizeInBytes() << " bytes" << std::endl;
  std::cout << "\t Memory requirement phenoypes: " << phenotypes.computeSizeInBytes() << " bytes" << std::endl;
  endTime = getTime();
  timeTakenForDataPreProcessing = getTimeDifference(startTime, endTime);

  std::cout << "#########EM-algorithm" << std::endl;
  startTime = getTime();
  double minEpsilon = .5 / static_cast<double>(haplotypes.getNumberDonors());
  if(minEpsilon - haplotypes.getEpsilon() < ZERO){
    std::cerr << "Chosen epsilon is larger than 0.5/n" <<std::endl;
    exit(EXIT_FAILURE);
  }
  else{
    haplotypes.initialiseFrequencies(phenotypes);
    haplotypes.EMAlgorithm(phenotypes);
    std::cout << "\t Summed haplotype frequencies: " << haplotypes.computeHaplotypeFrequencySum() << std::endl;
    std::cout << "\t Summed cutted haplotype frequencies: " << haplotypes.computeCuttedHaplotypeFrequencySum() << std::endl;
  }
  endTime = getTime();
  timeTakenForEMAlgorithm = getTimeDifference(startTime, endTime);

  startTime = getTime();
  haplotypes.writeFrequenciesToFile();
  haplotypes.deleteHaplotypesFile();
  endTime = getTime();
  timeTakenForWriting = getTimeDifference(startTime, endTime);

  std::cout << "#########Time" << std::endl;
  std::cout << "\t Data pre-processing time: " << timeTakenForDataPreProcessing << " mus" << std::endl;
  std::cout << "\t EM-algorithm time: " << timeTakenForEMAlgorithm << " mus" << std::endl;
  std::cout << "\t Writing time: " << timeTakenForWriting << " mus" << std::endl;
}
