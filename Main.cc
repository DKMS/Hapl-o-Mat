/*
 * Hapl-O-mat: A program for HLA haplotype frequency estimation
 *
 * Copyright (C) 2016, DKMS gGmbH 
 *
 * This file is part of Hapl-O-mat
 *
 * Hapl-O-mat is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * Hapl-O-mat is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Hapl-O-mat; see the file COPYING.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

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
    std::cerr << "Specify a input file format (MA, GL, GLC or READ)" << std::endl;
    exit(EXIT_FAILURE);
  }

  std::cout << std::endl;
  std::cout << "\tHapl-O-mat" << std::endl;
  std::cout << "\tCopyright (C) 2016 DKMS gGmbH" << std::endl;
  std::cout << std::endl;

  std::cout << "#########Initialization" << std::endl;
  timePoint startTime;
  timePoint endTime;
  double timeTakenForDataPreProcessing = 0.;
  double timeTakenForEMAlgorithm = 0.;
  double timeTakenForWriting = 0.;
  std::unique_ptr<Parameters> pParameters;
  std::unique_ptr<InputFile> pInputFile;
  if(inputFileFormat == "MA"){
    std::unique_ptr<ParametersMA> pParametersTmp = make_unique<ParametersMA>();
    std::unique_ptr<InputFile> pInputFileTmp = make_unique<MA>(*pParametersTmp);
    pParameters = std::move(pParametersTmp);
    pInputFile = std::move(pInputFileTmp);
  }
  else if(inputFileFormat == "GL"){
    std::unique_ptr<ParametersGL>pParametersTmp = make_unique<ParametersGL>();
    std::unique_ptr<InputFile> pInputFileTmp = make_unique<GL>(*pParametersTmp);
    pParameters = std::move(pParametersTmp);
    pInputFile = std::move(pInputFileTmp);
  }
  else if(inputFileFormat == "GLC"){
    std::unique_ptr<ParametersGLC>pParametersTmp = make_unique<ParametersGLC>();
    std::unique_ptr<InputFile> pInputFileTmp = make_unique<GLC>(*pParametersTmp);
    pParameters = std::move(pParametersTmp);
    pInputFile = std::move(pInputFileTmp);
  }
  else if(inputFileFormat == "READ"){
    std::unique_ptr<ParametersReadin>pParametersTmp = make_unique<ParametersReadin>();
    std::unique_ptr<InputFile> pInputFileTmp = make_unique<InputFileToRead>(*pParametersTmp);
    pParameters = std::move(pParametersTmp);
    pInputFile = std::move(pInputFileTmp);
  }
  else{
    std::cerr << "Specify one of the input file formats (MA, GL, READ)" << std::endl;
    exit(EXIT_FAILURE);
  }

  startTime = getTime();
  Phenotypes phenotypes;
  Haplotypes haplotypes(*pParameters);
  pInputFile->dataProcessing(phenotypes, haplotypes);
  pInputFile->printStatistics();
  std::cout << "\t Memory requirement haplotypes: " << haplotypes.computeSizeInBytes() << " bytes" << std::endl;
  std::cout << "\t Memory requirement genotypes: " << phenotypes.computeSizeInBytes() << " bytes" << std::endl;
  endTime = getTime();
  timeTakenForDataPreProcessing = getTimeDifference(startTime, endTime)/1000000.;

  std::cout << "#########EM algorithm" << std::endl;
  startTime = getTime();
  double minEpsilon = .5 / static_cast<double>(haplotypes.getNumberDonors());
  if(haplotypes.getEpsilon() - minEpsilon > ZERO){
    std::cout << "Chosen epsilon is larger than 0.5/n" <<std::endl;
  }
  else{
    haplotypes.initialiseFrequencies(phenotypes);
    haplotypes.EMAlgorithm(phenotypes);
    std::cout << "\t Sum haplotype frequencies: " << haplotypes.computeHaplotypeFrequencySum() << std::endl;
    std::cout << "\t Sum cutted haplotype frequencies: " << haplotypes.computeCuttedHaplotypeFrequencySum() << std::endl;
  }
  endTime = getTime();
  timeTakenForEMAlgorithm = getTimeDifference(startTime, endTime)/1000000.;

  startTime = getTime();
  haplotypes.writeFrequenciesToFile();
  haplotypes.deleteHaplotypesFile();
  endTime = getTime();
  timeTakenForWriting = getTimeDifference(startTime, endTime)/1000000.;

  std::cout << "#########Time" << std::endl;
  std::cout << "\t Data preprocessing time [s]: " << timeTakenForDataPreProcessing << std::endl;
  std::cout << "\t EM algorithm time [s]: " << timeTakenForEMAlgorithm << std::endl;
  std::cout << "\t Writing time [s]: " << timeTakenForWriting << std::endl;
}
