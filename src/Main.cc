/*
 * Hapl-o-Mat: A software for haplotype inference
 *
 * Copyright (C) 2016, DKMS gGmbH 
 *
 * Dr. Jürgen Sauter
 * Kressbach 1
 * 72072 Tuebingen, Germany
 *
 * T +49 7071 943-2060
 * F +49 7071 943-2090
 * sauter(at)dkms.de
 *
 * This file is part of Hapl-o-Mat
 *
 * Hapl-o-Mat is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * Hapl-o-Mat is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Hapl-o-Mat; see the file COPYING.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <memory>
#include <string>

#include "DataProcessing.h"
#include "Exceptions.h"
#include "Haplotype.h"
#include "Parameters.h"
#include "Phenotype.h"
#include "Utility.h"

int main(int argc, char *argv[]){

  std::cout << std::endl;
  std::cout << "\t Hapl-o-Mat" << std::endl;
  std::cout << "\t Copyright (C) 2016 DKMS gGmbH" << std::endl;
  std::cout << std::endl;

  try{
    std::cout << "#########Initialization" << std::endl;
    std::unique_ptr<Parameters> pParameters;
    std::unique_ptr<InputFile> pInputFile;

    std::string inputFileFormat;
    if(argc < 2)
      {
	throw InputFormatException();
      }
    inputFileFormat = argv[1];

    if(inputFileFormat == "MAC"){
      std::unique_ptr<ParametersMAC> pParametersTmp = std::make_unique<ParametersMAC>();
      std::unique_ptr<InputFile> pInputFileTmp = std::make_unique<MAC>(*pParametersTmp);
      pParameters = std::move(pParametersTmp);
      pInputFile = std::move(pInputFileTmp);
    }
    else if(inputFileFormat == "GLS"){
      std::unique_ptr<ParametersGLS>pParametersTmp = std::make_unique<ParametersGLS>();
      std::unique_ptr<InputFile> pInputFileTmp = std::make_unique<GLS>(*pParametersTmp);
      pParameters = std::move(pParametersTmp);
      pInputFile = std::move(pInputFileTmp);
    }
    else if(inputFileFormat == "GLSC"){
      std::unique_ptr<ParametersGLSC>pParametersTmp = std::make_unique<ParametersGLSC>();
      std::unique_ptr<InputFile> pInputFileTmp = std::make_unique<GLSC>(*pParametersTmp);
      pParameters = std::move(pParametersTmp);
      pInputFile = std::move(pInputFileTmp);
    }
    else if(inputFileFormat == "READ"){
      std::unique_ptr<ParametersReadin>pParametersTmp = std::make_unique<ParametersReadin>();
      std::unique_ptr<InputFile> pInputFileTmp = std::make_unique<InputFileToRead>(*pParametersTmp);
      pParameters = std::move(pParametersTmp);
    pInputFile = std::move(pInputFileTmp);
    }
    else{
      throw InputFormatException();
    }

    timePoint startTime = getTime();
    Phenotypes phenotypes;
    Haplotypes haplotypes(*pParameters);
    pInputFile->dataProcessing(phenotypes, haplotypes);
    pInputFile->printStatistics();
    std::cout << "\t Memory requirement haplotypes [MB]: " << haplotypes.computeSizeInBytes()/1000./1000. << std::endl;
    std::cout << "\t Memory requirement genotypes [MB]: " << phenotypes.computeSizeInBytes()/1000./1000. << std::endl;
    timePoint endTime = getTime();
    double timeTakenForDataPreProcessing = getTimeDifference(startTime, endTime)/1000000.;
    
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
    double timeTakenForEMAlgorithm = getTimeDifference(startTime, endTime)/1000000.;
    
    startTime = getTime();
    haplotypes.writeFrequenciesToFile();
    haplotypes.deleteHaplotypesFile();
    endTime = getTime();
    double timeTakenForWriting = getTimeDifference(startTime, endTime)/1000000.;
    
    std::cout << "#########Times" << std::endl;
    std::cout.precision(14);
    std::cout << "\t Data preprocessing [s]: " << timeTakenForDataPreProcessing << std::endl;
    std::cout << "\t EM algorithm [s]: " << timeTakenForEMAlgorithm << std::endl;
    std::cout << "\t Writing [s]: " << timeTakenForWriting << std::endl;

    return 0;
  }
  catch(const std::exception & e){
    std::cerr << e.what() << std::endl;
    std::cout << "Exit Hapl-o-Mat" << std::endl;
    return -1;
  }
}
