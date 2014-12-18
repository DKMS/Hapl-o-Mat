#include <iostream>
#include <cmath>

#include "Parameters.h"
#include "Utility.h"

void Parameters::val_assign(size_t & out, const std::string line){
  size_t pos = line.find("=");
  std::string value = line.substr(pos + 1);
  out = std::stoi(value);
}

void Parameters::val_assign(double & out, const std::string line){
  size_t pos = line.find("=");
  std::string value = line.substr(pos + 1);
  out = std::stod(value);
}

void Parameters::val_assign(std::string & out, const std::string line){
  size_t pos = line.find("=");
  out = line.substr(pos + 1);
}

void Parameters::initType_assign(const std::string line){
  size_t pos = line.find("=");
  std::string value = line.substr(pos + 1);
  if(value.compare("random") == 0)
    initType = random;
  else if(value.compare("perturbation") == 0)
    initType = perturbation;
  else if(value.compare("numberOccurence") == 0)
    initType = numberOccurence;
  else{
    std::cerr << "No initialisation routine for haplotype frequencies specified. Set routine to random" << std::endl;
    initType = random;
  }
}

void Parameters::precision_assign(const std::string line){
  size_t pos = line.find("=");
  std::string value = line.substr(pos + 1);
  if(value == "g")
    precision = Allele::codePrecision::g;
  else if(value == "4d")
    precision = Allele::codePrecision::fourDigit;
  else if(value =="G")
    precision = Allele::codePrecision::G;
  else if(value == "6d")
    precision = Allele::codePrecision::sixDigit;
  else if(value == "8d")
    precision = Allele::codePrecision::eightDigit;
  else if(value == "asItIs")
    precision = Allele::codePrecision::asItIs;
  else{
    std::cerr << "No code precision specified" << std::endl;
    exit(EXIT_FAILURE);
  }
}

void Parameters::computePrintPrecision(){
  if(epsilon < 1)
    printPrecision = size_t (- log10(epsilon));
  else
    printPrecision = 0;
}

std::string Parameters::printInitialisationHaplotypeFrequencies() const{

  std::string out;
  switch(initType){
  case random:
    {
      out = "random";
      break;
    }
  case perturbation:
    {
      out = "perturbation";
      break;
    }
  case numberOccurence:
    {
      out = "numberOccurence";
      break;
    }
  }
  return out;
}

void ParametersGL::init(){

  std::ifstream file;
  openFileToRead(parametersFileName, file);

  std::string line;
  while(std::getline(file, line)){
    if(line.find("#") != std::string::npos) continue;
    else if(line.find("FILENAME_PULL") != std::string::npos) val_assign(pullFileName, line);
    else if(line.find("FILENAME_GLID") != std::string::npos) val_assign(glidFileName, line);
    else if(line.find("FILENAME_HAPLOTYPES") != std::string::npos) val_assign(haplotypesFileName, line);
    else if(line.find("FILENAME_PHENOTYPES") != std::string::npos) val_assign(phenotypesFileName, line);
    else if(line.find("FILENAME_HAPLOTYPEFREQUENCIES") != std::string::npos) val_assign(haplotypeFrequenciesFileName, line);
    else if(line.find("FILENAME_EPSILON") != std::string::npos) val_assign(epsilonFileName, line);
    else if(line.find("CODE_PRECISION") != std::string::npos) precision_assign(line);
    else if(line.find("MINIMAL_FREQUENCY_PHENOTYPES") != std::string::npos) val_assign(minimalFrequency, line);
    else if(line.find("LOCI") != std::string::npos) loci_assign(line);
    else if(line.find("INITIALISATION_HAPLOTYPE_FREQUENCIES") != std::string::npos) initType_assign(line);
    else if(line.find("EPSILON") != std::string::npos) val_assign(epsilon, line);
    else if(line.find("SEED") != std::string::npos) val_assign(seed, line);
    else{
      std::cerr << "Could not match "
		<< line
		<< " of parametersGL file."
		<< std::endl;
    }
  }//while
  computePrintPrecision();
  file.close();
}

void ParametersGL::loci_assign(const std::string line){

  size_t pos = line.find("=");
  std::string value = line.substr(pos + 1);
  lociToDo = split(value, ',');
}

void ParametersGL::print() const {

  std::cout << "GL format" << std::endl;
  std::cout << "#########Parameters I/O" << std::endl;
  std::cout << "\t Read data from pull-file: " << pullFileName << std::endl; 
  std::cout << "\t Read data from glid-file: " << glidFileName << std::endl; 
  std::cout << "\t Write haplotypes to: " << haplotypesFileName << std::endl;
  std::cout << "\t Write phenotypes to: " << phenotypesFileName << std::endl;
  std::cout << "\t Write estimated haplotype frequencies to: " << haplotypeFrequenciesFileName << std::endl;
  std::cout << "\t Write epsilon vs steps to: " << epsilonFileName << std::endl;
  std::cout << "#########Parameters resolving reports" << std::endl;
  std::cout << "\t Minimal frequency of phenotypes= " << minimalFrequency << std::endl;
  std::cout << "\t Resolve codes to precision: " << Allele::printCodePrecision(precision) << std::endl;
  std::cout << "\t Consider loci: ";
  for(auto locus : lociToDo){std::cout << locus << " ";}
  std::cout << std::endl;
  std::cout << "#########Parameters EM-algorithm" << std::endl;
  std::cout << "\t Initialisation haplotype frequencies: " << printInitialisationHaplotypeFrequencies() << std::endl;
  std::cout << "\t Epsilon= " << epsilon << std::endl;
  std::cout << "\t Zero= " << ZERO << std::endl;
  std::cout << "\t PrintPrecision= " << printPrecision << std::endl;
  std::cout << "\t Seed= " << seed << std::endl;
  std::cout << std::endl;
}

void ParametersDKMS::init(){

  std::ifstream file;
  openFileToRead(parametersFileName, file);

  std::string line;
  while(std::getline(file, line)){
    if(line.find("#") != std::string::npos) continue;
    else if(line.find("FILENAME_INPUT") != std::string::npos) val_assign(inputFileName, line);
    else if(line.find("FILENAME_HAPLOTYPES") != std::string::npos) val_assign(haplotypesFileName, line);
    else if(line.find("FILENAME_PHENOTYPES") != std::string::npos) val_assign(phenotypesFileName, line);
    else if(line.find("FILENAME_HAPLOTYPEFREQUENCIES") != std::string::npos) val_assign(haplotypeFrequenciesFileName, line);
    else if(line.find("FILENAME_EPSILON") != std::string::npos) val_assign(epsilonFileName, line);
    else if(line.find("CODE_PRECISION") != std::string::npos) precision_assign(line);
    else if(line.find("MINIMAL_FREQUENCY_PHENOTYPES") != std::string::npos) val_assign(minimalFrequency, line);
    else if(line.find("INITIALISATION_HAPLOTYPE_FREQUENCIES") != std::string::npos) initType_assign(line);
    else if(line.find("EPSILON") != std::string::npos) val_assign(epsilon, line);
    else if(line.find("SEED") != std::string::npos) val_assign(seed, line);
    else{
      std::cerr << "Could not match "
		<< line
		<< " of parametersDKMS file."
		<< std::endl;
    }
  }//while
  computePrintPrecision();
  file.close();
}

void ParametersDKMS::print() const {

  std::cout << "DKMS format" << std::endl;
  std::cout << "#########Parameters I/O" << std::endl;
  std::cout << "\t Read data from: " << inputFileName << std::endl;
  std::cout << "\t Write haplotypes to: " << haplotypesFileName << std::endl;
  std::cout << "\t Write phenotypes to: " << phenotypesFileName << std::endl;
  std::cout << "\t Write estimated haplotype frequencies to: " << haplotypeFrequenciesFileName << std::endl;
  std::cout << "\t Write epsilon vs steps to: " << epsilonFileName << std::endl;
  std::cout << "#########Parameters resolving reports" << std::endl;
  std::cout << "\t Minimal frequency of phenotypes= " << minimalFrequency << std::endl;
  std::cout << "\t Resolve codes to precision: " << Allele::printCodePrecision(precision) << std::endl;
  std::cout << "#########Parameters EM-algorithm" << std::endl;
  std::cout << "\t Initialisation haplotype frequencies: " << printInitialisationHaplotypeFrequencies() << std::endl;
  std::cout << "\t Epsilon= " << epsilon << std::endl;
  std::cout << "\t Zero= " << ZERO << std::endl;
  std::cout << "\t PrintPrecision= " << printPrecision << std::endl;
  std::cout << "\t Seed= " << seed << std::endl;
  std::cout << std::endl;
}
