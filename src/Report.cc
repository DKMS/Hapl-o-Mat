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
#include <sstream>

#include "DataProcessing.h"
#include "Exceptions.h"
#include "Haplotype.h"
#include "Locus.h"
#include "Phenotype.h"
#include "Report.h"
#include "Utility.h"

std::unordered_map<std::string, std::shared_ptr<Locus>> ColumnReport::singleLocusGenotypesAlreadyDone;
double Report::numberNReports = 0.;
double Report::numberAReports = 0.;
double Report::numberMReports = 0.;
double Report::numberIReports = 0.;

std::string BasicReport::buildPhenotypeCode() const{

  std::string phenotypeCode = "";
  for(auto genotypeAtLocus : genotypeAtLoci){
    phenotypeCode += genotypeAtLocus.at(0);
    phenotypeCode += "+";
    phenotypeCode += genotypeAtLocus.at(1);
    phenotypeCode += "^";
  }
  phenotypeCode.pop_back();

  return phenotypeCode;
}

void BasicReport::buildHaploAndDiplotypes(Phenotypes::iterator itPhenotype,
				     Haplotypes & haplotypes,
				     std::ofstream & haplotypesFile,
				     const HaplotypeCombinations & haplotypeCombinations) const{

  auto i1end = haplotypeCombinations.getList().cend();
  for(auto i1 = haplotypeCombinations.getList().cbegin();
      i1 != i1end;
      i1++)
    {
      std::string codeHaplotype1;
      std::string codeHaplotype2;
      
      auto genotypeAtLocus = genotypeAtLoci.cbegin();
     
      auto i2end = i1->cend();
      for(auto i2 = i1->cbegin();
          i2 != i2end;
          i2++)
        {
          if(*i2){
            codeHaplotype2.append(genotypeAtLocus->at(0));
            codeHaplotype1.append(genotypeAtLocus->at(1));
          }
          else{
	    codeHaplotype1.append(genotypeAtLocus->at(0));
	    codeHaplotype2.append(genotypeAtLocus->at(1));
          }
	  codeHaplotype1.append("~");
          codeHaplotype2.append("~");
	  genotypeAtLocus ++;
        }
      
      codeHaplotype1.pop_back();
      codeHaplotype2.pop_back();
      
      //add haplotypes to list                                                                                                                     
      std::pair<Haplotypes::iterator, bool> inserted1 = haplotypes.add(codeHaplotype1);
      std::pair<Haplotypes::iterator, bool> inserted2 = haplotypes.add(codeHaplotype2);
      
      if(inserted1.second){
	haplotypesFile << codeHaplotype1 << std::endl;
      }
      if(inserted2.second){
	haplotypesFile << codeHaplotype2 << std::endl;
      }
      
      //build diplotype                                                                                                                          
      size_t id1 = inserted1.first->first;
      size_t id2 = inserted2.first->first;
      
      Diplotype diplotype;
      diplotype.id1 = id1;
      diplotype.id2 = id2;
      
      itPhenotype->second.addDiplotype(diplotype);
    }//haplotypeCombinations  
}

void ReadinReport::translateLine(const std::string line){

  std::stringstream ss(line);
  std::string type;
  std::string code;
  std::string frequencyAsText;
  if(ss >> id >> type >> frequencyAsText >> code){}

  frequency = std::stod(frequencyAsText);
  strVec_t genotypes = split(code, '^');
  for(auto genotype : genotypes){
    strVec_t alleles = split(genotype, '+');
    strArr_t genotypeAtLocus;
    genotypeAtLocus.at(0) = alleles.at(0);
    genotypeAtLocus.at(1) = alleles.at(1);
    genotypeAtLoci.push_back(genotypeAtLocus);
  }
}

void ColumnReport::resolveSingleLocusGenotype(const std::unique_ptr<Genotype> & genotype, 
					      const size_t positionWantedLocus){
					  
  std::vector<std::pair<strArr_t, double>> genotypesAtLocus;

  auto pos = singleLocusGenotypesAlreadyDone.find(genotype->getSingleLocusGenotype());
  if(pos == singleLocusGenotypesAlreadyDone.cend())
    {
      std::shared_ptr<Locus> pLocus = genotype->resolve(doAmbiguityFilter, expandAmbiguityLines);
    
      types.at(positionWantedLocus) = pLocus->getType();
      pLocus->reduce(genotypesAtLocus);
      
      singleLocusGenotypesAlreadyDone.emplace(genotype->getSingleLocusGenotype(), pLocus);
    }
  else
    {
      types.at(positionWantedLocus) = pos->second->getType();
      pos->second->reduce(genotypesAtLocus);
    }
  
  numberOfReports *= static_cast<double>(genotypesAtLocus.size());
  if(1./numberOfReports - minimalFrequency >= ZERO){
    genotypesWithFrequenciesAtLoci.at(positionWantedLocus) = genotypesAtLocus; 
  }
  else
    {
      throw SplittingGenotypeException(); 
    }
}

void Report::buildListOfReports(std::vector<std::shared_ptr<Report>> & listOfReports){

  std::vector<std::vector<std::pair<strArr_t, double>>> reports;
  cartesianProduct(reports, genotypesWithFrequenciesAtLoci);
  
  for(auto report : reports){
    strArrVec_t newGenotypeAtLoci;
    double newFrequency = 1.;
    for(auto locus : report){
      newGenotypeAtLoci.push_back(locus.first);
      newFrequency *= locus.second;
    }
    listOfReports.push_back(this->create(newGenotypeAtLoci, newFrequency, numberLoci, id, types));
  }//reports
}

std::string Report::evaluateReportType(const size_t numberReports) const{
  
  if(find(types.cbegin(),
	  types.cend(),
	  Locus::reportType::I) != types.cend())
    numberIReports += 1./static_cast<double>(numberReports);
  else if(find(types.cbegin(),
	       types.cend(),
	       Locus::reportType::M) != types.cend())
    numberMReports += 1./static_cast<double>(numberReports);
  else if(find(types.cbegin(),
	       types.cend(),
	       Locus::reportType::A) != types.cend())
    numberAReports += 1./static_cast<double>(numberReports);
  else if(find(types.cbegin(),
	       types.cend(),
	       Locus::reportType::N) != types.cend())
    numberNReports += 1./static_cast<double>(numberReports);

  std::string totalType = "";
  for(auto type : types){
    switch(type){
    case Locus::reportType::N:
      {
	totalType += "N";
	break;
      }
    case Locus::reportType::A:
      {	totalType += "A";
	break;
      }
    case Locus::reportType::M:
      {
	totalType += "M";
	break;
      }
    case Locus::reportType::I:
      {
	totalType += "I";
	break;
      }
    default:
      totalType += "?";
    }
  }
  return totalType;
}

void GLSReport::translateLine(const std::string line){

  id = leftOfFirstDelim(line, ';');

  std::string rightPartOfLine = rightOfFirstDelim(line, ';');
  strVec_t allGlids = split(rightPartOfLine, ':');

  glids.resize(numberLoci, 0);

  auto locusName = lociOrder.cbegin();
  for(auto glidNumber : allGlids)
    {
      auto locusAndResolution = lociAndResolutions.find(*locusName);
      if(locusAndResolution != lociAndResolutions.cend())
	{
	  size_t number = stoull(glidNumber);
	  size_t positionWantedLocus = std::distance(lociAndResolutions.begin(), locusAndResolution);
	  glids.at(positionWantedLocus) = number;
	}
      locusName ++;
    }
}
				
void GLSReport::resolve(std::vector<std::shared_ptr<Report>> & listOfReports,
		       const GlidFile & glid,
		       const double minimalFrequency,
		       const bool resolveUnknownGenotype){

  try{
    for(auto glidNumber = glids.cbegin();
	glidNumber != glids.cend();
	glidNumber ++){

      if(*glidNumber == 0){
	if(resolveUnknownGenotype){
	  size_t positionGlidNumber = glidNumber - glids.cbegin();
	  genotypesWithFrequenciesAtLoci.push_back(glid.getPossibleGenotypesForAllLoci().at(positionGlidNumber).getGenotypes());
	}
	else{
	  throw MissingGenotypeException();
	}
      }
      else{
	auto itGlid = glid.getList().find(*glidNumber);
	if(itGlid != glid.getList().cend())
	  {
	    std::shared_ptr<Locus> pLocus = itGlid->second;
	    std::vector<std::pair<strArr_t, double>> genotypesAtLocus;
	    pLocus->reduce(genotypesAtLocus);
	    genotypesWithFrequenciesAtLoci.push_back(genotypesAtLocus);
	    types.push_back(pLocus->getType());
	  }
	else
	  {
	    throw MissingGlidException(*glidNumber);
	  }
      }//else glidNumber=0
      
      numberOfReports *= static_cast<double>(genotypesWithFrequenciesAtLoci.rbegin()->size());
      if(1./numberOfReports - minimalFrequency < ZERO){
	throw SplittingGenotypeException();
      }
    }//for glids

    buildListOfReports(listOfReports);
  }
  catch(const FileException & e)
    {
      throw;
    }
  catch(const std::exception & e){
      std::cout << e.what() << std::endl;      
      std::cout << "Id "
		<< id
		<< " discarded."
		<< std::endl;      
  }
}

void GLSCReport::translateLine(const std::string line){
	
  std::stringstream ss(line);
  std::string entry;
  if(ss >> entry)
    id = entry;

  while(ss >> entry){  
    singleLocusGenotypes.push_back(entry);
  }
}

void GLSCReport::resolve(std::vector<std::shared_ptr<Report>> & listOfReports){

  try
    { 
      doLociMatch();

      for(auto singleLocusGenotype : singleLocusGenotypes)
	{
	  std::string locusName = split(singleLocusGenotype, '*')[0];
	  auto locusAndResolution = lociAndResolutions.find(locusName);
	  
	  if(locusAndResolution != lociAndResolutions.cend())
	    {
	      size_t positionWantedLocus = std::distance(lociAndResolutions.begin(), locusAndResolution);
	      std::unique_ptr<Genotype> genotype = std::make_unique<GLSGenotype>(singleLocusGenotype, locusAndResolution->second);
	      
	      resolveSingleLocusGenotype(genotype,
					 positionWantedLocus);
	    }
	}
      buildListOfReports(listOfReports);
    }
  catch(const FileException & e)
    {
      throw;
    }
  catch(const std::exception & e)
    {
      std::cout << e.what() << std::endl;      
      std::cout << "Id "
		<< id
		<< " discarded."
		<< std::endl;      
    }	    
}

void GLSCReport::doLociMatch() const{

  for(auto locusAndResolution : lociAndResolutions)
    {
      auto pos = find_if(singleLocusGenotypes.cbegin(),
			 singleLocusGenotypes.cend(),
			 [& locusAndResolution](const std::string singleLocusGenotype)
			 {
			   std::string locusName = split(singleLocusGenotype, '*')[0];
			   return locusName == locusAndResolution.first;
			 }); 

      if(pos == singleLocusGenotypes.cend())
	{
	  throw NotMatchingLociException_GLSC(locusAndResolution.first, id);
	}
    }
}

void MACReport::translateLine(const std::string line){

  std::stringstream ss(line);
  if(ss >> id){}
  
  for(auto locusName = lociNamesFromFile.cbegin();
      locusName != lociNamesFromFile.cend();
      locusName ++)
    {
      std::string entry;
      std::string entry2;
      if(ss >> entry >> entry2)
	{
	  strArr_t locus;
	  locus.at(0) =  *locusName + '*' + entry;
	  locusName ++;
	  locus.at(1) =  *locusName + '*' + entry2;

	  lociFromFile.push_back(locus);
	}
    }
}

void MACReport::resolve(std::vector<std::shared_ptr<Report>> & listOfReports){

  try
    {
      if(lociFromFile.size() >= numberLoci)
	{
	  auto locusNameFromFile = lociNamesFromFile.cbegin();
	  for(auto singleLocusGenotype : lociFromFile){
	    
	    auto locusAndResolution = lociAndResolutions.find(*locusNameFromFile);

	    if(locusAndResolution != lociAndResolutions.cend())
	      { 
		size_t positionWantedLocus = std::distance(lociAndResolutions.begin(), locusAndResolution);
		std::unique_ptr<Genotype> genotype = std::make_unique<MACGenotype>(singleLocusGenotype, locusAndResolution->second);
		
		resolveSingleLocusGenotype(genotype,
					   positionWantedLocus);
	      }
	    locusNameFromFile ++;
	    locusNameFromFile ++;
	  }
	  
	  buildListOfReports(listOfReports);
	}
      else
	{
	  throw InputLineException();
	}
    }

  catch(const FileException & e)
    {
      throw;
    }
  catch(const std::exception & e)
    {
      std::cout << e.what() << std::endl;      
      std::cout << "Id "
		<< id
		<< " discarded."
		<< std::endl;      
    }	    
}
