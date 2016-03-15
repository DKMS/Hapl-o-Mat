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

#include <sstream>
#include <fstream>
#include <iostream>

#include "Report.h"
#include "Locus.h"
#include "Glid.h"
#include "Utility.h"
#include "Phenotype.h"
#include "Haplotype.h"
#include "DataProcessing.h"

std::unordered_map<std::string, std::shared_ptr<Locus>> ColumnReport::singleLocusGenotypesAlreadyDone;
double Report::numberH0Reports = 0.;
double Report::numberH1Reports = 0.;
double Report::numberH2Reports = 0.;
double Report::numberH2MReports = 0.;
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
					      const size_t positionWantedLocus,
					      std::vector<std::vector<std::pair<strArr_t, double>>> & genotypesAtLoci){
					  
  std::vector<std::pair<strArr_t, double>> genotypesAtLocus;

  auto pos = singleLocusGenotypesAlreadyDone.find(genotype->getSingleLocusGenotype());
  if(pos == singleLocusGenotypesAlreadyDone.cend())
    {
      std::shared_ptr<Locus> pLocus = genotype->resolve(doH2Filter, expandH2Lines);
    
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
  if(1./numberOfReports - minimalFrequency < ZERO){
    std::cout << "Report "
	      << id
	      << " comes below allowed frequency. Report discarded."
	      << std::endl;
    discardReport = true;
  }
  else{
    genotypesAtLoci.at(positionWantedLocus) = genotypesAtLocus;
  }
}

void Report::buildListOfReports(std::vector<std::shared_ptr<Report>> & listOfReports,
				const std::vector<std::vector<std::pair<strArr_t,double>>> & genotypesAtLoci){

  std::vector<std::vector<std::pair<strArr_t, double>>> reports;
  cartesianProduct(reports, genotypesAtLoci);
  
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
	       Locus::reportType::H2M) != types.cend())
    numberH2MReports += 1./static_cast<double>(numberReports);
  else if(find(types.cbegin(),
	       types.cend(),
	       Locus::reportType::H2) != types.cend())
    numberH2Reports += 1./static_cast<double>(numberReports);
  else if(find(types.cbegin(),
	       types.cend(),
	       Locus::reportType::H1) != types.cend())
    numberH1Reports += 1./static_cast<double>(numberReports);
  else if(find(types.cbegin(),
	       types.cend(),
	       Locus::reportType::H0) != types.cend())
    numberH0Reports += 1./static_cast<double>(numberReports);

  std::string totalType = "";
  for(auto type : types){
    switch(type){
    case Locus::reportType::H0:
      {
	totalType += "H0";
	break;
      }
    case Locus::reportType::H1:
      {
	totalType += "H1";
	break;
      }
    case Locus::reportType::H2:
      {	totalType += "H2";
	break;
      }
    case Locus::reportType::H2M:
      {
	totalType += "H2M";
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

void GLReport::translateLine(const std::string line){

  id = leftOfFirstDelim(line, ';');

  std::string rightPartOfLine = rightOfFirstDelim(line, ';');
  strVec_t allGlids = split(rightPartOfLine, ':');

  glids.resize(numberLoci);

  auto locusName = lociOrder.cbegin();
  for(auto glidNumber : allGlids)
    {
      auto locusAndWantedAlleleGroup = lociAndWantedAlleleGroups.find(*locusName);
      if(locusAndWantedAlleleGroup != lociAndWantedAlleleGroups.cend())
	{
	  size_t number = stoull(glidNumber);
	  size_t positionWantedLocus = std::distance(lociAndWantedAlleleGroups.begin(), locusAndWantedAlleleGroup);
	  glids.at(positionWantedLocus) = number;
	}
      locusName ++;
    }
}
				
void GLReport::resolve(std::vector<std::shared_ptr<Report>> & listOfReports,
		       const GlidFile & glid,
		       const double minimalFrequency,
		       const bool resolveUnknownGenotype){

  std::vector<std::vector<std::pair<strArr_t, double>>> genotypesAtLoci;

  for(auto glidNumber = glids.cbegin();
      glidNumber != glids.cend();
      glidNumber ++){

    if(*glidNumber == 0){
      if(resolveUnknownGenotype){
	size_t positionGlidNumber = glidNumber - glids.cbegin();
	genotypesAtLoci.push_back(glid.getPossibleGenotypesForAllLoci().at(positionGlidNumber).getGenotypes());
      }
      else{
	discardReport = true;
	std::cout << "Report "
		  << id
		  << " contains Glid 0. Report discarded."
		  << std::endl;
	break;
      }
    }
    else{
      auto itGlid = glid.getList().find(*glidNumber);
      if(itGlid == glid.getList().cend()){
	std::cerr << "Key "
		  << *glidNumber
		  << " not in glid-file" << std::endl;
	exit(EXIT_FAILURE);
      }
      else{
	std::shared_ptr<Locus> pLocus = itGlid->second;
	std::vector<std::pair<strArr_t, double>> genotypesAtLocus;
	pLocus->reduce(genotypesAtLocus);
	genotypesAtLoci.push_back(genotypesAtLocus);
	types.push_back(pLocus->getType());
      }
    }//else glidNumber=0

    numberOfReports *= static_cast<double>(genotypesAtLoci.rbegin()->size());
    if(1./numberOfReports - minimalFrequency < ZERO){
      discardReport = true;
      std::cout << "Report "
		<< id
		<< " comes below allowed frequency. Report discarded."
		<< std::endl;
      break;
    }
    if(numberOfReports < ZERO){
      discardReport = true;	    
      std::cout << "Report "
		<< id
		<< " contains empty loci."
		<< std::endl;
      break;
    }
  }//for glids

  if(!discardReport)
    buildListOfReports(listOfReports, genotypesAtLoci);
}

void GLCReport::translateLine(const std::string line){
	
  std::stringstream ss(line);
  std::string entry;
  if(ss >> entry)
    id = entry;

  while(ss >> entry){  
    singleLocusGenotypes.push_back(entry);
  }
}

void GLCReport::resolve(std::vector<std::shared_ptr<Report>> & listOfReports){

  std::vector<std::vector<std::pair<strArr_t, double>>> genotypesAtLoci;
  genotypesAtLoci.resize(numberLoci);

  for(auto singleLocusGenotype : singleLocusGenotypes){

    if(!discardReport)
      {
	std::string locusName = split(singleLocusGenotype, '*')[0];
	auto locusAndWantedAlleleGroup = lociAndWantedAlleleGroups.find(locusName);
    
	if(locusAndWantedAlleleGroup != lociAndWantedAlleleGroups.cend())
	  {
	    size_t positionWantedLocus = std::distance(lociAndWantedAlleleGroups.begin(), locusAndWantedAlleleGroup);
	    std::unique_ptr<Genotype> genotype = make_unique<GLGenotype>(singleLocusGenotype, locusAndWantedAlleleGroup->second);
	    
	    resolveSingleLocusGenotype(genotype,
				       positionWantedLocus,
				       genotypesAtLoci);
	  }
      }
  }

  if(!discardReport){  
    buildListOfReports(listOfReports, genotypesAtLoci);
  }
}


void MAReport::translateLine(const std::string line){

  std::stringstream ss(line);
  std::string entry;
  if(ss >> entry)
    id = entry;

  auto locusName = lociNamesFromFile.cbegin();
  std::string entry2;
  while(ss >> entry >> entry2){
    std::string allele1 = *locusName + '*';
    locusName ++;
    std::string allele2 = *locusName + '*';
    locusName ++;
    allele1.append(entry);
    allele2.append(entry2);
    strArr_t locus;
    locus.at(0) = allele1;
    locus.at(1) = allele2;
    lociFromFile.push_back(locus);
  }
}

void MAReport::resolve(std::vector<std::shared_ptr<Report>> & listOfReports){

  std::vector<std::vector<std::pair<strArr_t, double>>> genotypesAtLoci;
  genotypesAtLoci.resize(numberLoci);

  auto locusNameFromFile = lociNamesFromFile.cbegin();
  for(auto singleLocusGenotype = lociFromFile.begin();
      singleLocusGenotype != lociFromFile.end();
      singleLocusGenotype ++){

    if(!discardReport)
      {
	auto locusAndWantedAlleleGroup = lociAndWantedAlleleGroups.find(*locusNameFromFile);

	if(locusAndWantedAlleleGroup != lociAndWantedAlleleGroups.cend())
	  { 
	    size_t positionWantedLocus = std::distance(lociAndWantedAlleleGroups.begin(), locusAndWantedAlleleGroup);
	    std::unique_ptr<Genotype> genotype = make_unique<MAGenotype>(*singleLocusGenotype, locusAndWantedAlleleGroup->second);
	    
	    resolveSingleLocusGenotype(genotype,
				       positionWantedLocus,
				       genotypesAtLoci);
	  }
	locusNameFromFile ++;
	locusNameFromFile ++;
      }
  }
  
  if(!discardReport)
    buildListOfReports(listOfReports, genotypesAtLoci);
}
