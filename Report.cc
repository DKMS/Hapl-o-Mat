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

FileNMDPCodes HReport::fileNMDPCodes("data/code2dna.txt");
std::unordered_map<std::string, std::shared_ptr<Locus>> HReport::lociAlreadyDone;
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

void BasicReport::buildHaploAndDiplotypes(PhenotypeList::iterator itPhenotype,
				     HaplotypeList & haplotypeList,
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
      std::pair<HaplotypeList::iterator, bool> inserted1 = haplotypeList.add(codeHaplotype1);
      std::pair<HaplotypeList::iterator, bool> inserted2 = haplotypeList.add(codeHaplotype2);
      
      if(inserted1.second){
	haplotypesFile << codeHaplotype1 << std::endl;
      }
      else{
	inserted1.first->second.incrementNumber();
      }
      if(inserted2.second){
	haplotypesFile << codeHaplotype2 << std::endl;
      }
      else{
	inserted2.first->second.incrementNumber();
      }
      
      //build diplotype                                                                                                                          
      size_t id1 = inserted1.first->first;
      size_t id2 = inserted2.first->first;
      
      Diplotype diplotype;
      diplotype.id1 = id1;
      diplotype.id2 = id2;
      if(diplotype.id1 == diplotype.id2)
	diplotype.sameHaplotype = true;
      else
	diplotype.sameHaplotype = false;
      
      itPhenotype->second.addDiplotype(diplotype);
    }//haplotypeCombinations  
}

void ReadinReport::translateLine(const std::string line){

  std::stringstream ss(line);
  std::string type;
  std::string code;
  if(ss >> id >> type >> frequency >> code){}

  strVec_t genotypes = split(code, '^');
  for(auto genotype : genotypes){
    strVec_t alleles = split(genotype, '+');
    sort(alleles.begin(), alleles.end());
    strArr_t genotypeAtLocus;
    genotypeAtLocus.at(0) = alleles.at(0);
    genotypeAtLocus.at(1) = alleles.at(1);
    genotypeAtLoci.push_back(genotypeAtLocus);
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
      {
	totalType += "H2";
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

void GLReport::translateLine(const std::string line, const std::vector<bool> & booleanLociToDo){

  id = leftOfFirstDelim(line, ';');

  std::string rightPartOfLine = rightOfFirstDelim(line, ';');
  strVec_t codes = split(rightPartOfLine, ':');
  inLoci.reserve(codes.size());
  size_t counter = 0;
  for(auto code : codes){
    if(booleanLociToDo.at(counter)){
      size_t number = stoull(code);
      inLoci.push_back(number);
    }
    counter ++;
  }
}
				
void GLReport::resolve(std::vector<std::shared_ptr<Report>> & listOfReports,
		       const GlidFile & glid,
		       const double minimalFrequency,
		       const bool resolveUnknownGenotype){

  std::vector<std::vector<std::pair<strArr_t, double>>> genotypesAtLoci;
  bool discardReport = false;

  for(auto code : inLoci){
    if(code == 0){
      if(resolveUnknownGenotype){
	genotypesAtLoci.push_back(glid.getPossibleGenotypesForAllLoci().find(genotypesAtLoci.size())->second.getGenotypes());
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
      auto itGlid = glid.getList().find(code);
      if(itGlid == glid.getList().cend()){
	std::cerr << "Key "
		  << code
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
    }//else code=0
  }//for inLoci

  for(auto locus : genotypesAtLoci){
    if(locus.empty()){
      discardReport = true;	    
      std::cout << "Report "
		<< id
		<< " contains empty loci."
		<< std::endl;
      break;
    }
  }
  
  if(!discardReport){
    double numberOfReports = 1.;
    for(auto locus : genotypesAtLoci)
      numberOfReports *= static_cast<double>(locus.size());
    
    if(1./numberOfReports - minimalFrequency < ZERO){
      std::cout << "Report "
		<< id
		<< " with "
		<< numberOfReports
		<< " phenotypes comes below allowed frequency. Report discarded."
		<< std::endl;
    }
    else
      buildListOfReports(listOfReports, genotypesAtLoci);
  }
}

void HReport::translateLine(const std::string line, const strVec_t lociNames){

  std::stringstream ss(line);
  std::string entry;
  if(ss >> entry)
    id = entry;

  auto locusName = lociNames.cbegin();
  std::string entry2;
  while(ss >> entry >> entry2){
    std::string code1 = *locusName;
    locusName ++;
    std::string code2 = *locusName;
    locusName ++;
    code1.append(entry);
    code2.append(entry2);
    strArr_t locus;
    locus.at(0) = code1;
    locus.at(1) = code2;
    inLoci.push_back(locus);
  }
}

void HReport::resolveNMDPCode(const std::string code, strVec_t & newCodes) const{

  std::string nmdpCode = findNMDPCode(code);
  auto itFileNMDPCodes = fileNMDPCodes.getList().find(nmdpCode);
  if(itFileNMDPCodes == fileNMDPCodes.getList().cend()){
    std::cerr << "Could not find NMDP-Code "
	      << nmdpCode
	      << std::endl;
    exit (EXIT_FAILURE);
  }
  else{
    std::string newCode = code;
    size_t positionAsterik = code.find('*') + 1;
    size_t positionNMDPCodeInCode = code.find(nmdpCode, positionAsterik);
    newCode.erase(positionNMDPCodeInCode);
    if(itFileNMDPCodes->second.find(':') != std::string::npos){
      std::size_t posLastColon = newCode.find_last_of(':');
      newCode.erase(posLastColon);
      posLastColon = newCode.find_last_of(':');
      if(posLastColon == std::string::npos)
	posLastColon = newCode.find_last_of('*');
      newCode.erase(posLastColon+1);
    }

    strVec_t splittedCode = split(itFileNMDPCodes->second, '/');
    for(auto itSplittedCode : splittedCode)
      {
	std::string newCode2 = newCode;
	newCode2.append(itSplittedCode);
	newCodes.push_back(newCode2);
      }//for splittedCode
  }//else

  if(newCodes.empty()){
    std::cerr << "Did not find allele from multi allele code "
	     << nmdpCode
	     << " in allAlleles.txt. Please update allAlleles.txt."
	     <<std::endl;
    exit(EXIT_FAILURE);
  }
}

void HReport::resolve(std::vector<std::shared_ptr<Report>> & listOfReports,
		      const double minimalFrequency,
		      const bool doH2Filter,
		      const bool expandH2Lines){

  std::vector<std::vector<std::pair<strArr_t, double>>> genotypesAtLoci;

  for(auto locus : inLoci){
    std::string locusCombination = "";
    for(auto code : locus){
      locusCombination += code;
      locusCombination += "+";    
    }
    locusCombination.pop_back();

    auto pos = lociAlreadyDone.find(locusCombination);
    if(pos == lociAlreadyDone.cend()){
      strVecArr_t locusPositions;
      size_t counter = 0;
      for(auto code : locus){
	strVec_t codes;
	if(checkNMDPCode(code)){
	  resolveNMDPCode(code, codes);
	}
	else{
	  codes.push_back(code);
	}
	locusPositions.at(counter) = codes;
	counter ++;
      }
      
      std::shared_ptr<Locus> pLocus (new UnphasedLocus(locusPositions, wantedPrecision, doH2Filter, expandH2Lines));
      pLocus->resolve();
      lociAlreadyDone.emplace(locusCombination, pLocus);

      types.push_back(pLocus->getType());
      std::vector<std::pair<strArr_t, double>> genotypesAtLocus;
      pLocus->reduce(genotypesAtLocus);
      genotypesAtLoci.push_back(genotypesAtLocus);
    }
    else{
      types.push_back(pos->second->getType());
      std::vector<std::pair<strArr_t, double>> genotypesAtLocus;
      pos->second->reduce(genotypesAtLocus);
      genotypesAtLoci.push_back(genotypesAtLocus);
    }
  }//for inLoci
  
  double numberOfReports = 1.;
  for(auto locus : genotypesAtLoci)
    numberOfReports *= static_cast<double>(locus.size());

  if(1./numberOfReports - minimalFrequency < ZERO){
    std::cout << "Report "
              << id
              << " with "
              << numberOfReports
              << " phenotypes comes below allowed frequency. Report discarded."
              << std::endl;
  }
  else
    buildListOfReports(listOfReports, genotypesAtLoci);
}

