#ifndef Report_header
#define Report_header

#include <array>
#include <vector>

#include "File.h"
#include "Typedefs.h"
#include "Allele.h"
#include "Phenotype.h"

class GlidFile;

class HaplotypeList;

class Report{

 public:
 Report(const Allele::codePrecision in_wantedPrecision, const size_t in_numberLoci)
   : genotypeAtLoci(),
    id(),
    frequency(),
    numberLoci(in_numberLoci),
    wantedPrecision(in_wantedPrecision)
    {
      findCombinations(numberLoci);
      genotypeAtLoci.reserve(numberLoci);
    }

  Report(const strArrVec_t & in_genotypeAtLoci,
	 const double in_frequency, 
	 const size_t in_numberLoci,
	 const std::string in_id)
    : genotypeAtLoci(in_genotypeAtLoci),
    id(in_id),
    frequency(in_frequency),
    numberLoci(in_numberLoci),
    wantedPrecision()
      {
	genotypeAtLoci.reserve(numberLoci);
      }

  void findCombinations(const size_t size);
  void writeCombinations() const;
  std::string buildPhenotypeCode() const;
  void buildHaploAndDiplotypes(PhenotypeList::iterator itPhenotype, HaplotypeList & haplotypeList) const;

  std::string getId() const {return id;}
  double getFrequency() const {return frequency;}
  const strArrVec_t & getGenotypeAtLoci() const {return genotypeAtLoci;}

 protected:
  strArrVec_t genotypeAtLoci;
  std::string id;
  double frequency;
  size_t numberLoci;
  Allele::codePrecision wantedPrecision;
  static std::vector<std::vector<bool>> haplotypeCombinations;
};

class GLReport : public Report{

 public:
 explicit GLReport(const std::string line,
		   const std::vector<bool> & booleanLociToDo,
		   const size_t numberLoci,
		   const Allele::codePrecision in_wantedPrecision) 
   : Report(in_wantedPrecision, numberLoci),
    inLoci()
      {
	translateLine(line, booleanLociToDo);
      }
  
  void translateLine(const std::string line, const std::vector<bool> & booleanLociToDo);
  void resolve(std::vector<GLReport> & listOfReports, const GlidFile & glid);
  
 private:
  std::vector<size_t> inLoci;
};

class HReport : public Report{
  
 public:
  explicit HReport(const std::string line,
		   const strVec_t & lociNames,
		   const size_t numberLoci,
		   const Allele::codePrecision in_wantedPrecision)
    : Report(in_wantedPrecision, numberLoci),
    inLoci()
      {
	translateLine(line, lociNames);
      }

  HReport(const strArrVec_t & in_genotypeAtLoci,
	  const double in_frequency,
	  const size_t in_numberLoci, 
	  const std::string in_id)
    : Report(in_genotypeAtLoci, in_frequency, in_numberLoci, in_id){}
  
  void translateLine(const std::string line, const strVec_t lociNames);
  void resolve(std::vector<HReport> & listOfReports);
  void resolveNMDPCode(const std::string code, strVec_t & newCodes) const;

 private:
  strArrVec_t inLoci;
  static FileNMDPCodes fileNMDPCodes;
};

#endif
