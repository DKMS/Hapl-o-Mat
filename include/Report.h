#ifndef Report_header
#define Report_header

#include <array>
#include <vector>
#include <fstream>

#include "File.h"
#include "Typedefs.h"
#include "Allele.h"
#include "Phenotype.h"
#include "Locus.h"

class GlidFile;
class HaplotypeList;
class HaplotypeCombinations;

class BasicReport{

 public:
  explicit BasicReport(const size_t in_numberLoci)
    : genotypeAtLoci(),
    id(),
    frequency(),
    numberLoci(in_numberLoci){}

 protected:
  strArrVec_t genotypeAtLoci;
  std::string id;
  double frequency;
  size_t numberLoci;

};

class Report{

 public:
 explicit Report(const Allele::codePrecision in_wantedPrecision,
		 const size_t in_numberLoci)
   : genotypeAtLoci(),
    id(),
    frequency(),
    numberLoci(in_numberLoci),
    wantedPrecision(in_wantedPrecision),
    types()
    {
      genotypeAtLoci.reserve(numberLoci);
    }
  explicit Report(const strArrVec_t & in_genotypeAtLoci,
		  const double in_frequency, 
		  const size_t in_numberLoci,
		  const std::string in_id,
		  const std::vector<Locus::reportType> & in_types)
    : genotypeAtLoci(in_genotypeAtLoci),
    id(in_id),
    frequency(in_frequency),
    numberLoci(in_numberLoci),
    wantedPrecision(),
    types(in_types)
      {
	genotypeAtLoci.reserve(numberLoci);
      }
  virtual std::shared_ptr<Report> create(const strArrVec_t & in_genotypeAtLoci,
					 const double in_frequency, 
					 const size_t in_numberLoci,
					 const std::string in_id,
					 const std::vector<Locus::reportType> & in_types) = 0;
  virtual ~Report(){}

  std::string buildPhenotypeCode() const;
  void buildHaploAndDiplotypes(PhenotypeList::iterator itPhenotype,
			       HaplotypeList & haplotypeList,
			       std::ofstream & haplotypesFile,
			       const HaplotypeCombinations & haplotypeCombinations) const; 
  void buildListOfReports(std::vector<std::shared_ptr<Report>> & listOfReports,
			  const std::vector<std::vector<std::pair<strArr_t, double>>> & genotypesAtLoci);
  std::string evaluateReportType(const size_t numberReports) const;
  
  std::string getId() const {return id;}
  double getFrequency() const {return frequency;}
  const strArrVec_t & getGenotypeAtLoci() const {return genotypeAtLoci;}
  static double getNumberH0Reports() {return numberH0Reports;}
  static double getNumberH1Reports() {return numberH1Reports;}
  static double getNumberH2Reports() {return numberH2Reports;}
  static double getNumberH2MReports() {return numberH2MReports;}
  static double getNumberIReports() {return numberIReports;}

 protected:
  strArrVec_t genotypeAtLoci;
  std::string id;
  double frequency;
  size_t numberLoci;
  Allele::codePrecision wantedPrecision;
  std::vector<Locus::reportType> types;
  static double numberH0Reports;
  static double numberH1Reports;
  static double numberH2Reports;
  static double numberH2MReports;
  static double numberIReports;
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
  explicit GLReport(const strArrVec_t & in_genotypeAtLoci,
		    const double in_frequency,
		    const size_t in_numberLoci, 
		    const std::string in_id,
		    const std::vector<Locus::reportType> & in_types)
    : Report(in_genotypeAtLoci, in_frequency, in_numberLoci, in_id, in_types){}
  virtual std::shared_ptr<Report> create(const strArrVec_t & in_genotypeAtLoci,
					 const double in_frequency, 
					 const size_t in_numberLoci,
					 const std::string in_id,
					 const std::vector<Locus::reportType> & in_types)
    {
      std::shared_ptr<Report> pReport = std::make_shared<GLReport> (in_genotypeAtLoci,
								    in_frequency, 
								    in_numberLoci,
								    in_id,
								    in_types);
      return pReport;
    }
  
  void translateLine(const std::string line, const std::vector<bool> & booleanLociToDo);
  void resolve(std::vector<std::shared_ptr<Report>> & listOfReports,
	       const GlidFile & glid,
	       const double minimalFrequency, 
	       const bool resolveUnknownGenotype);
  
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
  explicit HReport(const strArrVec_t & in_genotypeAtLoci,
		   const double in_frequency,
		   const size_t in_numberLoci, 
		   const std::string in_id,
		   const std::vector<Locus::reportType> & in_types)
    : Report(in_genotypeAtLoci, in_frequency, in_numberLoci, in_id, in_types){}
  virtual std::shared_ptr<Report> create(const strArrVec_t & in_genotypeAtLoci,
					 const double in_frequency, 
					 const size_t in_numberLoci,
					 const std::string in_id,
					 const std::vector<Locus::reportType> & in_types)
    {
      std::shared_ptr<Report> pReport = std::make_shared<HReport> (in_genotypeAtLoci,
								   in_frequency, 
								   in_numberLoci,
								   in_id,
								   in_types);
      return pReport;
    }
  
  void translateLine(const std::string line, const strVec_t lociNames);
  void resolve(std::vector<std::shared_ptr<Report>> & listOfReports,
	       const double minimalFrequency,
	       const bool doH2Filter,
	       const bool expandh2Lines);
  void resolveNMDPCode(const std::string code, strVec_t & newCodes) const;

 private:
  strArrVec_t inLoci;
  static std::unordered_map<std::string, std::shared_ptr<Locus>> lociAlreadyDone;
  static FileNMDPCodes fileNMDPCodes;
};

#endif
