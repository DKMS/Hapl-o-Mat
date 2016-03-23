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

#ifndef Report_header
#define Report_header

#include <array>
#include <vector>
#include <fstream>

#include "Typedefs.h"
#include "Allele.h"
#include "Phenotype.h"
#include "Locus.h"
#include "Genotypes.h"

class GlidFile;
class Haplotypes;
class HaplotypeCombinations;

class BasicReport{

 public:
  explicit BasicReport(const size_t in_numberLoci)
    :id(),
    frequency(),
    numberLoci(in_numberLoci),
    genotypeAtLoci(){}

  std::string buildPhenotypeCode() const;
  void buildHaploAndDiplotypes(Phenotypes::iterator itPhenotype,
			       Haplotypes & haplotypes,
			       std::ofstream & haplotypesFile,
			       const HaplotypeCombinations & haplotypeCombinations) const; 

  std::string getId() const {return id;}
  double getFrequency() const {return frequency;}
  const strArrVec_t & getGenotypeAtLoci() const {return genotypeAtLoci;}

 protected:
  std::string id;
  double frequency;
  size_t numberLoci;
  strArrVec_t genotypeAtLoci;
};

class ReadinReport : public BasicReport{

 public:
  ReadinReport(const std::string line,
	       const size_t in_numberLoci)
    : BasicReport(in_numberLoci)
  {
    translateLine(line);
  }

 private:
  void translateLine(const std::string line);
};

class Report : public BasicReport{

 public:
 explicit Report(const std::map<std::string, Allele::codePrecision> & in_lociAndResolutions)
   : BasicReport(in_lociAndResolutions.size()),
    lociAndResolutions(in_lociAndResolutions),
    types(),
    numberOfReports(1.),
    discardReport(false),
    genotypesWithFrequenciesAtLoci()
      {
	genotypeAtLoci.reserve(numberLoci);
      }
  explicit Report(const strArrVec_t & in_genotypeAtLoci,
		  const double in_frequency, 
		  const size_t in_numberLoci,
		  const std::string in_id,
		  const std::vector<Locus::reportType> & in_types)
    : BasicReport(in_numberLoci),
    lociAndResolutions(),
    types(in_types)
      {
	genotypeAtLoci = in_genotypeAtLoci;
	id = in_id;
	frequency = in_frequency;
      }
  virtual std::shared_ptr<Report> create(const strArrVec_t & in_genotypeAtLoci,
					 const double in_frequency, 
					 const size_t in_numberLoci,
					 const std::string in_id,
					 const std::vector<Locus::reportType> & in_types) = 0;
  virtual ~Report(){}

  void buildListOfReports(std::vector<std::shared_ptr<Report>> & listOfReports);
  std::string evaluateReportType(const size_t numberReports) const;

  static double getNumberNReports() {return numberNReports;}
  static double getNumberAReports() {return numberAReports;}
  static double getNumberMReports() {return numberMReports;}
  static double getNumberIReports() {return numberIReports;}

 protected:
  std::map<std::string, Allele::codePrecision> lociAndResolutions;
  std::vector<Locus::reportType> types;
  double numberOfReports;
  bool discardReport;
  std::vector<std::vector<std::pair<strArr_t, double>>> genotypesWithFrequenciesAtLoci;
  static double numberNReports;
  static double numberAReports;
  static double numberMReports;
  static double numberIReports;
};

class GLReport : public Report{

 public:
  explicit GLReport(const std::string line,
		    const strVec_t & in_lociOrder,
		    const std::map<std::string, Allele::codePrecision> & in_lociAndResolutions) 
    : Report(in_lociAndResolutions),
    lociOrder(in_lociOrder)
    {
      translateLine(line);
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
  
  void translateLine(const std::string line);
  void resolve(std::vector<std::shared_ptr<Report>> & listOfReports,
	       const GlidFile & glid,
	       const double minimalFrequency, 
	       const bool resolveUnknownGenotype);
  
 private:
  strVec_t lociOrder;
  std::vector<size_t> glids;
};

class ColumnReport : public Report{

 public:
  explicit ColumnReport(const std::map<std::string,
			Allele::codePrecision> & in_lociAndResolutions,
			const double in_minimalFrequency,
			const bool in_doAmbiguityFilter,
			const bool in_expandAmbiguityLines)
    : Report(in_lociAndResolutions),
    minimalFrequency(in_minimalFrequency),
    doAmbiguityFilter(in_doAmbiguityFilter),
    expandAmbiguityLines(in_expandAmbiguityLines)
    {
      genotypesWithFrequenciesAtLoci.resize(numberLoci);
    }

  explicit ColumnReport(const strArrVec_t & in_genotypeAtLoci,
			const double in_frequency,
			const size_t in_numberLoci, 
			const std::string in_id,
			const std::vector<Locus::reportType> & in_types)
    : Report(in_genotypeAtLoci, in_frequency, in_numberLoci, in_id, in_types){}
  
  virtual void translateLine(const std::string line) = 0;
  virtual void resolve(std::vector<std::shared_ptr<Report>> & listOfReports) = 0;

  void resolveSingleLocusGenotype(const std::unique_ptr<Genotype> & genotype,
				  const size_t positionWantedLocus);
  
 protected:
  double minimalFrequency;
  bool doAmbiguityFilter;
  bool expandAmbiguityLines;
  static std::unordered_map<std::string, std::shared_ptr<Locus>> singleLocusGenotypesAlreadyDone;
};

class GLCReport : public ColumnReport{

 public:
  explicit GLCReport(const std::string line,
		     const std::map<std::string, Allele::codePrecision> & in_lociAndResolutions,
		     const double in_minimalFrequency,
		     const bool in_doAmbiguityFilter,
		     const bool in_expandAmbiguityLines)
    : ColumnReport(in_lociAndResolutions, in_minimalFrequency, in_doAmbiguityFilter, in_expandAmbiguityLines),
    singleLocusGenotypes()
  {
    translateLine(line);
    types.resize(numberLoci);
  }
  explicit GLCReport(const strArrVec_t & in_genotypeAtLoci,
		     const double in_frequency,
		     const size_t in_numberLoci, 
		     const std::string in_id,
		     const std::vector<Locus::reportType> & in_types)
    : ColumnReport(in_genotypeAtLoci, in_frequency, in_numberLoci, in_id, in_types){}
  
  virtual std::shared_ptr<Report> create(const strArrVec_t & in_genotypeAtLoci,
					 const double in_frequency, 
					 const size_t in_numberLoci,
					 const std::string in_id,
					 const std::vector<Locus::reportType> & in_types)
    {
      std::shared_ptr<Report> pReport = std::make_shared<GLCReport> (in_genotypeAtLoci,
								    in_frequency, 
								    in_numberLoci,
								    in_id,
								    in_types);
      return pReport;
    }


  virtual void translateLine(const std::string line);
  virtual void resolve(std::vector<std::shared_ptr<Report>> & listOfReports);

 private:
  strVec_t singleLocusGenotypes;

};

class MAReport : public ColumnReport{
  
 public:
  explicit MAReport(const std::string line,
		   const strVec_t & in_lociNamesFromFile,
		   const std::map<std::string, Allele::codePrecision> & in_lociAndResolutions,
		   const double in_minimalFrequency,
		   const bool in_doAmbiguityFilter,
		   const bool in_expandAmbiguityLines)
    : ColumnReport(in_lociAndResolutions, in_minimalFrequency, in_doAmbiguityFilter, in_expandAmbiguityLines),
    lociFromFile(),
    lociNamesFromFile(in_lociNamesFromFile)
      {
	translateLine(line);
	types.resize(numberLoci);
      }
  explicit MAReport(const strArrVec_t & in_genotypeAtLoci,
		   const double in_frequency,
		   const size_t in_numberLoci, 
		   const std::string in_id,
		   const std::vector<Locus::reportType> & in_types)
    : ColumnReport(in_genotypeAtLoci, in_frequency, in_numberLoci, in_id, in_types){}

  virtual std::shared_ptr<Report> create(const strArrVec_t & in_genotypeAtLoci,
					 const double in_frequency, 
					 const size_t in_numberLoci,
					 const std::string in_id,
					 const std::vector<Locus::reportType> & in_types)
    {
      std::shared_ptr<Report> pReport = std::make_shared<MAReport> (in_genotypeAtLoci,
								   in_frequency, 
								   in_numberLoci,
								   in_id,
								   in_types);
      return pReport;
    }
  
  virtual void translateLine(const std::string line);
  virtual void resolve(std::vector<std::shared_ptr<Report>> & listOfReports);

 private:
  strArrVec_t lociFromFile;
  strVec_t lociNamesFromFile;
};



#endif
