#ifndef Allele_header
#define Allele_header

#include <string>
#include <memory>
#include <cmath>

#include "File.h"
#include "Typedefs.h"

class Allele{

 public:
  enum codePrecision{
    g,
    fourDigit,
    G,
    GForH2Filter,
    sixDigit,
    eightDigit,
    asItIs
  };

 explicit Allele(const std::string in_code,
		 const codePrecision in_precision,
		 const codePrecision in_wantedPrecision,
		 const double in_frequency)
   : code(in_code),
    precision(in_precision),
    wantedPrecision(in_wantedPrecision),
    frequency(in_frequency){}
  virtual std::shared_ptr<Allele> create(const std::string in_code,
					 const codePrecision in_precision,
					 const codePrecision in_wantedPrecision,
					 const double in_frequency) = 0;
  virtual ~Allele(){}

  virtual std::vector<std::shared_ptr<Allele>> translateTog() = 0; 
  virtual std::vector<std::shared_ptr<Allele>> translateTo4d() = 0;
  virtual std::vector<std::shared_ptr<Allele>> translateToG(const FileAllelesTogOrG & whichH1File) = 0;
  virtual std::vector<std::shared_ptr<Allele>> translateTo6d() = 0;
  //  virtual std::vector<std::shared_ptr<Allele>> translateTo8d() = 0;

  static std::shared_ptr<Allele> createAllele(const std::string code, const Allele::codePrecision wantedPrecision, const double alleleFrequency);
  static Allele::codePrecision identifyCodePrecision(const std::string code);
  std::vector<std::shared_ptr<Allele>> translate();
  std::vector<std::shared_ptr<Allele>>::iterator ispAlleleInList(std::vector<std::shared_ptr<Allele>> & listOfpAlleles) const;
  std::string allelesTog();
  strVec_t allelesToG(const FileAllelesTogOrG & whichH1File);
  strVec_t fourDigitOrgToG();
  strVec_t GToAlleles();
  strVec_t gToAlleles();
  strVec_t fourDigitToEightDigit();
  static std::string printCodePrecision(const codePrecision precision);
  double getFrequency() const {return frequency;}
  void multiplyFrequency(const double factor) {frequency *= factor;}
  void addFrequency(const double factor) {frequency += factor;}
  void sqrtFrequency() {frequency = sqrt(frequency);}
  std::string getCode() const {return code;}
  codePrecision getPrecision() const {return precision;}
  codePrecision getWantedPrecision() const {return wantedPrecision;}

 protected:
  std::string code;
  codePrecision precision;
  codePrecision wantedPrecision;
  double frequency;
  static FileAllelesTogOrG fileAllelesTog;
  static FileAllelesTogOrG fileAllelesToG;
  static FileAllelesTogOrG fileAllelesToGForH2Filter;
  static FilegToG filegToG;
  static FilegOrGOr4dToAlleles fileGToAlleles;
  static FilegOrGOr4dToAlleles filegToAlleles;
  static FilegOrGOr4dToAlleles file4dToAlleles;
};

class Alleleg : public Allele{

 public:
  explicit Alleleg(const std::string in_code,
		   const codePrecision in_precision,
		   const codePrecision in_wantedPrecision,
		   const double in_frequency)
    : Allele(in_code, in_precision, in_wantedPrecision, in_frequency){}
  explicit Alleleg(const std::string in_code,
		   const double in_frequency)
    : Allele(in_code, codePrecision::g, codePrecision::g, in_frequency){}
  virtual std::shared_ptr<Allele> create(const std::string in_code,
					 const codePrecision in_precision,
					 const codePrecision in_wantedPrecision,
					 const double in_frequency)
    {
      std::shared_ptr<Allele> pAllele = std::make_shared<Alleleg> (in_code,
								   in_precision,
								   in_wantedPrecision,
								   in_frequency);
      return pAllele;
    }

  virtual std::vector<std::shared_ptr<Allele>> translateTog();
  virtual std::vector<std::shared_ptr<Allele>> translateTo4d();
  virtual std::vector<std::shared_ptr<Allele>> translateToG(const FileAllelesTogOrG & whichH1File);
  virtual std::vector<std::shared_ptr<Allele>> translateTo6d(){
    std::vector<std::shared_ptr<Allele>> nix;
    return nix;};
  //  virtual std::vector<std::shared_ptr<Allele>> translateTo8d(){
  //  std::vector<std::shared_ptr<Allele>> nix;
  //  return nix;};  

 private:
};
 
class Allele4d : public Allele{

 public:
  explicit Allele4d(const std::string in_code,
		    const codePrecision in_precision,
		    const codePrecision in_wantedPrecision,
		    const double in_frequency)
    : Allele(in_code, in_precision, in_wantedPrecision, in_frequency){}
  explicit Allele4d(const std::string in_code,
		    const double in_frequency)
    : Allele(in_code, codePrecision::fourDigit, codePrecision::fourDigit, in_frequency){}
  virtual std::shared_ptr<Allele> create(const std::string in_code,
					 const codePrecision in_precision,
					 const codePrecision in_wantedPrecision,
					 const double in_frequency)
    {
      std::shared_ptr<Allele> pAllele = std::make_shared<Allele4d> (in_code,
								    in_precision,
								    in_wantedPrecision,
								    in_frequency);
      return pAllele;
    }

  virtual std::vector<std::shared_ptr<Allele>> translateTog();
  virtual std::vector<std::shared_ptr<Allele>> translateTo4d();
  virtual std::vector<std::shared_ptr<Allele>> translateToG(const FileAllelesTogOrG & whichH1File);
  virtual std::vector<std::shared_ptr<Allele>> translateTo6d();
  //  virtual std::vector<std::shared_ptr<Allele>> translateTo8d();  

 private:

};

class AlleleG : public Allele{

 public:
  explicit AlleleG(const std::string in_code,
		   const codePrecision in_precision,
		   const codePrecision in_wantedPrecision,
		   const double in_frequency)
    : Allele(in_code, in_precision, in_wantedPrecision, in_frequency){}
  explicit AlleleG(const std::string in_code,
		   const double in_frequency)
    : Allele(in_code, codePrecision::G, codePrecision::G, in_frequency){}
  virtual std::shared_ptr<Allele> create(const std::string in_code,
					 const codePrecision in_precision,
					 const codePrecision in_wantedPrecision,
					 const double in_frequency)
    {
      std::shared_ptr<Allele> pAllele = std::make_shared<AlleleG> (in_code,
								   in_precision,
								   in_wantedPrecision,
								   in_frequency);
      return pAllele;
    }
  
  virtual std::vector<std::shared_ptr<Allele>> translateTog();
  virtual std::vector<std::shared_ptr<Allele>> translateTo4d();
  virtual std::vector<std::shared_ptr<Allele>> translateToG(const FileAllelesTogOrG & whichH1File);
  virtual std::vector<std::shared_ptr<Allele>> translateTo6d(){
    std::vector<std::shared_ptr<Allele>> nix;
    return nix;};
  //  virtual std::vector<std::shared_ptr<Allele>> translateTo8d(){
  //    std::vector<std::shared_ptr<Allele>> nix;
  //    return nix;};  

 private:
  static FileGTog fileGTog;
};

class Allele6d : public Allele{
  
 public:
  explicit Allele6d(const std::string in_code,
		    const codePrecision in_precision,
		    const codePrecision in_wantedPrecision,
		    const double in_frequency)
    : Allele(in_code, in_precision, in_wantedPrecision, in_frequency){}
  explicit Allele6d(const std::string in_code,
		    const double in_frequency)
    : Allele(in_code, codePrecision::sixDigit, codePrecision::sixDigit, in_frequency){}
  virtual std::shared_ptr<Allele> create(const std::string in_code,
					 const codePrecision in_precision,
					 const codePrecision in_wantedPrecision,
					 const double in_frequency)
    {
      std::shared_ptr<Allele> pAllele = std::make_shared<Allele6d> (in_code,
								    in_precision,
								    in_wantedPrecision,
								    in_frequency);
      return pAllele;
    }
  
  virtual std::vector<std::shared_ptr<Allele>> translateTog();
  virtual std::vector<std::shared_ptr<Allele>> translateTo4d();
  virtual std::vector<std::shared_ptr<Allele>> translateToG(const FileAllelesTogOrG & whichH1File);
  virtual std::vector<std::shared_ptr<Allele>> translateTo6d(){
    std::vector<std::shared_ptr<Allele>> nix;
    return nix;};
  //  virtual std::vector<std::shared_ptr<Allele>> translateTo8d(){
  //    std::vector<std::shared_ptr<Allele>> nix;
  //  return nix;};  
  
 private:
};

class Allele8d : public Allele{

 public:
  explicit Allele8d(const std::string in_code,
		    const codePrecision in_precision,
		    const codePrecision in_wantedPrecision,
		    const double in_frequency)
    : Allele(in_code, in_precision, in_wantedPrecision, in_frequency){}
  virtual std::shared_ptr<Allele> create(const std::string in_code,
					 const codePrecision in_precision,
					 const codePrecision in_wantedPrecision,
					 const double in_frequency)
    {
      std::shared_ptr<Allele> pAllele = std::make_shared<Allele8d> (in_code,
								    in_precision,
								    in_wantedPrecision,
								    in_frequency);
      return pAllele;
    }

  virtual std::vector<std::shared_ptr<Allele>> translateTog();
  virtual std::vector<std::shared_ptr<Allele>> translateTo4d();
  virtual std::vector<std::shared_ptr<Allele>> translateToG(const FileAllelesTogOrG & whichH1File);
  virtual std::vector<std::shared_ptr<Allele>> translateTo6d(){
    std::vector<std::shared_ptr<Allele>> nix;
    return nix;};
  //  virtual std::vector<std::shared_ptr<Allele>> translateTo8d(){
  //  std::vector<std::shared_ptr<Allele>> nix;
  //  return nix;};  
  
 private:
};

#endif
