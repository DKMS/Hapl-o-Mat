#include <memory>

#include "DataProcessing.h"
#include "Allele.h"
#include "Locus.h"
#include "Typedefs.h"
#include "Phenotype.h"
#include "Haplotype.h"

int main(){
  
  strVec_t lociToDo;
  lociToDo.push_back("A");
  lociToDo.push_back("C");
  lociToDo.push_back("None");
  lociToDo.push_back("DRB1");
  lociToDo.push_back("DQB1");

  //  std::unique_ptr<DataProcessing> pDataProcessing (new GLDataProcessing("reports.pull", "reports.glid", lociToDo,   Allele::codePrecision::fourDigit));
  std::unique_ptr<DataProcessing> pDataProcessing (new DKMSDataProcessing("reports.txt", Allele::codePrecision::g,.0001));

  PhenotypeList pList;
  HaplotypeList hList;

  pDataProcessing->dataProcessing(pList, hList);

}
