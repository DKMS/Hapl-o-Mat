#include<memory>

#include "Typedefs.h"
#include "DataProcessing.h"

int main(){

  strVec_t lociToDo;
  lociToDo.push_back("A");
  lociToDo.push_back("NONE");
  lociToDo.push_back("B");
  lociToDo.push_back("NONE");
  lociToDo.push_back("NONE");

  std::unique_ptr<DataProcessing> pDataProcessing(new DataProcessingGL("reports.pull", lociToDo));
  pDataProcessing->dataProcessing();


    //  std::unique_ptr<DataProcessing> pDataProcessing(new DataProcessingGL("reports.txt"));
    //  pDataProcessing->dataProcessing();




  return 0;
}
