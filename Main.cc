#include<memory>

#include "DataProcessing.h"

int main(){

  std::unique_ptr<DataProcessing> pDataProcessing(new DataProcessingDKMS("reports.txt"));
  pDataProcessing->dataProcessing();



  return 0;
}
