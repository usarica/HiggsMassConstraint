#ifndef SCALAR_PDF_FACTORY_HVV_FAST
#define SCALAR_PDF_FACTORY_HVV_FAST

#include <ZZMatrixElement/MELA/interface/ScalarPdfFactory_HVV.h>
#include "RooSpinZero_7DComplex_withAccep_HVV_fast.h"
#include "TGraph.h"

class ScalarPdfFactory_HVV_fast : public ScalarPdfFactory_HVV{
public:

  ScalarPdfFactory_HVV_fast(RooSpin::modelMeasurables measurables_, bool acceptance_=false, RooSpin::VdecayType V1decay_=RooSpin::kVdecayType_Zll, RooSpin::VdecayType V2decay_=RooSpin::kVdecayType_Zll, Bool_t OnshellH_=true);
  ScalarPdfFactory_HVV_fast(RooSpin::modelMeasurables measurables_, double gRatio_[4][8], double gZGsRatio_[4][1], double gGsGsRatio_[3][1], bool pmf_applied_=false, bool acceptance_=false, RooSpin::VdecayType V1decay_=RooSpin::kVdecayType_Zll, RooSpin::VdecayType V2decay_=RooSpin::kVdecayType_Zll, Bool_t OnshellH_=true);
  ~ScalarPdfFactory_HVV_fast(){}

  void setFastPDFOverallIntegration(TGraph* tg);

protected:

  TGraph* tg_int;

  virtual void initPDF();

};


#endif



