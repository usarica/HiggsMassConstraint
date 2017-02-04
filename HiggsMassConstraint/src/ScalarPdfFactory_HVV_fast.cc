#include "ScalarPdfFactory_HVV_fast.h"


ScalarPdfFactory_HVV_fast::ScalarPdfFactory_HVV_fast(RooSpinZero::modelMeasurables measurables_, bool acceptance_, RooSpin::VdecayType V1decay_, RooSpin::VdecayType V2decay_, Bool_t OnshellH_) :
ScalarPdfFactory_HVV(measurables_, acceptance_, V1decay_, V2decay_, OnshellH_),
tg_int(0)
{
  destroyPDF();
  initPDF();
}
ScalarPdfFactory_HVV_fast::ScalarPdfFactory_HVV_fast(RooSpinZero::modelMeasurables measurables_, double gRatio_[4][8], double gZGsRatio_[4][1], double gGsGsRatio_[3][1], bool pmf_applied_, bool acceptance_, RooSpin::VdecayType V1decay_, RooSpin::VdecayType V2decay_, Bool_t OnshellH_) :
ScalarPdfFactory_HVV(measurables_, gRatio_, gZGsRatio_, gGsGsRatio_, pmf_applied_, acceptance_, V1decay_, V2decay_, OnshellH_),
tg_int(0)
{
  destroyPDF();
  initPDF();
}

void ScalarPdfFactory_HVV_fast::initPDF(){
  PDF = new RooSpinZero_7DComplex_withAccep_HVV_fast(
    "PDF", "PDF",
    measurables,
    parameters,
    couplings,
    accepParams,
    V1decay,V2decay
    );
  PDF_base = (RooSpin*)PDF;
}

void ScalarPdfFactory_HVV_fast::setFastPDFOverallIntegration(TGraph* tg){
  tg_int=tg;
  ((RooSpinZero_7DComplex_withAccep_HVV_fast*)PDF)->setOverallIntegralGraph(*tg_int);
}