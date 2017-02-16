#ifndef ROOSPINZERO_7DCOMPLEX_WITHACCEP_HVV_FAST
#define ROOSPINZERO_7DCOMPLEX_WITHACCEP_HVV_FAST

#include <vector>
#include <utility>
#include <ZZMatrixElement/MELA/interface/RooSpinZero_7DComplex_withAccep_HVV.h>
#include "TGraph.h"




class RooSpinZero_7DComplex_withAccep_HVV_fast : public RooSpinZero_7DComplex_withAccep_HVV{
public:

  RooSpinZero_7DComplex_withAccep_HVV_fast(){}
  RooSpinZero_7DComplex_withAccep_HVV_fast(
    const char *name, const char *title,
    modelMeasurables _measurables,
    modelParameters _parameters,
    modelCouplings _couplings,
    accepParameters _accepParams,
    RooSpin::VdecayType _Vdecay1=RooSpin::kVdecayType_Zll, RooSpin::VdecayType _Vdecay2=RooSpin::kVdecayType_Zll
    );
  RooSpinZero_7DComplex_withAccep_HVV_fast(const RooSpinZero_7DComplex_withAccep_HVV_fast& other, const char* name=0);
  virtual TObject* clone(const char* newname) const { return new RooSpinZero_7DComplex_withAccep_HVV_fast(*this, newname); }
  inline virtual ~RooSpinZero_7DComplex_withAccep_HVV_fast(){}

  Double_t evaluate() const;
  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const;

  void setOverallIntegralGraph(TGraph& tg){ tg_integral=tg; }

protected:

  TGraph tg_integral;

  Double_t evalOverallIntegral()const;

};

#endif
