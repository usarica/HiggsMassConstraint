#include "RooSpinZero_7DComplex_withAccep_HVV_fast.h"


RooSpinZero_7DComplex_withAccep_HVV_fast::RooSpinZero_7DComplex_withAccep_HVV_fast(
  const char *name, const char *title,
  modelMeasurables _measurables,
  modelParameters _parameters,
  modelCouplings _couplings,
  accepParameters _accepParams,
  RooSpin::VdecayType _Vdecay1, RooSpin::VdecayType _Vdecay2
  ) : RooSpinZero_7DComplex_withAccep_HVV(
  name, title,
  _measurables,
  _parameters,
  _couplings,
  _accepParams,
  _Vdecay1, _Vdecay2
  )
{}

RooSpinZero_7DComplex_withAccep_HVV_fast::RooSpinZero_7DComplex_withAccep_HVV_fast(
  const RooSpinZero_7DComplex_withAccep_HVV_fast& other, const char* name
  ) : RooSpinZero_7DComplex_withAccep_HVV(other, name),
  tg_integral(other.tg_integral)
{}

Double_t RooSpinZero_7DComplex_withAccep_HVV_fast::evaluate() const{
  Double_t value = RooSpinZero_7DComplex_withAccep_HVV::evaluate();
  if (
    value>1e-15
    && intCodeStart%prime_m1==0
    && intCodeStart%prime_m2==0
    && intCodeStart%prime_h1==0
    && intCodeStart%prime_h2==0
    && intCodeStart%prime_hs==0
    && intCodeStart%prime_Phi==0
    && intCodeStart%prime_Phi1==0
    ) value /= evalOverallIntegral();
  //std::cout << "value/code=" << value << " / " << intCodeStart << std::endl;
  return value;
}

Int_t RooSpinZero_7DComplex_withAccep_HVV_fast::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const{
  Int_t code = intCodeStart;
  if (checkFundamentalType(m1)){ if (matchArgs(allVars, analVars, m1) || Vdecay1==RooSpin::kVdecayType_GammaOnshell) code *= prime_m1; }
  if (checkFundamentalType(m2)){ if (matchArgs(allVars, analVars, m2) || Vdecay2==RooSpin::kVdecayType_GammaOnshell) code *= prime_m2; }
  if (checkFundamentalType(h1)){ if (matchArgs(allVars, analVars, h1) || Vdecay1==RooSpin::kVdecayType_GammaOnshell) code *= prime_h1; }
  if (checkFundamentalType(h2)){ if (matchArgs(allVars, analVars, h2) || Vdecay2==RooSpin::kVdecayType_GammaOnshell) code *= prime_h2; }
  if (checkFundamentalType(hs)){ if (matchArgs(allVars, analVars, hs)) code *= prime_hs; }
  if (checkFundamentalType(Phi)){ if (matchArgs(allVars, analVars, Phi) || Vdecay1==RooSpin::kVdecayType_GammaOnshell || Vdecay2==RooSpin::kVdecayType_GammaOnshell) code *= prime_Phi; }
  if (checkFundamentalType(Phi1)){ if (matchArgs(allVars, analVars, Phi1) || (Vdecay1==RooSpin::kVdecayType_GammaOnshell && Vdecay2==RooSpin::kVdecayType_GammaOnshell)) code *= prime_Phi1; }
  if (code==1) code=0;
  return code;
}
Double_t RooSpinZero_7DComplex_withAccep_HVV_fast::analyticalIntegral(Int_t code, const char* rangeName) const{
  //std::cout << "int/code=" << intCodeStart << " / " << code << std::endl;
  if (
    intCodeStart%prime_m1==0
    && intCodeStart%prime_m2==0
    && intCodeStart%prime_h1==0
    && intCodeStart%prime_h2==0
    && intCodeStart%prime_hs==0
    && intCodeStart%prime_Phi==0
    && intCodeStart%prime_Phi1==0
    ) return 1.;
  else return RooSpinZero_7DComplex_withAccep_HVV::analyticalIntegral(code, rangeName);
}
Double_t RooSpinZero_7DComplex_withAccep_HVV_fast::evalOverallIntegral() const{
  return tg_integral.Eval(m12);
}
