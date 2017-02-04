#include "RooHistSlicePdf.h" 
#include <cmath>
#include "Riostream.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
#include "TMath.h"
#include "TH3F.h"
#include "TAxis.h"
#include "RooDataHist.h"

using namespace TMath;
using namespace RooFit;
using namespace std;


RooHistSlicePdf::RooHistSlicePdf() :
RooAbsPdf(),
sliceVar("sliceVar", "sliceVar", this),
depList("depList", "List of dependents", this),
funcList("funcList", "List of functions", this),
coefList("coefList", "List of coefficients", this),
sliceVar_validitymin(0.),
sliceVar_validitymax(0.),
smoothAlgo(0),
smoothRegion(1.)
{}

RooHistSlicePdf::RooHistSlicePdf(
const char *name,
const char *title,
RooAbsReal& inSliceVar,
const RooArgList& inDepList,
const RooArgList& inFuncList,
const RooArgList& inCoefList,

Double_t inSliceVar_validitymin,
Double_t inSliceVar_validitymax,

Int_t inSmoothAlgo,
Double_t inSmoothRegion
) :
RooAbsPdf(name, title),
sliceVar("sliceVar", "sliceVar", this, inSliceVar),
depList("depList", "List of dependents", this),
funcList("funcList", "List of functions", this),
coefList("coefList", "List of coefficients", this),
sliceVar_validitymin(inSliceVar_validitymin),
sliceVar_validitymax(inSliceVar_validitymax),
smoothAlgo(inSmoothAlgo),
smoothRegion(inSmoothRegion)
{
  TIterator* coefIter = inCoefList.createIterator();
  RooAbsArg* coef;
  while ((coef = (RooAbsArg*)coefIter->Next())){
    if (!dynamic_cast<RooAbsReal*>(coef)){
      coutE(InputArguments) << "RooHistSlicePdf ERROR::RooHistSlicePdf(" << GetName() << ") coefficient " << coef->GetName() << " is not of type RooAbsReal" << endl;
      assert(0);
    }
    coefList.add(*coef);
  }
  delete coefIter;

  TIterator* funcIter = inFuncList.createIterator();
  RooAbsArg* func;
  while ((func = (RooAbsArg*)funcIter->Next())){
    if (!dynamic_cast<RooHistFunc*>(func)){
      coutE(InputArguments) << "RooHistSlicePdf ERROR::RooHistSlicePdf(" << GetName() << ") function " << func->GetName() << " is not of type RooHistFunc" << endl;
      assert(0);
    }
    funcList.add(*func);
  }
  delete funcIter;

  int nFuncs = (int)funcList.getSize();
  int nCoefs = (int)coefList.getSize();
  if (nFuncs%(2*nCoefs+1)==0){
    nSlices = nFuncs/(2*nCoefs+1);
  }
  else if (nFuncs%((int)(2*pow(nCoefs, 2)))==0){
    nSlices = nFuncs/((int)(2*pow(nCoefs, 2)));
  }
  else{
    coutE(InputArguments) << "RooHistSlicePdf ERROR::RooHistSlicePdf(" << GetName() << ") nFuncs " << nFuncs << " is not of divisible by 2*nCoefs+1 or 2*nCoefs**2. nCoefs " << nCoefs << endl;
    assert(0);
  }
  for (int mp=0; mp<nSlices; mp++) sliceIntegral.push_back(dynamic_cast<const RooHistFunc*>(funcList.at(mp))->analyticalIntegral(1000));
}

RooHistSlicePdf::RooHistSlicePdf(
  const RooHistSlicePdf& other,
  const char* name
  ) :
  RooAbsPdf(other, name),
  sliceVar("sliceVar", this, other.sliceVar),
  depList("depList", this, other.depList),
  funcList("funcList", this, other.funcList),
  coefList("funcList", this, other.coefList),
  nSlices(other.nSlices),
  sliceVar_validitymin(other.sliceVar_validitymin),
  sliceVar_validitymax(other.sliceVar_validitymax),
  sliceIntegral(other.sliceIntegral),
  smoothAlgo(other.smoothAlgo),
  smoothRegion(other.smoothRegion)
{}

Int_t RooHistSlicePdf::findNeighborBins() const{
  Double_t gridsize = (sliceVar_validitymax-sliceVar_validitymin)/((Double_t)(nSlices-1));
  Int_t bincode = (Int_t)(sliceVar-sliceVar_validitymin)/gridsize;
  return bincode;
}
Double_t RooHistSlicePdf::interpolateBin(Bool_t doIntegrate) const{
  Double_t result = 1e-100;

  Int_t bincode = findNeighborBins();
  Double_t gridsize = (sliceVar_validitymax-sliceVar_validitymin)/((Double_t)(nSlices-1));
  Int_t lowbin, highbin;
  if (bincode < 0){
    lowbin = 0; highbin = 1;
  }
  else if (bincode >= (nSlices - 1)){
    lowbin = nSlices - 2; highbin = nSlices - 1;
  }
  else{
    lowbin = bincode; highbin = lowbin + 1;
  }

  Double_t fdistance = (((sliceVar-sliceVar_validitymin)/gridsize) - ((Double_t)lowbin)) / ((Double_t)(highbin-lowbin));
  Double_t v1Nom = (!doIntegrate ? dynamic_cast<const RooHistFunc*>(funcList.at(lowbin))->getVal() : sliceIntegral[lowbin]);
  Double_t v2Nom = (!doIntegrate ? dynamic_cast<const RooHistFunc*>(funcList.at(highbin))->getVal() : sliceIntegral[highbin]);
  Double_t v1Inst = v1Nom;
  Double_t v2Inst = v2Nom;

  int nFuncs = funcList.getSize();
  int nCoefs = coefList.getSize();
  Bool_t fullVariation=!(nFuncs%(2*nCoefs+1)==0) && (nCoefs>0 && nFuncs%((int)(2*pow(nCoefs, 2)))==0);
  if (!fullVariation){
    for (Int_t ic=0; ic<nCoefs; ic++){
      RooAbsReal* coef = dynamic_cast<RooAbsReal*>(coefList.at(ic));
      Double_t coefVal = coef->getVal();

      int lowfunc, highfunc;
      lowfunc = lowbin+(2*ic+1)*nSlices;
      highfunc = highbin+(2*ic+1)*nSlices;
      Double_t v1Dn = (!doIntegrate ? dynamic_cast<const RooHistFunc*>(funcList.at(lowfunc))->getVal() : sliceIntegral[lowfunc]);
      Double_t v2Dn = (!doIntegrate ? dynamic_cast<const RooHistFunc*>(funcList.at(highfunc))->getVal() : sliceIntegral[highfunc]);
      lowfunc = lowbin+(2*ic+2)*nSlices;
      highfunc = highbin+(2*ic+2)*nSlices;
      Double_t v1Up = (!doIntegrate ? dynamic_cast<const RooHistFunc*>(funcList.at(lowfunc))->getVal() : sliceIntegral[lowfunc]);
      Double_t v2Up = (!doIntegrate ? dynamic_cast<const RooHistFunc*>(funcList.at(highfunc))->getVal() : sliceIntegral[highfunc]);

      v1Inst += interpolateVariation(coefVal, v1Nom, v1Up, v1Dn);
      v2Inst += interpolateVariation(coefVal, v2Nom, v2Up, v2Dn);
    }
  }
  else{
    coutE(InputArguments) << "RooHistSlicePdf::interpolateBin: Full variation is not supported at the moment." << endl;
    assert(0);
  }

  result = (1.-fdistance)*v1Inst + fdistance*v2Inst;
  return result;
}
Double_t RooHistSlicePdf::interpolateVariation(Double_t theta, Double_t valueCenter, Double_t valueHigh, Double_t valueLow) const {
  if (smoothAlgo<0) return 0;
  else{
    if (fabs(theta)>=smoothRegion) return theta * (theta > 0 ? valueHigh - valueCenter : valueCenter - valueLow);

    Double_t c_up  = 0;
    Double_t c_dn  = 0;
    Double_t c_cen = 0;
    Double_t addVal = 0;

    if (smoothAlgo != 1){
      // Quadratic interpolation null at zero and continuous at boundaries but not smooth at boundaries
      c_up  = +theta * (smoothRegion + theta) / (2 * smoothRegion);
      c_dn  = -theta * (smoothRegion - theta) / (2 * smoothRegion);
      c_cen = -theta * theta / smoothRegion;
    }
    else{
      // Quadratic interpolation that is everywhere differentiable but not null at zero
      c_up  = (smoothRegion + theta) * (smoothRegion + theta) / (4 * smoothRegion);
      c_dn  = (smoothRegion - theta) * (smoothRegion - theta) / (4 * smoothRegion);
      c_cen = -c_up - c_dn;
    }
    addVal = c_up * valueHigh + c_dn * valueLow + c_cen * valueCenter;
    return addVal;
  }
}

Double_t RooHistSlicePdf::evaluate() const{
  Double_t value = interpolateBin(false);
  if (value<=0.) return 1e-15;
  return value;
}
Int_t RooHistSlicePdf::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const{
  RooArgSet deps;
  for (Int_t ix=0; ix<depList.getSize(); ix++){
    RooRealVar* depvar=dynamic_cast<RooRealVar*>(depList.at(ix));
    if (depvar!=0) deps.add(*depvar);
  }
  if (matchArgs(allVars, analVars, deps)) return 1;
  else return 0;
}
Double_t RooHistSlicePdf::analyticalIntegral(Int_t code, const char* rangeName) const{
  switch (code){
  case 1:
  {
          Double_t value = interpolateBin(true);
          if (value<=0.) value = 1e-10;
          return value;
  }
  default:
    cerr << "getAnalyticalIntegral failed, so analytical integration did not complete!" << endl;
    assert(0);
  }
}
