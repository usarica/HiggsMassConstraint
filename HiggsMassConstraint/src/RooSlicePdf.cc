#include "RooSlicePdf.h" 
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


RooSlicePdf::RooSlicePdf() :
RooAbsPdf(),
sliceVar("sliceVar", "sliceVar", this),
depList("depList", "List of dependents", this),
sliceList("sliceList", "List of slices", this),
funcList("funcList", "List of functions", this),
coefList("coefList", "List of coefficients", this),
smoothAlgo(0),
smoothRegion(1.)
{}

RooSlicePdf::RooSlicePdf(
const char* name,
const char* title,
RooAbsReal& inSliceVar,
const RooArgList& inDepList,
const RooArgList& inSliceList,
const RooArgList& inFuncList,
const RooArgList& inCoefList,

Int_t inSmoothAlgo,
Double_t inSmoothRegion
) :
RooAbsPdf(name, title),
sliceVar("sliceVar", "sliceVar", this, inSliceVar),
depList("depList", "List of dependents", this),
sliceList("sliceList", "List of slices", this),
funcList("funcList", "List of functions", this),
coefList("coefList", "List of coefficients", this),
smoothAlgo(inSmoothAlgo),
smoothRegion(inSmoothRegion)
{

  setListProxyFromArgList(inDepList, depList);
  setListProxyFromArgList(inCoefList, coefList);
  setListProxyFromArgList(inSliceList, sliceList);
  setListProxyFromArgList(inFuncList, funcList);
  
  nSlices=sliceList.getSize();
  int nFuncs = (int)funcList.getSize();
  int nCoefs = (int)coefList.getSize();
  if (nFuncs!=(2*nCoefs+1)*nSlices){
    coutE(InputArguments) << "RooSlicePdf ERROR::RooSlicePdf(" << GetName() << ") nFuncs " << nFuncs << " nSlices*(2*nCoefs+1). nCoefs " << nCoefs << " nSlices " << nSlices << endl;
    assert(0);
  }
  for (int mp=0; mp<nFuncs; mp++) sliceIntegral.push_back(0);
}

RooSlicePdf::RooSlicePdf(
  const RooSlicePdf& other,
  const char* name
  ) :
  RooAbsPdf(other, name),
  sliceVar("sliceVar", this, other.sliceVar),
  depList("depList", this, other.depList),
  sliceList("sliceList", this, other.sliceList),
  funcList("funcList", this, other.funcList),
  coefList("coefList", this, other.coefList),
  nSlices(other.nSlices),
  sliceIntegral(other.sliceIntegral),
  smoothAlgo(other.smoothAlgo),
  smoothRegion(other.smoothRegion)
{}

Int_t RooSlicePdf::findNeighborBins() const{
  Int_t bin=-1;
  for (int is=0; is<nSlices-1; is++){
    if (sliceVar>=dynamic_cast<const RooAbsReal*>(sliceList.at(is))->getVal() && sliceVar<dynamic_cast<const RooAbsReal*>(sliceList.at(is+1))->getVal()){
      bin=is;
      break;
    }
  }
  if (bin==-1 && sliceVar>=dynamic_cast<const RooAbsReal*>(sliceList.at(nSlices-1))->getVal()) bin=nSlices-2;
  else if (bin==-1) bin=0;
  return bin;
}
Double_t RooSlicePdf::findBinWidth(Int_t bin) const{
  Double_t res=0;
  res = dynamic_cast<const RooAbsReal*>(sliceList.at(bin+1))->getVal() - dynamic_cast<const RooAbsReal*>(sliceList.at(bin))->getVal();
  return res;
}
Double_t RooSlicePdf::interpolateBin(Int_t intCode) const{
  Bool_t useSliceIntegral=(intCode>0 && intCode%2==0);
  //Bool_t integrateOverSlice=(intCode>0 && intCode%3==0);

  Int_t bincode = findNeighborBins();
  Double_t gridsize = findBinWidth(bincode);
  Int_t lowbin, highbin;
  lowbin = bincode; highbin = lowbin + 1;

  Double_t fdistance = (sliceVar - dynamic_cast<const RooAbsReal*>(sliceList.at(lowbin))->getVal())/gridsize;
  Double_t v1Nom = (!useSliceIntegral ? dynamic_cast<const RooAbsReal*>(funcList.at(lowbin))->getVal() : sliceIntegral[lowbin]);
  Double_t v2Nom = (!useSliceIntegral ? dynamic_cast<const RooAbsReal*>(funcList.at(highbin))->getVal() : sliceIntegral[highbin]);
  Double_t v1Inst = v1Nom;
  Double_t v2Inst = v2Nom;
  cout << "RooSlicePdf::interpolateBin: sliceVar=" << sliceVar << " is in bin " << bincode << " with fdistance=" << fdistance << endl;
  cout << "RooSlicePdf::interpolateBin: v1Inst=" << v1Inst << " v2Inst=" << v2Inst << endl;

  int nFuncs = funcList.getSize();
  int nCoefs = coefList.getSize();
  cout << "nFuncs=" << nFuncs << endl;
  cout << "nCoefs=" << nCoefs << endl;
  for (Int_t ic=0; ic<nCoefs; ic++){
    RooAbsReal* coef = dynamic_cast<RooAbsReal*>(coefList.at(ic));
    Double_t coefVal = coef->getVal();

    int lowfunc, highfunc;
    lowfunc = lowbin+(2*ic+1)*nSlices;
    highfunc = highbin+(2*ic+1)*nSlices;
    Double_t v1Dn = (!useSliceIntegral ? dynamic_cast<const RooAbsReal*>(funcList.at(lowfunc))->getVal() : sliceIntegral[lowfunc]);
    Double_t v2Dn = (!useSliceIntegral ? dynamic_cast<const RooAbsReal*>(funcList.at(highfunc))->getVal() : sliceIntegral[highfunc]);
    lowfunc = lowbin+(2*ic+2)*nSlices;
    highfunc = highbin+(2*ic+2)*nSlices;
    Double_t v1Up = (!useSliceIntegral ? dynamic_cast<const RooAbsReal*>(funcList.at(lowfunc))->getVal() : sliceIntegral[lowfunc]);
    Double_t v2Up = (!useSliceIntegral ? dynamic_cast<const RooAbsReal*>(funcList.at(highfunc))->getVal() : sliceIntegral[highfunc]);

    v1Inst += interpolateVariation(coefVal, v1Nom, v1Up, v1Dn);
    v2Inst += interpolateVariation(coefVal, v2Nom, v2Up, v2Dn);
  }

  Double_t A=v1Inst;
  Double_t B=v2Inst-v1Inst;
  Double_t result = A + B*fdistance;
  return result;
}
Double_t RooSlicePdf::interpolateVariation(Double_t theta, Double_t valueCenter, Double_t valueHigh, Double_t valueLow) const {
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

Double_t RooSlicePdf::evaluate() const{
  Double_t value = interpolateBin(0);
  if (value<=0.) return 1e-15;
  return value;
}
Int_t RooSlicePdf::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const{
  TIterator* parIter = allVars.createIterator();
  RooAbsArg* thePar;
  while ((thePar = (RooAbsArg*)parIter->Next())) cout << "allVars has " << thePar->GetName() << endl;
  cout << endl;
  delete parIter;

  parIter = analVars.createIterator();
  while ((thePar = (RooAbsArg*)parIter->Next())) cout << "analVars has " << thePar->GetName() << endl;
  cout << endl;
  delete parIter;

  return 0;
}
Double_t RooSlicePdf::analyticalIntegral(Int_t code, const char* rangeName) const{
  Double_t value = interpolateBin(0);
  if (value<=0.) return 1e-10;
  return value;
}

void RooSlicePdf::setListProxyFromArgList(const RooArgList& argList, RooListProxy& proxyList){
  TIterator* argIter = argList.createIterator();
  RooAbsArg* arg;
  while ((arg = (RooAbsArg*)argIter->Next())){
    if (!dynamic_cast<RooAbsReal*>(arg)){
      coutE(InputArguments) << "RooSlicePdf ERROR::RooSlicePdf(" << GetName() << ") argument " << arg->GetName() << " is not of type RooAbsReal" << endl;
      assert(0);
    }
    proxyList.add(*arg);
  }
  delete argIter;
}


ClassImp(RooSlicePdf)
