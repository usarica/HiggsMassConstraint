#include "RooNCSplinePdf_1D_fast.h" 
#include <cmath>
#include "TMath.h"
#include "Riostream.h" 
#include "RooAbsReal.h" 

using namespace TMath;
using namespace RooFit;
using namespace std;


RooNCSplinePdf_1D_fast::RooNCSplinePdf_1D_fast() :
RooNCSplinePdfCore()
{}

RooNCSplinePdf_1D_fast::RooNCSplinePdf_1D_fast(
const char* name,
const char* title,
RooAbsReal& inXVar,
std::vector<T>& inXList,
std::vector<T>& inFcnList
) :
RooNCSplinePdfCore(name, title, inXVar, inXList),
FcnList(inFcnList)
{
  if (npointsX()>1){
    int npoints;

    vector<vector<RooNCSplinePdfCore::T>> xA; getKappas(kappaX, 0); getAArray(kappaX, xA);
    npoints=kappaX.size();
    TMatrix_t xAtrans(npoints, npoints);
    for (int i=0; i<npoints; i++){ for (int j=0; j<npoints; j++){ xAtrans[i][j]=xA.at(i).at(j); } }
    Double_t det=0;
    TMatrix_t xAinv = xAtrans.Invert(&det);
    if (det==0.){
      coutE(InputArguments) << "RooNCSplinePdf_1D::interpolateFcn: Matrix xA could not be inverted. Something is wrong with the x coordinates of points!" << endl;
      assert(0);
    }

    coefficients = getCoefficientsAlongDirection(kappaX, xAinv, FcnList, -1);
  }
  else assert(0);

  emptyFcnList();
}

RooNCSplinePdf_1D_fast::RooNCSplinePdf_1D_fast(
  const RooNCSplinePdf_1D_fast& other,
  const char* name
  ) :
  RooNCSplinePdfCore(other, name),
  FcnList(other.FcnList),
  kappaX(other.kappaX),
  coefficients(other.coefficients)
{}

RooNCSplinePdfCore::T RooNCSplinePdf_1D_fast::interpolateFcn(Int_t code, const char* rangeName)const{
  DefaultAccumulator<RooNCSplinePdfCore::T> res=RooNCSplinePdfCore::T(0);

  // Get bins
  Int_t xbin=-1, xbinmin=-1, xbinmax=-1;
  RooNCSplinePdfCore::T tx=0, txmin=0, txmax=0;
  if (code==0 || code%2!=0){ // Case to just compute the value at x
    xbin = getWhichBin(theXVar, 0);
    tx = getTVar(kappaX, theXVar, xbin, 0);
  }
  else{ // Case to integrate along x
    xbinmin = getWhichBin(theXVar.min(rangeName), 0);
    txmin = getTVar(kappaX, theXVar.min(rangeName), xbinmin, 0);
    xbinmax = getWhichBin(theXVar.max(rangeName), 0);
    txmax = getTVar(kappaX, theXVar.max(rangeName), xbinmax, 0);
  }

  int nxbins = (int)coefficients.size();
  for (int ix=0; ix<nxbins; ix++){
    if (
      (xbin>=0 && ix!=xbin)
      ||
      (xbinmin>=0 && xbinmax>=xbinmin && !(xbinmin<=ix && ix<=xbinmax))
      ) continue;

    RooNCSplinePdfCore::T txlow=0, txhigh=1;
    if (code>0 && code%2==0){
      if (ix==xbinmin) txlow=txmin;
      if (ix==xbinmax) txhigh=txmax;
    }
    else txhigh=tx;

    // Get the x coefficients at bin ix and evaluate value of spline at x
    res += evalSplineSegment(coefficients.at(ix), kappaX.at(ix), txhigh, txlow, (code>0 && code%2==0));
  }

  return res.sum();
}

void RooNCSplinePdf_1D_fast::getKappas(vector<RooNCSplinePdfCore::T>& kappas, const Int_t /*whichDirection*/)const{
  kappas.clear();
  RooNCSplinePdfCore::T kappa=1;

  Int_t npoints;
  vector<RooNCSplinePdfCore::T> const* coord;
  npoints=npointsX();
  coord=&XList;

  for (Int_t j=0; j<npoints-1; j++){
    RooNCSplinePdfCore::T val_j = coord->at(j);
    RooNCSplinePdfCore::T val_jpo = coord->at(j+1);
    kappa = 1./(val_jpo-val_j);
    kappas.push_back(kappa);
  }
  kappas.push_back(kappa); // Push the same kappa_(N-1)=kappa_(N-2) at the end point
}
Int_t RooNCSplinePdf_1D_fast::getWhichBin(const RooNCSplinePdfCore::T& val, const Int_t /*whichDirection*/)const{
  Int_t bin=-1;
  RooNCSplinePdfCore::T valj, valjpo;
  Int_t npoints;
  vector<RooNCSplinePdfCore::T> const* coord;
  coord=&XList;
  npoints=npointsX();

  if (npoints<=1) bin=0;
  else{
    valjpo = coord->at(0);
    for (Int_t j=0; j<npoints-1; j++){
      valj = coord->at(j);
      valjpo = coord->at(j+1);
      if (val<valjpo && val>=valj){ bin=j; break; }
    }
    if (bin==-1 && val>=valjpo) bin=npoints-2;
    else if (bin==-1) bin=0;
  }

  return bin;
}
RooNCSplinePdfCore::T RooNCSplinePdf_1D_fast::getTVar(const vector<RooNCSplinePdfCore::T>& kappas, const RooNCSplinePdfCore::T& val, const Int_t& bin, const Int_t /*whichDirection*/)const{
  RooNCSplinePdfCore::T K;
  K=kappas.at(bin);
  return (val-XList.at(bin))*K;
}

Double_t RooNCSplinePdf_1D_fast::evaluate() const{
  Double_t value = interpolateFcn(0);
  if (value<=0.) return 1e-15;
  return value;
}
Int_t RooNCSplinePdf_1D_fast::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const{
  Int_t code=0;
  if (dynamic_cast<RooRealVar*>(theXVar.absArg())!=0){
    if (matchArgs(allVars, analVars, theXVar)) code=2;
  }
  return code;
}
Double_t RooNCSplinePdf_1D_fast::analyticalIntegral(Int_t code, const char* rangeName) const{
  Double_t value = interpolateFcn(code, rangeName);
  if (value<=0.) value = 1e-10;
  return value;
}

ClassImp(RooNCSplinePdf_1D_fast)
