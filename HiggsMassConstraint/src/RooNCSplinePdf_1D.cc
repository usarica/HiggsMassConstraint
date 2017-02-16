#include "RooNCSplinePdf_1D.h" 
#include <cmath>
#include "Riostream.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
#include "TH3F.h"
#include "TAxis.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "RooDataHist.h"

using namespace TMath;
using namespace RooFit;
using namespace std;


RooNCSplinePdf_1D::RooNCSplinePdf_1D() :
RooNCSplinePdfCore(),
theXVar("theXVar", "theXVar", this),
XList("XList", "List of X coordinates", this),
FcnList("FcnList", "List of function coordinates", this),
npointsX(0)
{}

RooNCSplinePdf_1D::RooNCSplinePdf_1D(
const char* name,
const char* title,
RooAbsReal& inXVar,
const RooArgList& inXList,
const RooArgList& inFcnList
) :
RooNCSplinePdfCore(name, title),
theXVar("theXVar", "theXVar", this, inXVar),
XList("XList", "List of X coordinates", this),
FcnList("FcnList", "List of function coordinates", this)
{
  TIterator* coefIter = inXList.createIterator();
  RooAbsArg* coef;
  while ((coef = (RooAbsArg*)coefIter->Next())){
    if (!dynamic_cast<RooAbsReal*>(coef)){
      coutE(InputArguments) << "RooNCSplinePdf_1D ERROR::RooNCSplinePdf_1D(" << GetName() << ") X variable " << coef->GetName() << " is not of type RooAbsReal" << endl;
      assert(0);
    }
    XList.add(*coef);
  }
  delete coefIter;

  npointsX = XList.getSize();
  if (npointsX<=1){
    coutE(InputArguments) << "RooNCSplinePdf_1D ERROR::RooNCSplinePdf_1D(" << GetName() << ") npointsX " << npointsX << "<=1 cannot have a spline." << endl;
    assert(0);
  }

  coefIter = inFcnList.createIterator();
  while ((coef = (RooAbsArg*)coefIter->Next())){
    if (!dynamic_cast<RooAbsReal*>(coef)){
      coutE(InputArguments) << "RooNCSplinePdf_1D ERROR::RooNCSplinePdf_1D(" << GetName() << ") Y variable " << coef->GetName() << " is not of type RooAbsReal" << endl;
      assert(0);
    }
    FcnList.add(*coef);
  }
  delete coefIter;
}

RooNCSplinePdf_1D::RooNCSplinePdf_1D(
  const RooNCSplinePdf_1D& other,
  const char* name
  ) :
  RooNCSplinePdfCore(other, name),
  theXVar("theXVar", this, other.theXVar),
  XList("XList", this, other.XList),
  FcnList("FcnList", this, other.FcnList),
  npointsX(other.npointsX)
{}

Double_t RooNCSplinePdf_1D::interpolateFcn(Int_t code, const char* rangeName)const{
  Double_t res=0;

  // Prepare A and kappa arrays for x coordinate
  int npoints;
  Double_t det;

  vector<vector<Double_t>> xA; vector<Double_t> kappaX; getKappa(kappaX, 0); getAArray(kappaX, xA);
  npoints=kappaX.size();
  TMatrixD xAtrans(npoints, npoints);
  for (int i=0; i<npoints; i++){ for (int j=0; j<npoints; j++){ xAtrans[i][j]=xA.at(i).at(j); } }
  det=0;
  TMatrixD xAinv = xAtrans.Invert(&det);
  if (det==0.){
    coutE(InputArguments) << "RooNCSplinePdf_1D::interpolateFcn: Matrix xA could not be inverted. Something is wrong with the x coordinates of points!" << endl;
    assert(0);
  }

  // Get bins
  Int_t xbin=-1, xbinmin=-1, xbinmax=-1;
  Double_t tx=0, txmin=0, txmax=0;
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

  // Get the grid of coefficients
  int nxbins=0;
  vector<Double_t> fcnList;
  for (int bin=0; bin<npointsX; bin++) fcnList.push_back((dynamic_cast<RooAbsReal*>(FcnList.at(bin)))->getVal());
  vector<vector<Double_t>> coefs = getCoefficientsAlongDirection(kappaX, xAinv, fcnList, -1);
  nxbins = (int)coefs.size();

  for (int ix=0; ix<nxbins; ix++){
    if (
      (xbin>=0 && ix!=xbin)
      ||
      (xbinmin>=0 && xbinmax>=xbinmin && !(xbinmin<=ix && ix<=xbinmax))
      ) continue;

    Double_t txlow=0, txhigh=1;
    if (code>0 && code%2==0){
      if (ix==xbinmin) txlow=txmin;
      if (ix==xbinmax) txhigh=txmax;
    }
    else txhigh=tx;

    // Get the x coefficients at bin ix
    vector<Double_t> xCoefs = coefs.at(ix);

    // Evaluate value of spline at x
    res += evalSplineSegment(xCoefs, txhigh, txlow, (code>0 && code%2==0));
  }

  return res;
}

void RooNCSplinePdf_1D::getKappa(vector<Double_t>& kappas, const Int_t /*whichDirection*/)const{
  kappas.clear();
  Double_t kappa;

  Int_t npoints;
  RooListProxy const* coord;
  npoints=npointsX;
  coord=&XList;

  for (int j=0; j<npoints-1; j++){
    Double_t val_j = (dynamic_cast<RooAbsReal*>(coord->at(j)))->getVal();
    Double_t val_jpo = (dynamic_cast<RooAbsReal*>(coord->at(j+1)))->getVal();
    kappa = 1./(val_jpo-val_j);
    kappas.push_back(kappa);
  }
  kappas.push_back(kappa); // Push the same kappa_(N-1)=kappa_(N-2) at the end point
}
Int_t RooNCSplinePdf_1D::getWhichBin(const Double_t& val, const Int_t /*whichDirection*/)const{
  Int_t bin=-1;
  Double_t valj, valjpo;
  int npoints;
  RooListProxy const* gridDir;
  gridDir=&XList;
  npoints=npointsX;

  if (npoints<=1) bin=0;
  else{
    valjpo = (dynamic_cast<RooAbsReal*>(gridDir->at(0)))->getVal();
    for (Int_t j=0; j<npoints-1; j++){
      valj = (dynamic_cast<RooAbsReal*>(gridDir->at(j)))->getVal();
      valjpo = (dynamic_cast<RooAbsReal*>(gridDir->at(j+1)))->getVal();
      if (val<valjpo && val>=valj){ bin=j; break; }
    }
    if (bin==-1 && val>=valjpo) bin=npoints-2;
    else if (bin==-1) bin=0;
  }

  return bin;
}
Double_t RooNCSplinePdf_1D::getTVar(const vector<Double_t>& kappas, const Double_t& val, const Int_t& bin, const Int_t /*whichDirection*/)const{
  Double_t K;
  K=kappas.at(bin);
  return (Double_t)((val-(dynamic_cast<RooAbsReal*>(XList.at(bin)))->getVal())*K);
}

Double_t RooNCSplinePdf_1D::evaluate() const{
  Double_t value = interpolateFcn(0);
  if (value<=0.) return 1e-15;
  return value;
}
Int_t RooNCSplinePdf_1D::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const{
  Int_t code=0;
  if (dynamic_cast<RooRealVar*>(theXVar.absArg())!=0){
    if (matchArgs(allVars, analVars, theXVar)) code=2;
  }
  return code;
}
Double_t RooNCSplinePdf_1D::analyticalIntegral(Int_t code, const char* rangeName) const{
  Double_t value = interpolateFcn(code, rangeName);
  if (value<=0.) value = 1e-10;
  return value;
}

ClassImp(RooNCSplinePdf_1D)
