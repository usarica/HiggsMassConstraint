#include "RooNCSplinePdf_1D_fast.h" 
#include <cmath>
#include "Riostream.h" 
#include "RooAbsReal.h" 
#include "RooConstVar.h" 
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


RooNCSplinePdf_1D_fast::RooNCSplinePdf_1D_fast() :
RooNCSplinePdf_1D()
{}

RooNCSplinePdf_1D_fast::RooNCSplinePdf_1D_fast(
const char* name,
const char* title,
RooAbsReal& inXVar,
const RooArgList& inXList,
const RooArgList& inFcnList
) :
RooNCSplinePdf_1D(name, title, inXVar, inXList, inFcnList)
{
  TIterator* coefIter = inXList.createIterator();
  RooAbsArg* coef;
  while ((coef = (RooAbsArg*)coefIter->Next())){
    if (!dynamic_cast<RooConstVar*>(coef)){
      coutE(InputArguments) << "RooNCSplinePdf_1D_fast ERROR::RooNCSplinePdf_1D_fast(" << GetName() << ") X variable " << coef->GetName() << " is not of type RooConstVar" << endl;
      assert(0);
    }
  }
  delete coefIter;

  coefIter = inFcnList.createIterator();
  while ((coef = (RooAbsArg*)coefIter->Next())){
    if (!dynamic_cast<RooConstVar*>(coef)){
      coutE(InputArguments) << "RooNCSplinePdf_1D_fast ERROR::RooNCSplinePdf_1D_fast(" << GetName() << ") function variable " << coef->GetName() << " is not of type RooConstVar" << endl;
      assert(0);
    }
  }
  delete coefIter;

  if (npointsX>1){
    int npoints;

    vector<vector<Double_t>> xA; getKappa(kappaX, 0); getAArray(kappaX, xA);
    npoints=kappaX.size();
    TMatrixD xAtrans(npoints, npoints);
    for (int i=0; i<npoints; i++){ for (int j=0; j<npoints; j++){ xAtrans[i][j]=xA.at(i).at(j); } }
    Double_t det=0;
    TMatrixD xAinv = xAtrans.Invert(&det);
    if (det==0.){
      coutE(InputArguments) << "RooNCSplinePdf_1D::interpolateFcn: Matrix xA could not be inverted. Something is wrong with the x coordinates of points!" << endl;
      assert(0);
    }

    vector<Double_t> fcnList;
    for (int bin=0; bin<npointsX; bin++) fcnList.push_back((dynamic_cast<RooAbsReal*>(FcnList.at(bin)))->getVal());
    coefficients = getCoefficientsAlongDirection(kappaX, xAinv, fcnList, -1);
  }
  else assert(0);
}

RooNCSplinePdf_1D_fast::RooNCSplinePdf_1D_fast(
  const RooNCSplinePdf_1D_fast& other,
  const char* name
  ) : RooNCSplinePdf_1D(other,name),
  kappaX(other.kappaX)
{
  for (unsigned int ibin=0; ibin<other.coefficients.size(); ibin++){
    vector<Double_t> coef;
    for (unsigned int ic=0; ic<other.coefficients.at(ibin).size(); ic++) coef.push_back(other.coefficients.at(ibin).at(ic));
    coefficients.push_back(coef);
  }
}

Double_t RooNCSplinePdf_1D_fast::interpolateFcn(Int_t code, const char* rangeName)const{
  Double_t res = 0;

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

  int nxbins = (int)coefficients.size();
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

    // Get the x coefficients at bin ix and evaluate value of spline at x
    res += evalSplineSegment(coefficients.at(ix), kappaX.at(ix), txhigh, txlow, (code>0 && code%2==0));
  }

  return res;
}

Double_t RooNCSplinePdf_1D_fast::evaluate() const{
  Double_t value = interpolateFcn(0);
  if (value<=0.) return 1e-15;
  return value;
}

Double_t RooNCSplinePdf_1D_fast::analyticalIntegral(Int_t code, const char* rangeName) const{
  Double_t value = interpolateFcn(code, rangeName);
  if (value<=0.) value = 1e-10;
  return value;
}

ClassImp(RooNCSplinePdf_1D_fast)
