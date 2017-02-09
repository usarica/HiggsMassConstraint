#include "RooNCSplinePdf_fast.h" 
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


RooNCSplinePdf_fast::RooNCSplinePdf_fast() :
RooNCSplinePdf()
{}

RooNCSplinePdf_fast::RooNCSplinePdf_fast(
const char* name,
const char* title,
RooAbsReal& inVar,
const RooArgList& inXList,
const RooArgList& inYList
) :
RooNCSplinePdf(name, title, inVar, inXList, inYList)
{
  TIterator* coefIter = inXList.createIterator();
  RooAbsArg* coef;
  while ((coef = (RooAbsArg*)coefIter->Next())){
    if (!dynamic_cast<RooConstVar*>(coef)){
      coutE(InputArguments) << "RooNCSplinePdf_fast ERROR::RooNCSplinePdf_fast(" << GetName() << ") X variable " << coef->GetName() << " is not of type RooConstVar" << endl;
      assert(0);
    }
  }
  delete coefIter;

  coefIter = inYList.createIterator();
  while ((coef = (RooAbsArg*)coefIter->Next())){
    if (!dynamic_cast<RooConstVar*>(coef)){
      coutE(InputArguments) << "RooNCSplinePdf_fast ERROR::RooNCSplinePdf_fast(" << GetName() << ") Y variable " << coef->GetName() << " is not of type RooConstVar" << endl;
      assert(0);
    }
  }
  delete coefIter;

  if (npoints>1){
    vector<Double_t> BArray;
    vector<vector<Double_t>> AArray;
    getKappa(kappas);
    getBArray(kappas, BArray);
    getAArray(kappas, AArray);

    TVectorD Btrans(npoints);
    TMatrixD Atrans(npoints, npoints);
    for (int i=0; i<npoints; i++){
      Btrans[i]=BArray.at(i);
      for (int j=0; j<npoints; j++){
        Atrans[i][j]=AArray.at(i).at(j);
      }
    }
    Double_t det=0;
    TMatrixD Ainv = Atrans.Invert(&det);
    if (det==0.){
      coutE(InputArguments) << "RooNCSplinePdf_fast::interpolateFcn: Matrix A could not be inverted. Something is wrongwith the x coordinates of points!" << endl;
      assert(0);
    }
    TVectorD Strans = Ainv*Btrans;
    for (int ib=0; ib<npoints-1; ib++){
      vector<Double_t> coef = getCoefficients(Strans, kappas, ib);
      //cout << "Coefficients at bin " << ib << ": " << endl;
      //for (unsigned int ic=0; ic<coef.size(); ic++) cout << coef.at(ic) << " ";
      //cout << endl;
      coefficients.push_back(coef);
    }
  }
  else{
    vector<Double_t> coef;
    coef.push_back((dynamic_cast<RooAbsReal*>(YList.at(0)))->getVal());
    coefficients.push_back(coef);
  }
}

RooNCSplinePdf_fast::RooNCSplinePdf_fast(
  const RooNCSplinePdf_fast& other,
  const char* name
  ) : RooNCSplinePdf(other,name),
  kappas(other.kappas)
{
  for (unsigned int ibin=0; ibin<other.coefficients.size(); ibin++){
    vector<Double_t> coef;
    for (unsigned int ic=0; ic<other.coefficients.at(ibin).size(); ic++) coef.push_back(other.coefficients.at(ibin).at(ic));
    coefficients.push_back(coef);
  }
}

Double_t RooNCSplinePdf_fast::interpolateFcn(Int_t code, const char* rangeName)const{
  const Double_t d_epsilon = 0;
  Double_t res = d_epsilon;

  if (npoints>1){
    if (code==0){
      Int_t bin = getWhichBin(theVar);
      Double_t tv = getTVar(kappas, theVar, bin);

      res=0;
      vector<Double_t> coef = coefficients.at(bin);
      for (int ip=0; ip<4; ip++) res += coef.at(ip)*pow(tv, ip);
    }
    else{
      Double_t varMin = theVar.min(rangeName); Int_t binlow = getWhichBin(varMin); Double_t tvlow = getTVar(kappas, varMin, binlow);
      Double_t varMax = theVar.max(rangeName); Int_t binhigh = getWhichBin(varMax); Double_t tvhigh = getTVar(kappas, varMax, binhigh);

      res=0;
      for (int ib=binlow; ib<=binhigh; ib++){
        Double_t tmax=1., tmin=0.;
        if (ib==binlow) tmin=tvlow;
        else if (ib==binhigh) tmax=tvhigh;
        vector<Double_t> coef = coefficients.at(ib);
        for (int ip=1; ip<=4; ip++) res += coef.at(ip-1)*(pow(tmax, ip)-pow(tmin, ip))/((Double_t)ip)/kappas.at(ib);
      }
    }
  }
  else if(npoints==1){
    res = coefficients.at(0).at(0);
    if (code!=0 && code%2==0) res *= (theVar.max(rangeName)-theVar.min(rangeName));
  }

  return res;
}

Double_t RooNCSplinePdf_fast::evaluate() const{
  Double_t value = interpolateFcn(0);
  if (value<=0.) return 1e-15;
  return value;
}

Double_t RooNCSplinePdf_fast::analyticalIntegral(Int_t code, const char* rangeName) const{
  Double_t value = interpolateFcn(code, rangeName);
  if (value<=0.) value = 1e-10;
  return value;
}

ClassImp(RooNCSplinePdf_fast)
