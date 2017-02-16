#include "RooNCSplinePdfCore.h" 
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


RooNCSplinePdfCore::RooNCSplinePdfCore() :
RooAbsPdf()
{}

RooNCSplinePdfCore::RooNCSplinePdfCore(
  const char* name,
  const char* title
  ) :
  RooAbsPdf(name, title),
  verbosity(RooNCSplinePdfCore::kSilent)
{}

RooNCSplinePdfCore::RooNCSplinePdfCore(
  const RooNCSplinePdfCore& other,
  const char* name
  ) :
  RooAbsPdf(other, name),
  verbosity(other.verbosity)
{}


void RooNCSplinePdfCore::getBArray(const std::vector<Double_t>& kappas, const vector<Double_t>& fcnList, std::vector<Double_t>& BArray)const{
  BArray.clear();
  int npoints=kappas.size();
  if (npoints!=(int)fcnList.size()){
    coutE(InputArguments) << "RooNCSplinePdfCore::getBArray: Dim(kappas)=" << npoints << " != Dim(fcnList)=" << fcnList.size() << endl;
    assert(0);
  }
  BArray.push_back(3.*(fcnList.at(1)-fcnList.at(0)));
  for (int j=1; j<npoints-1; j++){
    Double_t val_j = fcnList.at(j);
    Double_t val_jpo = fcnList.at(j+1);
    Double_t val_jmo = fcnList.at(j-1);
    Double_t kappa_j = kappas.at(j);
    Double_t kappa_jmo = kappas.at(j-1);
    Double_t rsq = pow(kappa_j/kappa_jmo, 2);
    Double_t Bval = val_jpo*rsq + val_j*(1.-rsq) - val_jmo;
    Bval *= 3.;
    BArray.push_back(Bval);
  }
  BArray.push_back(3.*(fcnList.at(npoints-1)-fcnList.at(npoints-2)));
}
void RooNCSplinePdfCore::getAArray(const vector<Double_t>& kappas, vector<vector<Double_t>>& AArray)const{
  AArray.clear();
  Int_t npoints = kappas.size();
  for (int i=0; i<npoints; i++){
    vector<Double_t> Ai;
    Double_t* Aiarray = new Double_t[npoints];
    for (int j=0; j<npoints; j++) Aiarray[j]=0;

    if (i==0){ Aiarray[0]=2; Aiarray[1]=kappas.at(1)/kappas.at(0); }
    else if (i==npoints-1){ Aiarray[npoints-2]=1; Aiarray[npoints-1]=2.*kappas.at(npoints-1)/kappas.at(npoints-2); }
    else{
      Double_t kappa_j = kappas.at(i);
      Double_t kappa_jmo = kappas.at(i-1);
      Double_t kappa_jpo = kappas.at(i+1);

      Aiarray[i-1]=1;
      Aiarray[i]=2.*kappa_j/kappa_jmo*(1.+kappa_j/kappa_jmo);
      Aiarray[i+1]=kappa_j*kappa_jpo/pow(kappa_jmo, 2);
    }
    for (int j=0; j<npoints; j++) Ai.push_back(Aiarray[j]);
    delete[] Aiarray;
    AArray.push_back(Ai);
  }
}

vector<Double_t> RooNCSplinePdfCore::getCoefficients(const TVectorD& S, const vector<Double_t>& kappas, const vector<Double_t>& fcnList, const Int_t& bin)const{
  Double_t A, B, C, D;
  vector<Double_t> res;

  A=fcnList.at(bin);
  B=S[bin];

  Double_t dFcn=fcnList.at(bin+1)-A;
  C=3.*dFcn - 2.*B - S[bin+1]*kappas.at(bin+1)/kappas.at(bin);
  D=-2.*dFcn + B + S[bin+1]*kappas.at(bin+1)/kappas.at(bin);

  res.push_back(A);
  res.push_back(B);
  res.push_back(C);
  res.push_back(D);
  return res;
}
vector<vector<Double_t>> RooNCSplinePdfCore::getCoefficientsAlongDirection(const std::vector<Double_t>& kappas, const TMatrixD& Ainv, const vector<Double_t>& fcnList, const Int_t pickBin)const{
  vector<Double_t> BArray;
  getBArray(kappas, fcnList, BArray);

  Int_t npoints = BArray.size();
  TVectorD Btrans(npoints);
  for (int i=0; i<npoints; i++) Btrans[i]=BArray.at(i);
  TVectorD Strans = Ainv*Btrans;

  vector<vector<Double_t>> coefs;
  for (int bin=0; bin<npoints-1; bin++){
    if (pickBin>=0 && bin!=pickBin) continue;
    vector<Double_t> coef = getCoefficients(Strans, kappas, fcnList, bin);
    coefs.push_back(coef);
  }
  return coefs;
}

Double_t RooNCSplinePdfCore::evalSplineSegment(const std::vector<Double_t>& coefs, const Double_t kappa, Double_t tup, Double_t tdn, Bool_t doIntegrate)const{
  Double_t res=0;
  for (unsigned int ic=0; ic<coefs.size(); ic++){
    if (doIntegrate) res += coefs.at(ic)*(pow(tup, (int)(ic+1))-pow(tdn, (int)(ic+1)))/((Double_t)(ic+1))/kappa;
    else res += coefs.at(ic)*pow(tup, (int)ic);
  }
  return res;
}

void RooNCSplinePdfCore::setVerbosity(VerbosityLevel flag){ verbosity=flag; }


ClassImp(RooNCSplinePdfCore)
