#include "RooNCSplinePdf.h" 
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


RooNCSplinePdf::RooNCSplinePdf() :
RooAbsPdf(),
theVar("theVar", "theVar", this),
XList("XList", "List of X coordinates", this),
YList("YList", "List of Y coordinates", this),
npoints(0)
{}

RooNCSplinePdf::RooNCSplinePdf(
const char* name,
const char* title,
RooAbsReal& inVar,
const RooArgList& inXList,
const RooArgList& inYList
) :
RooAbsPdf(name, title),
theVar("theVar", "theVar", this, inVar),
XList("XList", "List of X coordinates", this),
YList("YList", "List of Y coordinates", this)
{
  TIterator* coefIter = inXList.createIterator();
  RooAbsArg* coef;
  while ((coef = (RooAbsArg*)coefIter->Next())){
    if (!dynamic_cast<RooAbsReal*>(coef)){
      coutE(InputArguments) << "RooNCSplinePdf ERROR::RooNCSplinePdf(" << GetName() << ") X variable " << coef->GetName() << " is not of type RooAbsReal" << endl;
      assert(0);
    }
    XList.add(*coef);
  }
  delete coefIter;

  coefIter = inYList.createIterator();
  while ((coef = (RooAbsArg*)coefIter->Next())){
    if (!dynamic_cast<RooAbsReal*>(coef)){
      coutE(InputArguments) << "RooNCSplinePdf ERROR::RooNCSplinePdf(" << GetName() << ") Y variable " << coef->GetName() << " is not of type RooAbsReal" << endl;
      assert(0);
    }
    YList.add(*coef);
  }
  delete coefIter;

  Int_t nX = XList.getSize();
  Int_t nY = YList.getSize();
  if (nX!=nY){
    coutE(InputArguments) << "RooNCSplinePdf ERROR::RooNCSplinePdf(" << GetName() << ") nX " << nX << " != nY " << nY << endl;
    assert(0);
  }
  else npoints=nX;

  if (npoints<1){
    coutE(InputArguments) << "RooNCSplinePdf ERROR::RooNCSplinePdf(" << GetName() << ") npoints " << npoints << "<1 cannot have a spline." << endl;
    assert(0);
  }
}

RooNCSplinePdf::RooNCSplinePdf(
  const RooNCSplinePdf& other,
  const char* name
  ) :
  RooAbsPdf(other, name),
  theVar("theVar", this, other.theVar),
  XList("XList", this, other.XList),
  YList("YList", this, other.YList),
  npoints(other.npoints)
{}

Double_t RooNCSplinePdf::interpolateFcn(Int_t code, const char* rangeName)const{
  const Double_t d_epsilon = 0;
  Double_t res = d_epsilon;

  if (npoints>1){
    vector<Double_t> kappas;
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

    //cout << "A array:\n";
    //for (int ix=0; ix<npoints; ix++){
    //  for (int iy=0; iy<npoints; iy++) cout << Atrans[ix][iy] << " ";
    //  cout << endl;
    //}

    Double_t det=0;
    TMatrixD Ainv = Atrans.Invert(&det);
    if (det==0.){
      coutE(InputArguments) << "RooNCSplinePdf::interpolateFcn: Matrix A could not be inverted. Something is wrongwith the x coordinates of points!" << endl;
      assert(0);
    }
    TVectorD Strans = Ainv*Btrans;

    if (code==0){
      Int_t bin = getWhichBin(theVar);
      Double_t tv = getTVar(kappas, theVar, bin);

      res=0;
      vector<Double_t> coef = getCoefficients(Strans, kappas, bin);
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
        vector<Double_t> coef = getCoefficients(Strans, kappas, ib);
        for (int ip=1; ip<=4; ip++) res += coef.at(ip-1)*(pow(tmax, ip)-pow(tmin, ip))/((Double_t)ip)/kappas.at(ib);
      }
    }
  }
  else if(npoints==1){
    res = (dynamic_cast<RooAbsReal*>(YList.at(0)))->getVal();
    if (code!=0 && code%2==0) res *= (theVar.max(rangeName)-theVar.min(rangeName));
  }

  return res;
}

void RooNCSplinePdf::getKappa(vector<Double_t>& kappas)const{
  kappas.clear();
  Double_t kappa;
  for (int j=0; j<npoints-1; j++){
    Double_t val_j = (dynamic_cast<RooAbsReal*>(XList.at(j)))->getVal();
    Double_t val_jpo = (dynamic_cast<RooAbsReal*>(XList.at(j+1)))->getVal();
    kappa = 1./(val_jpo-val_j);
    kappas.push_back(kappa);
  }
  kappas.push_back(kappa); // Push the same kappa_(N-1)=kappa_(N-2) at the end point
}
void RooNCSplinePdf::getBArray(const vector<Double_t>& kappas, vector<Double_t>& BArray)const{
  BArray.clear();
  BArray.push_back(
    3.*((dynamic_cast<RooAbsReal*>(YList.at(1)))->getVal()-(dynamic_cast<RooAbsReal*>(YList.at(0)))->getVal())
    );
  for (int j=1; j<npoints-1; j++){
    Double_t val_j = (dynamic_cast<RooAbsReal*>(YList.at(j)))->getVal();
    Double_t val_jpo = (dynamic_cast<RooAbsReal*>(YList.at(j+1)))->getVal();
    Double_t val_jmo = (dynamic_cast<RooAbsReal*>(YList.at(j-1)))->getVal();
    Double_t kappa_j = kappas.at(j);
    Double_t kappa_jmo = kappas.at(j-1);
    Double_t rsq = pow(kappa_j/kappa_jmo, 2);
    Double_t Bval = val_jpo*rsq + val_j*(1.-rsq) - val_jmo;
    Bval *= 3.;
    BArray.push_back(Bval);
  }
  BArray.push_back(
    3.*((dynamic_cast<RooAbsReal*>(YList.at(npoints-1)))->getVal()-(dynamic_cast<RooAbsReal*>(YList.at(npoints-2)))->getVal())
    );
}
void RooNCSplinePdf::getAArray(const vector<Double_t>& kappas, vector<vector<Double_t>>& AArray)const{
  AArray.clear();
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
Int_t RooNCSplinePdf::getWhichBin(const Double_t& xval)const{
  Int_t bin=-1;
  Double_t xj, xjpo;
  if (npoints<=1) return 0;
  else{
    xjpo = (dynamic_cast<RooAbsReal*>(XList.at(0)))->getVal();
    for (Int_t j=0; j<npoints-1; j++){
      xj = (dynamic_cast<RooAbsReal*>(XList.at(j)))->getVal();
      xjpo = (dynamic_cast<RooAbsReal*>(XList.at(j+1)))->getVal();
      if (xval<xjpo && xval>=xj){ bin=j; break; }
    }
    if (bin==-1 && xval>=xjpo) bin=npoints-2;
    else if (bin==-1) bin=0;
    return bin;
  }
}
Double_t RooNCSplinePdf::getTVar(const vector<Double_t>& kappas, const Double_t& xval, const Int_t& bin)const{
  Double_t K;
  K=kappas.at(bin);
  return (Double_t)((xval-(dynamic_cast<RooAbsReal*>(XList.at(bin)))->getVal())*K);
}
vector<Double_t> RooNCSplinePdf::getCoefficients(const TVectorD& S, const vector<Double_t>& kappas, const Int_t& bin)const{
  Double_t A, B, C, D;
  vector<Double_t> res;

  A=(dynamic_cast<RooAbsReal*>(YList.at(bin)))->getVal();
  B=S[bin];

  Double_t dY=(dynamic_cast<RooAbsReal*>(YList.at(bin+1)))->getVal()-A;
  C=3.*dY - 2.*B - S[bin+1]*kappas.at(bin+1)/kappas.at(bin);
  D=-2.*dY + B + S[bin+1]*kappas.at(bin+1)/kappas.at(bin);

  //cout << "A,B,C,D in bin " << bin << ":\n";
  //cout << "\t" << A << " " << B << " " << C << " " << D << endl;

  res.push_back(A);
  res.push_back(B);
  res.push_back(C);
  res.push_back(D);
  return res;
}

Double_t RooNCSplinePdf::evaluate() const{
  Double_t value = interpolateFcn(0);
  if (value<=0.) return 1e-15;
  return value;
}
Int_t RooNCSplinePdf::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const{
  Int_t code=0;
  if (dynamic_cast<RooRealVar*>(theVar.absArg())!=0){
    if (matchArgs(allVars, analVars, theVar)) code=2;
  }
  return code;
}
Double_t RooNCSplinePdf::analyticalIntegral(Int_t code, const char* rangeName) const{
  Double_t value = interpolateFcn(code, rangeName);
  if (value<=0.) value = 1e-10;
  return value;
}

ClassImp(RooNCSplinePdf)
