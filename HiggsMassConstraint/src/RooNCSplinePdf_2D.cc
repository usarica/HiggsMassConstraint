#include "RooNCSplinePdf_2D.h" 
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


RooNCSplinePdf_2D::RooNCSplinePdf_2D() :
RooAbsPdf(),
theXVar("theXVar", "theXVar", this),
XList("XList", "List of X coordinates", this),
YList("YList", "List of Y coordinates", this),
npointsX(0),
npointsY(0)
{}

RooNCSplinePdf_2D::RooNCSplinePdf_2D(
  const char* name,
  const char* title,
  RooAbsReal& inXVar,
  RooAbsReal& inYVar,
  const RooArgList& inXList, // X and Y define the grid
  const RooArgList& inYList,
  const RooArgList& inZList // Z has dimension dim(X)*dim(Y) with Z[i][j] corresponding to X[i], Y[j]
  ) :
  RooAbsPdf(name, title),
  theXVar("theXVar", "theXVar", this, inXVar),
  theYVar("theYVar", "theYVar", this, inYVar),
  XList("XList", "List of X coordinates", this),
  YList("YList", "List of Y coordinates", this),
  ZList("ZList", "List of Z coordinates", this)
{
  TIterator* coefIter = inXList.createIterator();
  RooAbsArg* coef;
  while ((coef = (RooAbsArg*)coefIter->Next())){
    if (!dynamic_cast<RooAbsReal*>(coef)){
      coutE(InputArguments) << "RooNCSplinePdf_2D ERROR::RooNCSplinePdf_2D(" << GetName() << ") X variable " << coef->GetName() << " is not of type RooAbsReal" << endl;
      assert(0);
    }
    XList.add(*coef);
  }
  delete coefIter;

  coefIter = inYList.createIterator();
  while ((coef = (RooAbsArg*)coefIter->Next())){
    if (!dynamic_cast<RooAbsReal*>(coef)){
      coutE(InputArguments) << "RooNCSplinePdf_2D ERROR::RooNCSplinePdf_2D(" << GetName() << ") Y variable " << coef->GetName() << " is not of type RooAbsReal" << endl;
      assert(0);
    }
    YList.add(*coef);
  }
  delete coefIter;

  coefIter = inZList.createIterator();
  while ((coef = (RooAbsArg*)coefIter->Next())){
    if (!dynamic_cast<RooAbsReal*>(coef)){
      coutE(InputArguments) << "RooNCSplinePdf_2D ERROR::RooNCSplinePdf_2D(" << GetName() << ") Y variable " << coef->GetName() << " is not of type RooAbsReal" << endl;
      assert(0);
    }
    ZList.add(*coef);
  }
  delete coefIter;

  npointsX = XList.getSize();
  npointsY = YList.getSize();
  Int_t nZ = ZList.getSize();
  if (nZ!=npointsX*npointsY){
    coutE(InputArguments) << "RooNCSplinePdf_2D ERROR::RooNCSplinePdf_2D(" << GetName() << ") nX " << npointsX << " * nY " << npointsY << " != nZ " << nZ << endl;
    assert(0);
  }

  if (npointsX<=1 || npointsY<=1){
    coutE(InputArguments) << "RooNCSplinePdf_2D ERROR::RooNCSplinePdf_2D(" << GetName() << ") npointsX " << npointsX << " or npointsY " << npointsY << "<=1 cannot have a spline." << endl;
    assert(0);
  }
}

RooNCSplinePdf_2D::RooNCSplinePdf_2D(
  const RooNCSplinePdf_2D& other,
  const char* name
  ) :
  RooAbsPdf(other, name),
  theXVar("theXVar", this, other.theXVar),
  theYVar("theYVar", this, other.theYVar),
  XList("XList", this, other.XList),
  YList("YList", this, other.YList),
  ZList("ZList", this, other.ZList),
  npointsX(other.npointsX),
  npointsY(other.npointsY)
{}


Double_t RooNCSplinePdf_2D::interpolateFcn(Int_t code, const char* rangeName)const{
  Double_t res=0;

  // Prepare A and kappa arrays for x and y coordinates
  int npoints;
  Double_t det;

  vector<vector<Double_t>> xA; vector<Double_t> kappaX; getKappa(kappaX, 0); getAArray(kappaX, xA);
  npoints=kappaX.size();
  TMatrixD xAtrans(npoints, npoints);
  for (int i=0; i<npoints; i++){ for (int j=0; j<npoints; j++){ xAtrans[i][j]=xA.at(i).at(j); } }
  det=0;
  TMatrixD xAinv = xAtrans.Invert(&det);
  if (det==0.){
    coutE(InputArguments) << "RooNCSplinePdf_2D::interpolateFcn: Matrix xA could not be inverted. Something is wrong with the x coordinates of points!" << endl;
    assert(0);
  }

  vector<vector<Double_t>> yA; vector<Double_t> kappaY; getKappa(kappaY, 1); getAArray(kappaY, yA);
  TMatrixD yAtrans(npoints, npoints);
  for (int i=0; i<npoints; i++){ for (int j=0; j<npoints; j++){ yAtrans[i][j]=yA.at(i).at(j); } }
  det=0;
  TMatrixD yAinv = yAtrans.Invert(&det);
  if (det==0.){
    coutE(InputArguments) << "RooNCSplinePdf_2D::interpolateFcn: Matrix yA could not be inverted. Something is wrong with the y coordinates of points!" << endl;
    assert(0);
  }

  // Get bins
  Int_t xbin=-1, xbinmin=-1, xbinmax=-1, ybin=-1, ybinmin=-1, ybinmax=-1;
  Double_t tx=0, txmin=0, txmax=0, ty=0, tymin=0, tymax=0;
  if (code==0 || code%2!=0){ // Just compute the value
    xbin = getWhichBin(theXVar, 0);
    tx = getTVar(kappaX, theXVar, xbin, 0);
  }
  else{ // Integrate along x
    xbinmin = getWhichBin(theXVar.min(rangeName), 0);
    txmin = getTVar(kappaX, theXVar.min(rangeName), xbinmin, 0);
    xbinmax = getWhichBin(theXVar.max(rangeName), 0);
    txmax = getTVar(kappaX, theXVar.max(rangeName), xbinmax, 0);
  }
  if (code==0 || code%3!=0){ // Just compute the value
    ybin = getWhichBin(theYVar, 0);
    ty = getTVar(kappaY, theYVar, ybin, 1);
  }
  else{ // Integrate along y
    ybinmin = getWhichBin(theYVar.min(rangeName), 0);
    tymin = getTVar(kappaY, theYVar.min(rangeName), ybinmin, 1);
    ybinmax = getWhichBin(theYVar.max(rangeName), 0);
    tymax = getTVar(kappaY, theYVar.max(rangeName), ybinmax, 1);
  }

  // Get the grid of coefficients
  vector<vector<vector<Double_t>>> coefsAlongY; // [A,B,C,D][xbin][ybin]
  int npoldim=0;
  int nxbins=0;
  for (Int_t j=0; j<npointsY; j++){
    vector<vector<Double_t>> xcoefsAtYj = getCoefficientsPerY(kappaX, xAinv, j, -1); // [ix][A,B,C,D] at each y_j
    if (j==0){
      nxbins=xcoefsAtYj.size();
      npoldim=xcoefsAtYj.at(0).size();
      for (int ipow=0; ipow<npoldim; ipow++){
        vector<vector<Double_t>> dum_xycoefarray;
        for (int ix=0; ix<nxbins; ix++){
          vector<Double_t> dum_ycoefarray;
          dum_xycoefarray.push_back(dum_ycoefarray);
        }
        coefsAlongY.push_back(dum_xycoefarray);
      }
    }
    if (nxbins!=(int)xcoefsAtYj.size() || npoldim!=(int)xcoefsAtYj.at(0).size()) assert(0);
    for (int ix=0; ix<nxbins; ix++){
      for (int ipow=0; ipow<npoldim; ipow++) coefsAlongY.at(ipow).at(ix).push_back(xcoefsAtYj.at(ix).at(ipow));
    }
  }

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

    vector<Double_t> xCoefs;
    for (int ic=0; ic<npoldim; ic++){
      vector<vector<Double_t>> yCoefs = getCoefficientsAlongDirection(kappaY, yAinv, coefsAlongY.at(ic).at(ix), -1); // [iy][A,B,C,D]

      Double_t theCoef=0;
      for (int iy=0; iy<(int)yCoefs.size(); iy++){
        if (
          (ybin>=0 && iy!=ybin)
          ||
          (ybinmin>=0 && ybinmax>=ybinmin && !(ybinmin<=iy && iy<=ybinmax))
          ) continue;
        Double_t tylow=0, tyhigh=1;
        if (code>0 && code%3==0){
          if (iy==ybinmin) tylow=tymin;
          if (iy==ybinmax) tyhigh=tymax;
        }
        for (unsigned int icy=0; icy<yCoefs.at(iy).size(); icy++){
          if (code>0 && code%3==0) theCoef += yCoefs.at(iy).at(icy)*(pow(tyhigh, (int)(icy+1))-pow(tylow, (int)(icy+1)))/((Double_t)(icy+1));
          else theCoef += yCoefs.at(iy).at(icy)*pow(ty, (int)icy);
        }
      }
      xCoefs.push_back(theCoef);
    }
    for (unsigned int icx=0; icx<xCoefs.size(); icx++){
      if (code>0 && code%3==0) res += xCoefs.at(icx)*(pow(txhigh, (int)(icx+1))-pow(txlow, (int)(icx+1)))/((Double_t)(icx+1));
      else res += xCoefs.at(icx)*pow(tx, (int)icx);
    }
  }

  return res;
}

void RooNCSplinePdf_2D::getKappa(vector<Double_t>& kappas, const Int_t whichDirection)const{
  kappas.clear();
  Double_t kappa;

  Int_t npoints;
  RooListProxy const* coord;
  if (whichDirection==0){
    npoints=npointsX;
    coord=&XList;
  }
  else{
    npoints=npointsY;
    coord=&YList;
  }
  for (int j=0; j<npoints-1; j++){
    Double_t val_j = (dynamic_cast<RooAbsReal*>(coord->at(j)))->getVal();
    Double_t val_jpo = (dynamic_cast<RooAbsReal*>(coord->at(j+1)))->getVal();
    kappa = 1./(val_jpo-val_j);
    kappas.push_back(kappa);
  }
  kappas.push_back(kappa); // Push the same kappa_(N-1)=kappa_(N-2) at the end point
}
void RooNCSplinePdf_2D::getBArray(const std::vector<Double_t>& kappas, const vector<Double_t>& fcnList, std::vector<Double_t>& BArray)const{
  BArray.clear();
  int npoints=kappas.size();
  if (npoints!=(int)fcnList.size()){
    coutE(InputArguments) << "RooNCSplinePdf_2D::getBArray: Dim(kappas)=" << npoints << " != Dim(fcnList)=" << fcnList.size() << endl;
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
void RooNCSplinePdf_2D::getAArray(const vector<Double_t>& kappas, vector<vector<Double_t>>& AArray)const{
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
Int_t RooNCSplinePdf_2D::getWhichBin(const Double_t& val, const Int_t whichDirection)const{
  Int_t bin=-1;
  Double_t valj, valjpo;
  int npoints;
  RooListProxy const* gridDir;
  if (whichDirection==0){
    gridDir=&XList;
    npoints=npointsX;
  }
  else{
    gridDir=&YList;
    npoints=npointsY;
  }

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
Double_t RooNCSplinePdf_2D::getTVar(const vector<Double_t>& kappas, const Double_t& val, const Int_t& bin, const Int_t whichDirection)const{
  Double_t K;
  K=kappas.at(bin);
  if (whichDirection==0) return (Double_t)((val-(dynamic_cast<RooAbsReal*>(XList.at(bin)))->getVal())*K);
  else return (Double_t)((val-(dynamic_cast<RooAbsReal*>(YList.at(bin)))->getVal())*K);
}

vector<Double_t> RooNCSplinePdf_2D::getCoefficients(const TVectorD& S, const vector<Double_t>& kappas, const vector<Double_t>& fcnList, const Int_t& bin)const{
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
vector<vector<Double_t>> RooNCSplinePdf_2D::getCoefficientsPerY(const std::vector<Double_t>& kappaX, const TMatrixD& xAinv, const Int_t& ybin, const Int_t xbin)const{
  vector<Double_t> fcnList;
  for (int bin=0; bin<npointsX; bin++) fcnList.push_back((dynamic_cast<RooAbsReal*>(ZList.at(npointsX*ybin+bin)))->getVal());
  vector<vector<Double_t>> coefs = getCoefficientsAlongDirection(kappaX, xAinv, fcnList, xbin);
  return coefs;
}
vector<vector<Double_t>> RooNCSplinePdf_2D::getCoefficientsAlongDirection(const std::vector<Double_t>& kappas, const TMatrixD& Ainv, const vector<Double_t>& fcnList, const Int_t pickBin)const{
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

Double_t RooNCSplinePdf_2D::evaluate() const{
  Double_t value = interpolateFcn(0);
  if (value<=0.) return 1e-15;
  return value;
}
Int_t RooNCSplinePdf_2D::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const{
  Int_t code=1;
  if (dynamic_cast<RooRealVar*>(theXVar.absArg())!=0){
    if (matchArgs(allVars, analVars, theXVar)) code*=2;
  }
  if (dynamic_cast<RooRealVar*>(theYVar.absArg())!=0){
    if (matchArgs(allVars, analVars, theYVar)) code*=3;
  }
  if (code==1) code=0;
  return code;
}
Double_t RooNCSplinePdf_2D::analyticalIntegral(Int_t code, const char* rangeName) const{
  Double_t value = interpolateFcn(code, rangeName);
  if (value<=0.) value = 1e-10;
  return value;
}

ClassImp(RooNCSplinePdf_2D)
