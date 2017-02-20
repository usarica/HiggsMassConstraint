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
RooNCSplinePdf_1D(),
YList("YList", "List of Y coordinates", this),
npointsY(0)
{}

RooNCSplinePdf_2D::RooNCSplinePdf_2D(
  const char* name,
  const char* title,
  RooAbsReal& inXVar,
  RooAbsReal& inYVar,
  const RooArgList& inXList, // X and Y define the grid
  const RooArgList& inYList,
  const RooArgList& inFcnList // Z has dimension dim(X)*dim(Y) with Z[i][j] corresponding to X[i], Y[j]
  ) :
  RooNCSplinePdf_1D(name, title, inXVar, inXList, inFcnList),
  theYVar("theYVar", "theYVar", this, inYVar),
  YList("YList", "List of Y coordinates", this)
{
  TIterator* coefIter = inYList.createIterator();
  RooAbsArg* coef;
  while ((coef = (RooAbsArg*)coefIter->Next())){
    if (!dynamic_cast<RooAbsReal*>(coef)){
      coutE(InputArguments) << "RooNCSplinePdf_2D ERROR::RooNCSplinePdf_2D(" << GetName() << ") Y variable " << coef->GetName() << " is not of type RooAbsReal" << endl;
      assert(0);
    }
    YList.add(*coef);
  }
  delete coefIter;

  npointsY = YList.getSize();

  if (npointsY<=1){
    coutE(InputArguments) << "RooNCSplinePdf_2D ERROR::RooNCSplinePdf_2D(" << GetName() << ") npointsY " << npointsY << "<=1 cannot have a spline." << endl;
    assert(0);
  }
}

RooNCSplinePdf_2D::RooNCSplinePdf_2D(
  const RooNCSplinePdf_2D& other,
  const char* name
  ) :
  RooNCSplinePdf_1D(other, name),
  theYVar("theYVar", this, other.theYVar),
  YList("YList", this, other.YList),
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
  npoints=kappaY.size();
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
  if (code==0 || code%3!=0){ // Case to just compute the value at y
    ybin = getWhichBin(theYVar, 0);
    ty = getTVar(kappaY, theYVar, ybin, 1);
  }
  else{ // Case to integrate along y
    ybinmin = getWhichBin(theYVar.min(rangeName), 0);
    tymin = getTVar(kappaY, theYVar.min(rangeName), ybinmin, 1);
    ybinmax = getWhichBin(theYVar.max(rangeName), 0);
    tymax = getTVar(kappaY, theYVar.max(rangeName), ybinmax, 1);
  }

  // Get the grid of coefficients
  vector<vector<vector<Double_t>>> coefsAlongY; // [Ax(y),Bx(y),Cx(y),Dx(y)][xbin][ybin]
  int npoldim=0;
  int nxbins=0;
  for (Int_t j=0; j<npointsY; j++){
    vector<vector<Double_t>> xcoefsAtYj = getCoefficientsPerY(kappaX, xAinv, j, -1); // [ix][Ax,Bx,Cx,Dx] at each y_j
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
    if (nxbins!=(int)xcoefsAtYj.size() || npoldim!=(int)xcoefsAtYj.at(0).size()){
      coutE(InputArguments) << "RooNCSplinePdf_2D::interpolateFcn: nxbins!=(int)xcoefsAtYj.size() || npoldim!=(int)xcoefsAtYj.at(0).size()!" << endl;
      assert(0);
    }
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
    else txhigh=tx;

    // Get the x coefficients interpolated across y
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
        else tyhigh=ty;

        theCoef += evalSplineSegment(yCoefs.at(iy), kappaY.at(iy), tyhigh, tylow, (code>0 && code%3==0));
      }
      xCoefs.push_back(theCoef);
    }

    // Evaluate value of spline at x with coefficients evaluated at y
    res += evalSplineSegment(xCoefs, kappaX.at(ix), txhigh, txlow, (code>0 && code%2==0));
  }

  return res;
}

void RooNCSplinePdf_2D::getKappa(vector<Double_t>& kappas, const Int_t whichDirection)const{
  kappas.clear();
  Double_t kappa=1;

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
  for (Int_t j=0; j<npoints-1; j++){
    Double_t val_j = (dynamic_cast<RooAbsReal*>(coord->at(j)))->getVal();
    Double_t val_jpo = (dynamic_cast<RooAbsReal*>(coord->at(j+1)))->getVal();
    kappa = 1./(val_jpo-val_j);
    kappas.push_back(kappa);
  }
  kappas.push_back(kappa); // Push the same kappa_(N-1)=kappa_(N-2) at the end point
}
Int_t RooNCSplinePdf_2D::getWhichBin(const Double_t& val, const Int_t whichDirection)const{
  Int_t bin=-1;
  Double_t valj, valjpo;
  Int_t npoints;
  RooListProxy const* coord;
  if (whichDirection==0){
    coord=&XList;
    npoints=npointsX;
  }
  else{
    coord=&YList;
    npoints=npointsY;
  }

  if (npoints<=1) bin=0;
  else{
    valjpo = (dynamic_cast<RooAbsReal*>(coord->at(0)))->getVal();
    for (Int_t j=0; j<npoints-1; j++){
      valj = (dynamic_cast<RooAbsReal*>(coord->at(j)))->getVal();
      valjpo = (dynamic_cast<RooAbsReal*>(coord->at(j+1)))->getVal();
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
  RooListProxy const* coord;
  if (whichDirection==0) coord=&XList;
  else coord=&YList;
  return (Double_t)((val-(dynamic_cast<RooAbsReal*>(coord->at(bin)))->getVal())*K);
}

vector<vector<Double_t>> RooNCSplinePdf_2D::getCoefficientsPerY(const std::vector<Double_t>& kappaX, const TMatrixD& xAinv, const Int_t& ybin, const Int_t xbin)const{
  vector<Double_t> fcnList;
  for (int bin=0; bin<npointsX; bin++){ fcnList.push_back((dynamic_cast<RooAbsReal*>(FcnList.at(npointsY*bin+ybin)))->getVal()); }
  vector<vector<Double_t>> coefs = getCoefficientsAlongDirection(kappaX, xAinv, fcnList, xbin);
  //for (unsigned int bin=0; bin<fcnList.size(); bin++) cout << "RooNCSplinePdf_2D::getCoefficientsPerY: fcnList[" << bin << " | ybin = " << ybin << ", xbin = " << xbin << "] = " << fcnList[bin] << endl;
  //for (unsigned int bin=0; bin<coefs.size(); bin++){
  //  cout << "RooNCSplinePdf_2D::getCoefficientsPerY: coefs[" << bin << "] = ";
  //  for (unsigned int ic=0; ic<coefs.at(bin).size(); ic++){
  //    if (ic<coefs.at(bin).size()-1) cout << coefs.at(bin).at(ic) << " ";
  //    else cout << coefs.at(bin).at(ic) << endl;
  //  }
  //}
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
