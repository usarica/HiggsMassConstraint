#include "RooNCSplinePdf_2D_fast.h" 
#include <cmath>
#include "TMath.h"
#include "Riostream.h" 
#include "RooAbsReal.h" 

using namespace TMath;
using namespace RooFit;
using namespace std;


RooNCSplinePdf_2D_fast::RooNCSplinePdf_2D_fast() :
RooNCSplinePdfCore(),
theYVar("theYVar", "theYVar", this)
{}

RooNCSplinePdf_2D_fast::RooNCSplinePdf_2D_fast(
  const char* name,
  const char* title,
  RooAbsReal& inXVar,
  RooAbsReal& inYVar,
  std::vector<T>& inXList,
  std::vector<T>& inYList,
  std::vector<std::vector<T>>& inFcnList
  ) :
  RooNCSplinePdfCore(name, title, inXVar, inXList),
  theYVar("theYVar", "theYVar", this, inYVar),
  YList(inYList),
  FcnList(inFcnList)
{
  if (npointsX()>1 && npointsY()>1){
    // Prepare A and kappa arrays for x and y coordinates
    int npoints;
    Double_t det;

    vector<vector<RooNCSplinePdfCore::T>> xA; getKappas(kappaX, 0); getAArray(kappaX, xA);
    npoints=kappaX.size();
    TMatrix_t xAtrans(npoints, npoints);
    for (int i=0; i<npoints; i++){ for (int j=0; j<npoints; j++){ xAtrans[i][j]=xA.at(i).at(j); } }
    det=0;
    TMatrix_t xAinv = xAtrans.Invert(&det);
    if (det==0.){
      coutE(InputArguments) << "RooNCSplinePdf_2D_fast::interpolateFcn: Matrix xA could not be inverted. Something is wrong with the x coordinates of points!" << endl;
      assert(0);
    }

    vector<vector<RooNCSplinePdfCore::T>> yA; getKappas(kappaY, 1); getAArray(kappaY, yA);
    npoints=kappaY.size();
    TMatrix_t yAtrans(npoints, npoints);
    for (int i=0; i<npoints; i++){ for (int j=0; j<npoints; j++){ yAtrans[i][j]=yA.at(i).at(j); } }
    det=0;
    TMatrix_t yAinv = yAtrans.Invert(&det);
    if (det==0.){
      coutE(InputArguments) << "RooNCSplinePdf_2D_fast::interpolateFcn: Matrix yA could not be inverted. Something is wrong with the y coordinates of points!" << endl;
      assert(0);
    }

    // Get the grid of coefficients
    vector<vector<vector<RooNCSplinePdfCore::T>>> coefsAlongY; // [Ax(y),Bx(y),Cx(y),Dx(y)][xbin][ybin]
    int npoldim=0;
    int nxbins=0;
    for (unsigned int j=0; j<npointsY(); j++){
      vector<vector<RooNCSplinePdfCore::T>> xcoefsAtYj = getCoefficientsPerY(kappaX, xAinv, j, -1); // [ix][Ax,Bx,Cx,Dx] at each y_j
      if (j==0){
        nxbins=xcoefsAtYj.size();
        npoldim=xcoefsAtYj.at(0).size();
        for (int ipow=0; ipow<npoldim; ipow++){
          vector<vector<RooNCSplinePdfCore::T>> dum_xycoefarray;
          for (int ix=0; ix<nxbins; ix++){
            vector<RooNCSplinePdfCore::T> dum_ycoefarray;
            dum_xycoefarray.push_back(dum_ycoefarray);
          }
          coefsAlongY.push_back(dum_xycoefarray);
        }
      }
      if (nxbins!=(int)xcoefsAtYj.size() || npoldim!=(int)xcoefsAtYj.at(0).size()){
        coutE(InputArguments) << "RooNCSplinePdf_2D_fast::interpolateFcn: nxbins!=(int)xcoefsAtYj.size() || npoldim!=(int)xcoefsAtYj.at(0).size()!" << endl;
        assert(0);
      }
      for (int ix=0; ix<nxbins; ix++){
        for (int ipow=0; ipow<npoldim; ipow++) coefsAlongY.at(ipow).at(ix).push_back(xcoefsAtYj.at(ix).at(ipow));
      }
    }

    for (int ix=0; ix<nxbins; ix++){
      // Get the x coefficients interpolated across y
      vector<vector<vector<RooNCSplinePdfCore::T>>> xCoefs;
      for (int ic=0; ic<npoldim; ic++){
        vector<vector<RooNCSplinePdfCore::T>> yCoefs = getCoefficientsAlongDirection(kappaY, yAinv, coefsAlongY.at(ic).at(ix), -1); // [iy][A,B,C,D]
        xCoefs.push_back(yCoefs);
      }
      coefficients.push_back(xCoefs);
    }
  }
  else assert(0);

  emptyFcnList();
}

RooNCSplinePdf_2D_fast::RooNCSplinePdf_2D_fast(
  const RooNCSplinePdf_2D_fast& other,
  const char* name
  ) :
  RooNCSplinePdfCore(other, name),
  theYVar("theYVar", this, other.theYVar),
  YList(other.YList),
  FcnList(other.FcnList),
  kappaX(other.kappaX),
  kappaY(other.kappaY),
  coefficients(other.coefficients)
{}


RooNCSplinePdfCore::T RooNCSplinePdf_2D_fast::interpolateFcn(Int_t code, const char* rangeName)const{
  DefaultAccumulator<RooNCSplinePdfCore::T> res=RooNCSplinePdfCore::T(0);

  if (verbosity==RooNCSplinePdfCore::kVerbose){ cout << "RooNCSplinePdf_2D_fast(" << GetName() << ")::interpolateFcn begin with code: " << code << endl; }

  // Get bins
  Int_t xbin=-1, xbinmin=-1, xbinmax=-1, ybin=-1, ybinmin=-1, ybinmax=-1;
  RooNCSplinePdfCore::T tx=0, txmin=0, txmax=0, ty=0, tymin=0, tymax=0;
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
    ybin = getWhichBin(theYVar, 1);
    ty = getTVar(kappaY, theYVar, ybin, 1);
  }
  else{ // Case to integrate along y
    ybinmin = getWhichBin(theYVar.min(rangeName), 1);
    tymin = getTVar(kappaY, theYVar.min(rangeName), ybinmin, 1);
    ybinmax = getWhichBin(theYVar.max(rangeName), 1);
    tymax = getTVar(kappaY, theYVar.max(rangeName), ybinmax, 1);
  }

  for (int ix=0; ix<(int)coefficients.size(); ix++){
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

    if (verbosity==RooNCSplinePdfCore::kVerbose){
      if (code==0 || code%2!=0) cout << "Evaluating tx=" << txhigh << " in bin " << ix << endl;
      else cout << "Evaluating tx[" << txlow << ", " << txhigh << "] in bin " << ix << endl;
    }

    // Get the x coefficients interpolated across y
    vector<RooNCSplinePdfCore::T> xCoefs;
    for (int ic=0; ic<(int)coefficients.at(ix).size(); ic++){
      const vector<vector<RooNCSplinePdfCore::T>>& yCoefs = coefficients.at(ix).at(ic);

      if (verbosity==RooNCSplinePdfCore::kVerbose) cout << "\tCoefficient " << ic << ":\n";

      DefaultAccumulator<RooNCSplinePdfCore::T> theCoef=RooNCSplinePdfCore::T(0);
      for (int iy=0; iy<(int)yCoefs.size(); iy++){
        if (
          (ybin>=0 && iy!=ybin)
          ||
          (ybinmin>=0 && ybinmax>=ybinmin && !(ybinmin<=iy && iy<=ybinmax))
          ) continue;

        RooNCSplinePdfCore::T tylow=0, tyhigh=1;
        if (code>0 && code%3==0){
          if (iy==ybinmin) tylow=tymin;
          if (iy==ybinmax) tyhigh=tymax;
        }
        else tyhigh=ty;

        if (verbosity==RooNCSplinePdfCore::kVerbose){
          if (code==0 || code%3!=0) cout << "\tEvaluating ty=" << tyhigh << " in bin " << iy << endl;
          else cout << "\tEvaluating ty[" << tylow << ", " << tyhigh << "] in bin " << iy << endl;
        }

        theCoef += evalSplineSegment(yCoefs.at(iy), kappaY.at(iy), tyhigh, tylow, (code>0 && code%3==0));
      }

      //if (code==0) cout << "\tCoefficient is " << theCoef << endl;

      xCoefs.push_back(theCoef.sum());
    }

    // Evaluate value of spline at x with coefficients evaluated at y
    res += evalSplineSegment(xCoefs, kappaX.at(ix), txhigh, txlow, (code>0 && code%2==0));
  }

  return res.sum();
}


void RooNCSplinePdf_2D_fast::getKappas(vector<RooNCSplinePdfCore::T>& kappas, const Int_t whichDirection)const{
  kappas.clear();
  RooNCSplinePdfCore::T kappa=1;

  Int_t npoints;
  vector<RooNCSplinePdfCore::T> const* coord;
  if (whichDirection==0){
    npoints=npointsX();
    coord=&XList;
  }
  else{
    npoints=npointsY();
    coord=&YList;
  }

  for (Int_t j=0; j<npoints-1; j++){
    RooNCSplinePdfCore::T val_j = coord->at(j);
    RooNCSplinePdfCore::T val_jpo = coord->at(j+1);
    kappa = 1./(val_jpo-val_j);
    kappas.push_back(kappa);
  }
  kappas.push_back(kappa); // Push the same kappa_(N-1)=kappa_(N-2) at the end point
}
Int_t RooNCSplinePdf_2D_fast::getWhichBin(const RooNCSplinePdfCore::T& val, const Int_t whichDirection)const{
  Int_t bin=-1;
  RooNCSplinePdfCore::T valj, valjpo;
  Int_t npoints;
  vector<RooNCSplinePdfCore::T> const* coord;
  if (whichDirection==0){
    coord=&XList;
    npoints=npointsX();
  }
  else{
    coord=&YList;
    npoints=npointsY();
  }

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
RooNCSplinePdfCore::T RooNCSplinePdf_2D_fast::getTVar(const vector<RooNCSplinePdfCore::T>& kappas, const RooNCSplinePdfCore::T& val, const Int_t& bin, const Int_t whichDirection)const{
  RooNCSplinePdfCore::T K;
  K=kappas.at(bin);
  vector<RooNCSplinePdfCore::T> const* coord;
  if (whichDirection==0) coord=&XList;
  else coord=&YList;
  return (val-coord->at(bin))*K;
}

vector<vector<RooNCSplinePdfCore::T>> RooNCSplinePdf_2D_fast::getCoefficientsPerY(const std::vector<RooNCSplinePdfCore::T>& kappaX, const TMatrix_t& xAinv, const Int_t& ybin, const Int_t xbin)const{
  vector<RooNCSplinePdfCore::T> fcnList;
  for (unsigned int bin=0; bin<npointsX(); bin++){ fcnList.push_back(FcnList.at(ybin).at(bin)); }
  vector<vector<RooNCSplinePdfCore::T>> coefs = getCoefficientsAlongDirection(kappaX, xAinv, fcnList, xbin);
  return coefs;
}

Double_t RooNCSplinePdf_2D_fast::evaluate() const{
  Double_t value = interpolateFcn(0);
  if (value<=0.){
    coutE(Eval) << "RooNCSplinePdf_2D_fast ERROR::RooNCSplinePdf_2D_fast(" << GetName() << ") evaluation returned " << value << " at (x, y) = (" << theXVar << ", " << theYVar << ")" << endl;
    value = 1e-15;
  }
  if (verbosity==RooNCSplinePdfCore::kVerbose){ cout << "RooNCSplinePdf_2D_fast(" << GetName() << ")::evaluate = " << value << " at (x, y) = (" << theXVar << ", " << theYVar << ")" << endl; }
  return value;
}
Int_t RooNCSplinePdf_2D_fast::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const{
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
Double_t RooNCSplinePdf_2D_fast::analyticalIntegral(Int_t code, const char* rangeName) const{
  Double_t value = interpolateFcn(code, rangeName);
  if (value<=0.){
    coutE(Integration) << "RooNCSplinePdf_2D_fast ERROR::RooNCSplinePdf_2D_fast(" << GetName() << ") integration returned " << value << " for code = " << code << endl;
    value = 1e-10;
  }
  if (verbosity==RooNCSplinePdfCore::kVerbose){ cout << "RooNCSplinePdf_2D_fast(" << GetName() << ")::analyticalIntegral = " << value << " for code = " << code << endl; }
  return value;
}

ClassImp(RooNCSplinePdf_2D_fast)
