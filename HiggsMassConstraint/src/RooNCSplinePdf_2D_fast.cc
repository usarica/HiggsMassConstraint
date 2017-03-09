#include "RooNCSplinePdf_2D_fast.h" 
#include <cmath>
#include "Riostream.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
#include "RooConstVar.h" 
#include "TH3F.h"
#include "TAxis.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "RooDataHist.h"

using namespace TMath;
using namespace RooFit;
using namespace std;


RooNCSplinePdf_2D_fast::RooNCSplinePdf_2D_fast() : RooNCSplinePdf_2D()
{}

RooNCSplinePdf_2D_fast::RooNCSplinePdf_2D_fast(
  const char* name,
  const char* title,
  RooAbsReal* inXVar,
  RooAbsReal* inYVar,
  const RooArgList* inXList,
  const RooArgList* inYList,
  std::vector<const RooArgList*>& inFcnList
  ) : RooNCSplinePdf_2D(name, title, inXVar, inYVar, inXList, inYList, inFcnList, true)
{
  if (npointsX>1 && npointsY>1){
    // Prepare A and kappa arrays for x and y coordinates
    int npoints;
    Double_t det;

    vector<vector<Double_t>> xA; getKappa(kappaX, 0); getAArray(kappaX, xA);
    npoints=kappaX.size();
    TMatrixD xAtrans(npoints, npoints);
    for (int i=0; i<npoints; i++){ for (int j=0; j<npoints; j++){ xAtrans[i][j]=xA.at(i).at(j); } }
    det=0;
    TMatrixD xAinv = xAtrans.Invert(&det);
    if (det==0.){
      coutE(InputArguments) << "RooNCSplinePdf_2D_fast::interpolateFcn: Matrix xA could not be inverted. Something is wrong with the x coordinates of points!" << endl;
      assert(0);
    }

    vector<vector<Double_t>> yA; getKappa(kappaY, 1); getAArray(kappaY, yA);
    npoints=kappaY.size();
    TMatrixD yAtrans(npoints, npoints);
    for (int i=0; i<npoints; i++){ for (int j=0; j<npoints; j++){ yAtrans[i][j]=yA.at(i).at(j); } }
    det=0;
    TMatrixD yAinv = yAtrans.Invert(&det);
    if (det==0.){
      coutE(InputArguments) << "RooNCSplinePdf_2D_fast::interpolateFcn: Matrix yA could not be inverted. Something is wrong with the y coordinates of points!" << endl;
      assert(0);
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
        coutE(InputArguments) << "RooNCSplinePdf_2D_fast::interpolateFcn: nxbins!=(int)xcoefsAtYj.size() || npoldim!=(int)xcoefsAtYj.at(0).size()!" << endl;
        assert(0);
      }
      for (int ix=0; ix<nxbins; ix++){
        for (int ipow=0; ipow<npoldim; ipow++) coefsAlongY.at(ipow).at(ix).push_back(xcoefsAtYj.at(ix).at(ipow));
      }
    }

    for (int ix=0; ix<nxbins; ix++){
      // Get the x coefficients interpolated across y
      vector<vector<vector<Double_t>>> xCoefs;
      for (int ic=0; ic<npoldim; ic++){
        vector<vector<Double_t>> yCoefs = getCoefficientsAlongDirection(kappaY, yAinv, coefsAlongY.at(ic).at(ix), -1); // [iy][A,B,C,D]
        xCoefs.push_back(yCoefs);
      }
      coefficients.push_back(xCoefs);
    }
  }
  else assert(0);
}

RooNCSplinePdf_2D_fast::RooNCSplinePdf_2D_fast(
  const RooNCSplinePdf_2D_fast& other,
  const char* name
  ) :
  RooNCSplinePdf_2D(other, name),
  kappaX(other.kappaX),
  kappaY(other.kappaY)
{
  for (unsigned int i1=0; i1<other.coefficients.size(); i1++){
    vector<vector<vector<Double_t>>> coefarray;
    for (unsigned int i2=0; i2<other.coefficients.at(i1).size(); i2++){
      vector<vector<Double_t>> coefs;
      for (unsigned int i3=0; i3<other.coefficients.at(i1).at(i2).size(); i3++){
        vector<Double_t> coef;
        for (unsigned int i4=0; i4<other.coefficients.at(i1).at(i2).at(i3).size(); i4++) coef.push_back(other.coefficients.at(i1).at(i2).at(i3).at(i4));
        coefs.push_back(coef);
      }
      coefarray.push_back(coefs);
    }
    coefficients.push_back(coefarray);
  }
}


Double_t RooNCSplinePdf_2D_fast::interpolateFcn(Int_t code, const char* rangeName)const{
  Double_t res=0;

  if (verbosity==RooNCSplinePdfCore::kVerbose){ cout << "RooNCSplinePdf_2D_fast(" << GetName() << ")::interpolateFcn begin with code: " << code << endl; }

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

    Double_t txlow=0, txhigh=1;
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
    vector<Double_t> xCoefs;
    for (int ic=0; ic<(int)coefficients.at(ix).size(); ic++){
      const vector<vector<Double_t>>& yCoefs = coefficients.at(ix).at(ic);

      if (verbosity==RooNCSplinePdfCore::kVerbose) cout << "\tCoefficient " << ic << ":\n";

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

        if (verbosity==RooNCSplinePdfCore::kVerbose){
          if (code==0 || code%3!=0) cout << "\tEvaluating ty=" << tyhigh << " in bin " << iy << endl;
          else cout << "\tEvaluating ty[" << tylow << ", " << tyhigh << "] in bin " << iy << endl;
        }

        theCoef += evalSplineSegment(yCoefs.at(iy), kappaY.at(iy), tyhigh, tylow, (code>0 && code%3==0));
      }
      
      //if (code==0) cout << "\tCoefficient is " << theCoef << endl;
      
      xCoefs.push_back(theCoef);
    }

    // Evaluate value of spline at x with coefficients evaluated at y
    res += evalSplineSegment(xCoefs, kappaX.at(ix), txhigh, txlow, (code>0 && code%2==0));
  }

  return res;
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
