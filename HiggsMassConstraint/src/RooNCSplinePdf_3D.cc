#include "RooNCSplinePdf_3D.h" 
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


RooNCSplinePdf_3D::RooNCSplinePdf_3D() :
RooNCSplinePdfCore(),
theYVar("theYVar", "theYVar", this),
theZVar("theZVar", "theZVar", this),
YList("YList", "List of Y coordinates", this),
ZList("ZList", "List of Z coordinates", this),
npointsY(0),
npointsZ(0)
{}

RooNCSplinePdf_3D::RooNCSplinePdf_3D(
  const char* name,
  const char* title,
  RooAbsReal* inXVar,
  RooAbsReal* inYVar,
  RooAbsReal* inZVar,
  const RooArgList* inXList,
  const RooArgList* inYList,
  const RooArgList* inZList,
  std::vector<std::vector<const RooArgList*>>& inFcnList,
  bool inUseConst
  ) :
  RooNCSplinePdfCore(name, title, inXVar, inXList, inUseConst),
  theYVar("theYVar", "theYVar", this, *inYVar),
  theZVar("theZVar", "theZVar", this, *inZVar),
  YList("YList", "List of Y coordinates", this),
  ZList("ZList", "List of Z coordinates", this)
{
  setProxyList(inYList, YList);
  setProxyList(inZList, ZList);
  npointsY = YList.getSize();
  npointsZ = ZList.getSize();

  for (unsigned int iz=0; iz<inFcnList.size(); iz++){
    vector<RooListProxy> dum;
    FcnList.push_back(dum);
    for (unsigned int iy=0; iy<inFcnList.at(iz).size(); iy++){
      const RooArgList* flist = inFcnList.at(iz).at(iy);
      FcnList.at(iz).emplace_back(Form("FcnListZ%iY%i", (int)iz, (int)iy), Form("List of function Z%iY%i coordinates", (int)iz, (int)iy), this);
      setProxyList(flist, FcnList.at(iz).at(iy));
    }
  }
}

RooNCSplinePdf_3D::RooNCSplinePdf_3D(
  const RooNCSplinePdf_3D& other,
  const char* name
  ) :
  RooNCSplinePdfCore(other, name),
  theYVar("theYVar", this, other.theYVar),
  theZVar("theZVar", this, other.theZVar),
  YList("YList", this, other.YList),
  ZList("ZList", this, other.ZList),
  npointsY(other.npointsY),
  npointsZ(other.npointsZ)
{
  for (unsigned int iz=0; iz<other.FcnList.size(); iz++){
    vector<RooListProxy> dum;
    FcnList.push_back(dum);
    for (unsigned int iy=0; iy<other.FcnList.at(iz).size(); iy++) FcnList.at(iz).emplace_back(Form("FcnListZ%iY%i", (int)iz, (int)iy), this, other.FcnList.at(iz).at(iy));
  }
}


Double_t RooNCSplinePdf_3D::interpolateFcn(Int_t code, const char* rangeName)const{
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
    coutE(InputArguments) << "RooNCSplinePdf_3D::interpolateFcn: Matrix xA could not be inverted. Something is wrong with the x coordinates of points!" << endl;
    assert(0);
  }

  vector<vector<Double_t>> yA; vector<Double_t> kappaY; getKappa(kappaY, 1); getAArray(kappaY, yA);
  npoints=kappaY.size();
  TMatrixD yAtrans(npoints, npoints);
  for (int i=0; i<npoints; i++){ for (int j=0; j<npoints; j++){ yAtrans[i][j]=yA.at(i).at(j); } }
  det=0;
  TMatrixD yAinv = yAtrans.Invert(&det);
  if (det==0.){
    coutE(InputArguments) << "RooNCSplinePdf_3D::interpolateFcn: Matrix yA could not be inverted. Something is wrong with the y coordinates of points!" << endl;
    assert(0);
  }

  vector<vector<Double_t>> zA; vector<Double_t> kappaZ; getKappa(kappaZ, 2); getAArray(kappaZ, zA);
  npoints=kappaZ.size();
  TMatrixD zAtrans(npoints, npoints);
  for (int i=0; i<npoints; i++){ for (int j=0; j<npoints; j++){ zAtrans[i][j]=zA.at(i).at(j); } }
  det=0;
  TMatrixD zAinv = zAtrans.Invert(&det);
  if (det==0.){
    coutE(InputArguments) << "RooNCSplinePdf_3D::interpolateFcn: Matrix zA could not be inverted. Something is wrong with the z coordinates of points!" << endl;
    assert(0);
  }

  // Get bins
  vector<Int_t> varprime; varprime.push_back(2); varprime.push_back(3); varprime.push_back(5);
  const Int_t ndims=(const Int_t)varprime.size();
  vector<Int_t> varbin;
  vector<Int_t> varbinmin;
  vector<Int_t> varbinmax;
  vector<Float_t> tvar;
  vector<Float_t> tvarmin;
  vector<Float_t> tvarmax;
  vector<const vector<Double_t>*> varkappa; varkappa.push_back(&kappaX); varkappa.push_back(&kappaY); varkappa.push_back(&kappaZ);
  vector<const RooRealProxy*> varcoord; varcoord.push_back(&theXVar); varcoord.push_back(&theYVar); varcoord.push_back(&theZVar);
  for (Int_t idim=0; idim<ndims; idim++){
    varbin.at(idim)=-1;
    varbinmin.at(idim)=-1;
    varbinmax.at(idim)=-1;
    tvar.at(idim)=0;
    tvarmin.at(idim)=0;
    tvarmax.at(idim)=0;
    if (code==0 || code%varprime.at(idim)!=0){
      varbin.at(idim) = getWhichBin(*(varcoord.at(idim)), idim);
      tvar.at(idim) = getTVar(*(varkappa.at(idim)), *(varcoord.at(idim)), varbin.at(idim), idim);
    }
    else{
      varbinmin.at(idim) = getWhichBin(varcoord.at(idim)->min(rangeName), idim);
      tvarmin.at(idim) = getTVar(*(varkappa.at(idim)), varcoord.at(idim)->min(rangeName), varbinmin.at(idim), idim);
      varbinmax.at(idim) = getWhichBin(varcoord.at(idim)->max(rangeName), idim);
      tvarmax.at(idim) = getTVar(*(varkappa.at(idim)), varcoord.at(idim)->max(rangeName), varbinmax.at(idim), idim);
    }
  }

  // Get the grid of coefficients
  std::vector<std::vector<
    std::vector<std::vector<
    std::vector<std::vector<Double_t>>
    >>
    >> coefficients; // [ix][A_x,B_x,C_x,D_x][iy][A_x_y,B_x_y,C_x_y,D_x_y][iz][A_x_y_z,B_x_y_z,C_x_y_z,D_x_y_z]
  int nxpoldim=0;
  int nxbins=0;
  int nybins=0;
  int nypoldim=0;
  std::vector<
    std::vector<std::vector<
    std::vector<std::vector<Double_t>>
    >>
    > coefsAlongZ; // [ix][A_x,B_x,C_x,D_x][iy][A_x_y,B_x_y,C_x_y,D_x_y][z_k] at each z_k
  for (Int_t k=0; k<npointsZ; k++){
    std::vector<std::vector<
      std::vector<std::vector<Double_t>>
      >> coefficients_perZ; // [ix][A_x,B_x,C_x,D_x][iy][A_x_y,B_x_y,C_x_y,D_x_y] at each z_k

    vector<vector<vector<Double_t>>> coefsAlongY; // [Ax(y),Bx(y),Cx(y),Dx(y)][xbin][ybin] in each z
    for (Int_t j=0; j<npointsY; j++){
      vector<vector<Double_t>> xcoefsAtYjZk = getCoefficientsPerYPerZ(kappaX, xAinv, j, k, -1); // [ix][Ax,Bx,Cx,Dx] at each y_j z_k
      if (j==0){
        if (k==0){
          nxbins=xcoefsAtYjZk.size();
          nxpoldim=xcoefsAtYjZk.at(0).size();
        }
        for (int ix=0; ix<nxbins; ix++){
          vector<vector<Double_t>> dum_xycoefarray;
          for (int icx=0; icx<nxpoldim; icx++){
            vector<Double_t> dum_ycoefarray;
            dum_xycoefarray.push_back(dum_ycoefarray);
          }
          coefsAlongY.push_back(dum_xycoefarray);
        }
      }
      if (nxbins!=(int)xcoefsAtYjZk.size() || nxpoldim!=(int)xcoefsAtYjZk.at(0).size()){
        coutE(InputArguments) << "RooNCSplinePdf_3D::interpolateFcn: nxbins!=(int)xcoefsAtYjZk.size() || nxpoldim!=(int)xcoefsAtYjZk.at(0).size()!" << endl;
        assert(0);
      }
      for (int ix=0; ix<nxbins; ix++){
        for (int icx=0; icx<nxpoldim; icx++) coefsAlongY.at(ix).at(icx).push_back(xcoefsAtYjZk.at(ix).at(icx));
      }
    }

    for (int ix=0; ix<nxbins; ix++){
      // Get the x coefficients interpolated across y
      vector<vector<vector<Double_t>>> xCoefs;
      for (int icx=0; icx<nxpoldim; icx++){
        vector<vector<Double_t>> yCoefs = getCoefficientsAlongDirection(kappaY, yAinv, coefsAlongY.at(icx).at(ix), -1); // [iy][A,B,C,D]
        xCoefs.push_back(yCoefs);
      }
      coefficients_perZ.push_back(xCoefs);
    }

    if (k==0){
      nybins = coefficients_perZ.at(0).at(0).size();
      nypoldim = coefficients_perZ.at(0).at(0).at(0).size();
      for (int ix=0; ix<nxbins; ix++){
        vector<vector<vector<vector<Double_t>>>> xCoefs;
        for (int icx=0; icx<nxpoldim; icx++){
          vector<vector<vector<Double_t>>> xCoefsAlongY;
          for (int iy=0; iy<nybins; iy++){
            vector<vector<Double_t>> yCoefs;
            for (int icy=0; icy<nypoldim; icy++){
              vector<Double_t> yCoefAtZj;
              yCoefs.push_back(yCoefAtZj);
            }
            xCoefsAlongY.push_back(yCoefs);
          }
          xCoefs.push_back(xCoefsAlongY);
        }
        coefsAlongZ.push_back(xCoefs);
      }
    }
    for (int ix=0; ix<nxbins; ix++){
      for (int icx=0; icx<nxpoldim; icx++){
        for (int iy=0; iy<nybins; iy++){
          for (int icy=0; icy<nypoldim; icy++){
            coefsAlongZ.at(ix).at(icx).at(iy).at(icy).push_back(coefficients_perZ.at(ix).at(icx).at(iy).at(icy));
          }
        }
      }
    }
  }

  for (int ix=0; ix<nxbins; ix++){
    vector<vector<vector<vector<vector<Double_t>>>>> xCoefs;
    for (int icx=0; icx<nxpoldim; icx++){
      vector<vector<vector<vector<Double_t>>>> xCoefsAlongY;
      for (int iy=0; iy<nybins; iy++){
        vector<vector<vector<Double_t>>> yCoefs;
        for (int icy=0; icy<nypoldim; icy++){
          vector<vector<Double_t>> yCoefsAlongZ = getCoefficientsAlongDirection(kappaZ, zAinv, coefsAlongZ.at(ix).at(icx).at(iy).at(icy), -1); // [iz][A,B,C,D]
          yCoefs.push_back(yCoefsAlongZ);
        }
        xCoefsAlongY.push_back(yCoefs);
      }
      xCoefs.push_back(xCoefsAlongY);
    }
    coefficients.push_back(xCoefs);
  }

  // Evaluate over the coefficiencts
  for (int ix=0; ix<(int)coefficients.size(); ix++){
    if (
      (varbin[0]>=0 && ix!=varbin[0])
      ||
      (varbinmin[0]>=0 && varbinmax[0]>=varbinmin[0] && !(varbinmin[0]<=ix && ix<=varbinmax[0]))
      ) continue;
    Double_t txlow=0, txhigh=1;
    if (code>0 && code%varprime[0]==0){
      if (ix==varbinmin[0]) txlow=tvarmin[0];
      if (ix==varbinmax[0]) txhigh=tvarmax[0];
    }
    else txhigh=tvar[0];
    // Get the x coefficients interpolated across y
    vector<Double_t> xCoefs;
    for (int icx=0; icx<(int)coefficients.at(ix).size(); icx++){
      Double_t theXCoef=0;
      for (int iy=0; iy<(int)coefficients.at(ix).at(icx).size(); iy++){
        if (
          (varbin[1]>=0 && iy!=varbin[1])
          ||
          (varbinmin[1]>=0 && varbinmax[1]>=varbinmin[1] && !(varbinmin[1]<=iy && iy<=varbinmax[1]))
          ) continue;
        Double_t tylow=0, tyhigh=1;
        if (code>0 && code%varprime[1]==0){
          if (iy==varbinmin[1]) tylow=tvarmin[1];
          if (iy==varbinmax[1]) tyhigh=tvarmax[1];
        }
        else tyhigh=tvar[1];
        // Get the y coefficients interpolated across z
        vector<Double_t> yCoefs;
        for (int icy=0; icy<(int)coefficients.at(ix).at(icx).at(iy).size(); icy++){
          Double_t theYCoef=0;
          for (int iz=0; iz<(int)coefficients.at(ix).at(icx).at(iy).at(icy).size(); iz++){
            if (
              (varbin[2]>=0 && iz!=varbin[2])
              ||
              (varbinmin[2]>=0 && varbinmax[2]>=varbinmin[2] && !(varbinmin[2]<=iz && iz<=varbinmax[2]))
              ) continue;

            Double_t tzlow=0, tzhigh=1;
            if (code>0 && code%varprime[2]==0){
              if (iz==varbinmin[2]) tzlow=tvarmin[2];
              if (iz==varbinmax[2]) tzhigh=tvarmax[2];
            }
            else tzhigh=tvar[2];

            theYCoef += evalSplineSegment(coefficients.at(ix).at(icx).at(iy).at(icy).at(iz), varkappa[2]->at(iz), tzhigh, tzlow, (code>0 && code%varprime[2]==0));
          }
          yCoefs.push_back(theYCoef);
        }
        theXCoef += evalSplineSegment(yCoefs, varkappa[1]->at(iy), tyhigh, tylow, (code>0 && code%varprime[1]==0));
      }
      xCoefs.push_back(theXCoef);
    }
    // Evaluate value of spline at x with coefficients evaluated at y
    res += evalSplineSegment(xCoefs, varkappa[0]->at(ix), txhigh, txlow, (code>0 && code%varprime[0]==0));
  }

  return res;
}

void RooNCSplinePdf_3D::getKappa(vector<Double_t>& kappas, const Int_t whichDirection)const{
  kappas.clear();
  Double_t kappa=1;

  Int_t npoints;
  RooListProxy const* coord;
  if (whichDirection==0){
    npoints=npointsX;
    coord=&XList;
  }
  else if (whichDirection==1){
    npoints=npointsY;
    coord=&YList;
  }
  else{
    npoints=npointsZ;
    coord=&ZList;
  }
  for (Int_t j=0; j<npoints-1; j++){
    Double_t val_j = (dynamic_cast<RooAbsReal*>(coord->at(j)))->getVal();
    Double_t val_jpo = (dynamic_cast<RooAbsReal*>(coord->at(j+1)))->getVal();
    kappa = 1./(val_jpo-val_j);
    kappas.push_back(kappa);
  }
  kappas.push_back(kappa); // Push the same kappa_(N-1)=kappa_(N-2) at the end point
}
Int_t RooNCSplinePdf_3D::getWhichBin(const Double_t& val, const Int_t whichDirection)const{
  //cout << "RooNCSplinePdf_3D::getWhichBin: begin with val = " << val << " dir = " << whichDirection << endl;
  Int_t bin=-1;
  Double_t valj, valjpo;
  Int_t npoints;
  RooListProxy const* coord;
  if (whichDirection==0){
    npoints=npointsX;
    //cout << "RooNCSplinePdf_3D::getWhichBin: npoints=" << npoints << endl;
    coord=&XList;
    //cout << "RooNCSplinePdf_3D::getWhichBin: coord=" << coord << endl;
  }
  else if (whichDirection==1){
    npoints=npointsY;
    //cout << "RooNCSplinePdf_3D::getWhichBin: npoints=" << npoints << endl;
    coord=&YList;
    //cout << "RooNCSplinePdf_3D::getWhichBin: coord=" << coord << endl;
  }
  else{
    npoints=npointsZ;
    //cout << "RooNCSplinePdf_3D::getWhichBin: npoints=" << npoints << endl;
    coord=&ZList;
    //cout << "RooNCSplinePdf_3D::getWhichBin: coord=" << coord << endl;
  }
  //if (coord!=0) cout << "RooNCSplinePdf_3D::getWhichBin: coordN=" << coord->getSize() << endl;

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
Double_t RooNCSplinePdf_3D::getTVar(const vector<Double_t>& kappas, const Double_t& val, const Int_t& bin, const Int_t whichDirection)const{
  Double_t K;
  K=kappas.at(bin);
  RooListProxy const* coord;
  if (whichDirection==0) coord=&XList;
  else if (whichDirection==1) coord=&YList;
  else coord=&ZList;
  return (Double_t)((val-(dynamic_cast<RooAbsReal*>(coord->at(bin)))->getVal())*K);
}

vector<vector<Double_t>> RooNCSplinePdf_3D::getCoefficientsPerYPerZ(
  const std::vector<Double_t>& kappaX, const TMatrixD& xAinv,
  const Int_t& ybin, const Int_t& zbin,
  const Int_t xbin
  )const{
  vector<Double_t> fcnList;
  for (int bin=0; bin<npointsX; bin++){ fcnList.push_back((dynamic_cast<RooAbsReal*>(FcnList.at(zbin).at(ybin).at(bin)))->getVal()); }
  vector<vector<Double_t>> coefs = getCoefficientsAlongDirection(kappaX, xAinv, fcnList, xbin);
  return coefs;
}

Double_t RooNCSplinePdf_3D::evaluate() const{
  Double_t value = interpolateFcn(0);
  if (value<=0.) return 1e-15;
  return value;
}
Int_t RooNCSplinePdf_3D::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const{
  Int_t code=1;
  if (dynamic_cast<RooRealVar*>(theXVar.absArg())!=0){
    if (matchArgs(allVars, analVars, theXVar)) code*=2;
  }
  if (dynamic_cast<RooRealVar*>(theYVar.absArg())!=0){
    if (matchArgs(allVars, analVars, theYVar)) code*=3;
  }
  if (dynamic_cast<RooRealVar*>(theZVar.absArg())!=0){
    if (matchArgs(allVars, analVars, theZVar)) code*=5;
  }
  if (code==1) code=0;
  return code;
}
Double_t RooNCSplinePdf_3D::analyticalIntegral(Int_t code, const char* rangeName) const{
  Double_t value = interpolateFcn(code, rangeName);
  if (value<=0.) value = 1e-10;
  return value;
}

ClassImp(RooNCSplinePdf_3D)
