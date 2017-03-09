#include "RooNCSplinePdf_3D_fast.h" 
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


RooNCSplinePdf_3D_fast::RooNCSplinePdf_3D_fast() : RooNCSplinePdf_3D()
{}

RooNCSplinePdf_3D_fast::RooNCSplinePdf_3D_fast(
  const char* name,
  const char* title,
  RooAbsReal* inXVar,
  RooAbsReal* inYVar,
  RooAbsReal* inZVar,
  const RooArgList* inXList,
  const RooArgList* inYList,
  const RooArgList* inZList,
  std::vector<std::vector<const RooArgList*>>& inFcnList
  ) : RooNCSplinePdf_3D(name, title, inXVar, inYVar, inZVar, inXList, inYList, inZList, inFcnList)
{
  if (npointsX>1 && npointsY>1 && npointsZ>1){
    // Prepare A and kappa arrays for x, y and z coordinates
    int npoints;
    Double_t det;

    vector<vector<Double_t>> xA; getKappa(kappaX, 0); getAArray(kappaX, xA);
    npoints=kappaX.size();
    TMatrixD xAtrans(npoints, npoints);
    for (int i=0; i<npoints; i++){ for (int j=0; j<npoints; j++){ xAtrans[i][j]=xA.at(i).at(j); } }
    det=0;
    TMatrixD xAinv = xAtrans.Invert(&det);
    if (det==0.){
      coutE(InputArguments) << "RooNCSplinePdf_3D_fast::interpolateFcn: Matrix xA could not be inverted. Something is wrong with the x coordinates of points!" << endl;
      assert(0);
    }

    vector<vector<Double_t>> yA; getKappa(kappaY, 1); getAArray(kappaY, yA);
    npoints=kappaY.size();
    TMatrixD yAtrans(npoints, npoints);
    for (int i=0; i<npoints; i++){ for (int j=0; j<npoints; j++){ yAtrans[i][j]=yA.at(i).at(j); } }
    det=0;
    TMatrixD yAinv = yAtrans.Invert(&det);
    if (det==0.){
      coutE(InputArguments) << "RooNCSplinePdf_3D_fast::interpolateFcn: Matrix yA could not be inverted. Something is wrong with the y coordinates of points!" << endl;
      assert(0);
    }

    vector<vector<Double_t>> zA; getKappa(kappaZ, 2); getAArray(kappaZ, zA);
    npoints=kappaZ.size();
    TMatrixD zAtrans(npoints, npoints);
    for (int i=0; i<npoints; i++){ for (int j=0; j<npoints; j++){ zAtrans[i][j]=zA.at(i).at(j); } }
    det=0;
    TMatrixD zAinv = zAtrans.Invert(&det);
    if (det==0.){
      coutE(InputArguments) << "RooNCSplinePdf_3D_fast::interpolateFcn: Matrix zA could not be inverted. Something is wrong with the z coordinates of points!" << endl;
      assert(0);
    }

    //cout << "RooNCSplinePdf_3D_fast::RooNCSplinePdf_3D_fast: Initial kappa arrays are set up" << endl;

    // Get the grid of coefficients
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
      //cout << "Finding coefficients along z line " << k << endl;

      std::vector<std::vector<
        std::vector<std::vector<Double_t>>
        >> coefficients_perZ; // [ix][A_x,B_x,C_x,D_x][iy][A_x_y,B_x_y,C_x_y,D_x_y] at each z_k

      vector<vector<vector<Double_t>>> coefsAlongY; // [xbin][Ax(y),Bx(y),Cx(y),Dx(y)][ybin] in each z
      for (Int_t j=0; j<npointsY; j++){
        vector<vector<Double_t>> xcoefsAtYjZk = getCoefficientsPerYPerZ(kappaX, xAinv, j, k, -1); // [ix][Ax,Bx,Cx,Dx] at each y_j z_k
        //cout << "\tCoefficients in y line " << j << " are found" << endl;
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
          coutE(InputArguments) << "RooNCSplinePdf_3D_fast::interpolateFcn: nxbins!=(int)xcoefsAtYjZk.size() || nxpoldim!=(int)xcoefsAtYjZk.at(0).size()!" << endl;
          assert(0);
        }
        for (int ix=0; ix<nxbins; ix++){
          for (int icx=0; icx<nxpoldim; icx++) coefsAlongY.at(ix).at(icx).push_back(xcoefsAtYjZk.at(ix).at(icx));
        }
      }

      //cout << "\tCoefficients in each y line are found" << endl;

      for (int ix=0; ix<nxbins; ix++){
        // Get the x coefficients interpolated across y
        vector<vector<vector<Double_t>>> xCoefs;
        for (int icx=0; icx<nxpoldim; icx++){
          vector<vector<Double_t>> yCoefs = getCoefficientsAlongDirection(kappaY, yAinv, coefsAlongY.at(ix).at(icx), -1); // [iy][A,B,C,D]
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

    //cout << "Finding the final coefficients" << endl;
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

  }
  else assert(0);
}

RooNCSplinePdf_3D_fast::RooNCSplinePdf_3D_fast(
  const RooNCSplinePdf_3D_fast& other,
  const char* name
  ) :
  RooNCSplinePdf_3D(other, name),
  kappaX(other.kappaX),
  kappaY(other.kappaY),
  kappaZ(other.kappaZ)
{
  for (unsigned int i1=0; i1<other.coefficients.size(); i1++){
    vector<vector<vector<vector<vector<Double_t>>>>> c1;
    for (unsigned int i2=0; i2<other.coefficients.at(i1).size(); i2++){
      vector<vector<vector<vector<Double_t>>>> c2;
      for (unsigned int i3=0; i3<other.coefficients.at(i1).at(i2).size(); i3++){
        vector<vector<vector<Double_t>>> c3;
        for (unsigned int i4=0; i4<other.coefficients.at(i1).at(i2).at(i3).size(); i4++){
          vector<vector<Double_t>> c4;
          for (unsigned int i5=0; i5<other.coefficients.at(i1).at(i2).at(i3).at(i4).size(); i5++){
            vector<Double_t> c5;
            for (unsigned int i6=0; i6<other.coefficients.at(i1).at(i2).at(i3).at(i4).at(i5).size(); i6++){
              c5.push_back(other.coefficients.at(i1).at(i2).at(i3).at(i4).at(i5).at(i6));
            }
            c4.push_back(c5);
          }
          c3.push_back(c4);
        }
        c2.push_back(c3);
      }
      c1.push_back(c2);
    }
    coefficients.push_back(c1);
  }
}

Double_t RooNCSplinePdf_3D_fast::interpolateFcn(Int_t code, const char* rangeName)const{
  Double_t res=0;

  if (verbosity==RooNCSplinePdfCore::kVerbose){
    cout << "RooNCSplinePdf_3D_fast(" << GetName() << ")::interpolateFcn begin with code: " << code << endl;
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
    Int_t binval=-1;
    Int_t binmin=-1;
    Int_t binmax=-1;
    Float_t tval=0;
    Float_t tmin=0;
    Float_t tmax=0;
    if (code==0 || code%varprime.at(idim)!=0){
      binval = getWhichBin(*(varcoord.at(idim)), idim);
      tval = getTVar(*(varkappa.at(idim)), *(varcoord.at(idim)), binval, idim);
    }
    else{
      binmin = getWhichBin(varcoord.at(idim)->min(rangeName), idim);
      tmin = getTVar(*(varkappa.at(idim)), varcoord.at(idim)->min(rangeName), binmin, idim);
      binmax = getWhichBin(varcoord.at(idim)->max(rangeName), idim);
      tmax = getTVar(*(varkappa.at(idim)), varcoord.at(idim)->max(rangeName), binmax, idim);
    }
    varbin.push_back(binval);
    varbinmin.push_back(binmin);
    varbinmax.push_back(binmax);
    tvar.push_back(tval);
    tvarmin.push_back(tmin);
    tvarmax.push_back(tmax);
  }

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

Double_t RooNCSplinePdf_3D_fast::evaluate() const{
  Double_t value = interpolateFcn(0);
  if (value<=0.){
    coutE(Eval) << "RooNCSplinePdf_3D_fast ERROR::RooNCSplinePdf_3D_fast(" << GetName() << ") evaluation returned " << value << " at (x, y) = (" << theXVar << ", " << theYVar << ")" << endl;
    value = 1e-15;
  }
  if (verbosity==RooNCSplinePdfCore::kVerbose){ cout << "RooNCSplinePdf_3D_fast(" << GetName() << ")::evaluate = " << value << " at (x, y) = (" << theXVar << ", " << theYVar << ")" << endl; }
  return value;
}
Double_t RooNCSplinePdf_3D_fast::analyticalIntegral(Int_t code, const char* rangeName) const{
  Double_t value = interpolateFcn(code, rangeName);
  if (value<=0.){
    coutE(Integration) << "RooNCSplinePdf_3D_fast ERROR::RooNCSplinePdf_3D_fast(" << GetName() << ") integration returned " << value << " for code = " << code << endl;
    value = 1e-10;
  }
  if (verbosity==RooNCSplinePdfCore::kVerbose){ cout << "RooNCSplinePdf_3D_fast(" << GetName() << ")::analyticalIntegral = " << value << " for code = " << code << endl; }
  return value;
}

ClassImp(RooNCSplinePdf_3D_fast)
