#include <cassert> 
#include "RooGaussianMomConstraint.h"

using namespace TMath;
using namespace std;


RooGaussianMomConstraint::RooGaussianMomConstraint(
  const char* name, const char* title,
  const RooArgList& variables_,
  const RooArgList& means_,
  const RooArgList& matrixElement_,
  RooGaussianMomConstraint::CoordinateSystem coordinates_,
  RooGaussianMomConstraint::VerbosityLevel verbosity_
  ) : RooAbsPdf(name, title),

  variables("variables", "List of variables", this),
  means("means", "List of means", this),
  matrixElement("matrixElement", "List of inverse covariance matrix elements", this),
  coordinates(coordinates_),
  verbosity(verbosity_)
{
  setProxyList(variables_, variables, 3);
  setProxyList(means_, means, 3);
  setProxyList(matrixElement_, matrixElement, 3*3);
}


RooGaussianMomConstraint::RooGaussianMomConstraint(const RooGaussianMomConstraint& other, const char* name) : RooAbsPdf(other, name),
variables("variables", "List of variables", this),
means("means", "List of means", this),
matrixElement("matrixElement", "List of inverse covariance matrix elements", this),
coordinates(other.coordinates),
verbosity(other.verbosity)
{
  setProxyList(other.variables, variables, 3);
  setProxyList(other.means, means, 3);
  setProxyList(other.matrixElement, matrixElement, 3*3);
}

void RooGaussianMomConstraint::setVerbosity(RooGaussianMomConstraint::VerbosityLevel verbosity_){ verbosity=verbosity_; }

void RooGaussianMomConstraint::setProxyList(const RooArgList& args, RooListProxy& target, Int_t checkDim){
  TIterator* argIter = args.createIterator();
  RooAbsArg* arg;
  Int_t nargs=0;
  while ((arg = (RooAbsArg*)argIter->Next())){
    if (nargs==1) break;
    if (arg!=0){
      if (!dynamic_cast<RooRealVar*>(arg)){
        coutE(InputArguments) << "ERROR::RooGaussianMomConstraint(" << GetName() << ")::setProxyList: Proxy " << arg->GetName() << " is not of type RooRealVar." << endl;
        assert(0);
      }
      target.add(*arg);
      nargs++;
    }
  }
  delete argIter;
  /*
  if (checkDim>0 && checkDim!=nargs){
    coutE(InputArguments) << "ERROR::RooGaussianMomConstraint(" << GetName() << ")::setProxyList: Proxies (size=" << nargs << ") do not have the required size " << checkDim << "." << endl;
    assert(0);
  }
  */
  if (nargs==0){
    coutE(InputArguments) << "ERROR::RooGaussianMomConstraint(" << GetName() << ")::setProxyList: Number of proxy arguments are 0!" << endl;
    assert(0);
  }
}


Double_t RooGaussianMomConstraint::computeGaussian(const Int_t code) const{ // The code can be either for var1, 2 or 3, but it can only be for one of them.
  if (coordinates==RooGaussianMomConstraint::kXYZ) return computeCaseXYZ(code);
  else if (coordinates==RooGaussianMomConstraint::kRhoLambdaPhi) return computeCaseRhoLambdaPhi(code);
  else if (code==1) return 1e-15;
  else return 1e-10;
}


Double_t RooGaussianMomConstraint::evaluate() const{
  if (verbosity==RooGaussianMomConstraint::kVerbose)
    cout << "Begin RooGaussianMomConstraint::evaluate in PDF " << GetName() << endl;

  Double_t result = computeGaussian(1);
  if (verbosity==RooGaussianMomConstraint::kVerbose)
    cout << "End RooGaussianMomConstraint::evaluate in PDF " << GetName() << '\n' << endl;

  return result;
}

Double_t RooGaussianMomConstraint::analyticalIntegral(Int_t code, const char* /*rangeName*/) const{
  if (verbosity==RooGaussianMomConstraint::kVerbose)
    cout << "Begin RooGaussianMomConstraint::analyticalIntegral in PDF " << GetName() << endl;

  Double_t result = computeGaussian(code);
  if (verbosity==RooGaussianMomConstraint::kVerbose)
    cout << "End RooGaussianMomConstraint::analyticalIntegral in PDF " << GetName() << '\n' << endl;

  return result;
}

Int_t RooGaussianMomConstraint::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const{
  Int_t code = 1;

  if (verbosity==RooGaussianMomConstraint::kVerbose){
    cout << "Begin RooGaussianMomConstraint::getAnalyticalIntegral." << endl;
    cout << "\tInitial allVars : " << endl;
    allVars.Print("v");
    cout << "\tInitial analVars: " << endl;
    analVars.Print("v");
  }

  if (coordinates==RooGaussianMomConstraint::kXYZ){
    if (matchArgs(allVars, analVars, RooArgSet(*(variables.at(0))))) code *= prime_var1;
    //else if (matchArgs(allVars, analVars, RooArgSet(*(variables.at(1))))) code *= prime_var2;
    //else if (matchArgs(allVars, analVars, RooArgSet(*(variables.at(2))))) code *= prime_var3;
    else if (matchArgs(allVars, analVars, RooArgSet(*(means.at(0))))) code *= prime_mean1;
    //else if (matchArgs(allVars, analVars, RooArgSet(*(means.at(1))))) code *= prime_mean2;
    //else if (matchArgs(allVars, analVars, RooArgSet(*(means.at(2))))) code *= prime_mean3;
  }
  else if (coordinates==RooGaussianMomConstraint::kRhoLambdaPhi){
    if (matchArgs(allVars, analVars, RooArgSet(*(variables.at(0))))) code *= prime_var1;
    //else if (matchArgs(allVars, analVars, RooArgSet(*(variables.at(1))))) code *= prime_var2; // No analytical integration over lambda possible; an implementation of Erf with complex arguments is needed!
    //else if (matchArgs(allVars, analVars, RooArgSet(*(variables.at(2))))) code *= prime_var3;
    else if (matchArgs(allVars, analVars, RooArgSet(*(means.at(0))))) code *= prime_mean1;
    //else if (matchArgs(allVars, analVars, RooArgSet(*(means.at(1))))) code *= prime_mean2; // No analytical integration over lambda possible; an implementation of Erf with complex arguments is needed!
    //else if (matchArgs(allVars, analVars, RooArgSet(*(means.at(2))))) code *= prime_mean3;
  }
  if (code==1) code=0;

  if (verbosity==RooGaussianMomConstraint::kVerbose){
    cout << "\tFinal allVars: " << endl;
    allVars.Print("v");
    cout << "\tFinal analVars: " << endl;
    analVars.Print("v");
    cout << "\tCode: " << code << endl;
    cout << "End RooGaussianMomConstraint::getAnalyticalIntegral." << endl;
  }

  return code;
}


Double_t RooGaussianMomConstraint::computeCaseXYZ(const Int_t code) const{ // The code can be either for var1, 2 or 3, but it can only be for one of them.
  Int_t intVar=-1;
  const int nDims=1;
  if (code>0){
    if ((code%prime_var1)==0)intVar=0;
    else if ((code%prime_var2)==0)intVar=1;
    else if ((code%prime_var3)==0)intVar=2;
    else if ((code%prime_mean1)==0)intVar=3;
    else if ((code%prime_mean2)==0)intVar=4;
    else if ((code%prime_mean3)==0)intVar=5;
  }
  if (verbosity==RooGaussianMomConstraint::kVerbose)
    cout << "Begin RooGaussianMomConstraint::computeCaseXYZ with code = " << code << " and intVar = " << intVar << endl;

  Double_t value=0;
  if (intVar<0){ // Just do regular computation
    const Double_t epsilon = 1e-15;

    if (verbosity==RooGaussianMomConstraint::kVerbose)
      cout << "RooGaussianMomConstraint::computeCaseXYZ case 0" << endl;

    for (int ix=0; ix<nDims; ix++){
      /*
      if (
      (((fixCode%prime_var1)==0 || (fixCode%prime_mean1)==0) && ix==0) ||
      (((fixCode%prime_var2)==0 || (fixCode%prime_mean2)==0) && ix==1) ||
      (((fixCode%prime_var3)==0 || (fixCode%prime_mean3)==0) && ix==2)
      ) continue;
      */
      Double_t val_ix = dynamic_cast<const RooRealVar*>(variables.at(ix))->getVal();
      Double_t valbar_ix = dynamic_cast<const RooRealVar*>(means.at(ix))->getVal();
      Double_t diff_ix=val_ix-valbar_ix;

      if (verbosity==RooGaussianMomConstraint::kVerbose){
        cout << "ix = " << ix << endl;
        cout << "x, xb = " << val_ix << " " << valbar_ix << endl;
      }

      for (int iy=0; iy<nDims; iy++){
        /*
        if (
        (((fixCode%prime_var1)==0 || (fixCode%prime_mean1)==0) && iy==0) ||
        (((fixCode%prime_var2)==0 || (fixCode%prime_mean2)==0) && iy==1) ||
        (((fixCode%prime_var3)==0 || (fixCode%prime_mean3)==0) && iy==2)
        ) continue;
        */
        Double_t invCovMat = dynamic_cast<const RooRealVar*>(matrixElement.at(3*ix+iy))->getVal();

        if (verbosity==RooGaussianMomConstraint::kVerbose){
          cout << "iy = " << iy << endl;
          cout << "invCovMat(ix, iy) = " << invCovMat << endl;
        }

        if (invCovMat!=0.){
          Double_t val_iy = dynamic_cast<const RooRealVar*>(variables.at(iy))->getVal();
          Double_t valbar_iy = dynamic_cast<const RooRealVar*>(means.at(iy))->getVal();
          Double_t diff_iy=val_iy-valbar_iy;
          value += -0.5*invCovMat*diff_ix*diff_iy;

          if (verbosity==RooGaussianMomConstraint::kVerbose)
            cout << "y, yb = " << val_iy << " " << valbar_iy << endl;
        }
        if (verbosity==RooGaussianMomConstraint::kVerbose)
          cout << endl;
      }
    }
    value = exp(value);

    if (verbosity==RooGaussianMomConstraint::kVerbose)
      cout << "RooGaussianMomConstraint::computeCaseXYZ intVar==" << intVar << " result=" << value << " epsilon=" << epsilon << endl;

    if (!(value==value) || value<epsilon) return epsilon;
  }
  else{
    const Double_t epsilon = 1e-10;

    if (verbosity==RooGaussianMomConstraint::kVerbose)
      cout << "RooGaussianMomConstraint::computeCaseXYZ case 1" << endl;

    Int_t intVarpr=intVar%3;
    //Double_t xmin, xmax, y, z, xb, yb, zb, dxmax, dxmin, dy, dz;
    Double_t xmin, xmax, xb, dxmax, dxmin, dy, dz;
    const RooListProxy* vararray;
    const RooListProxy* meanarray;
    if (intVar<3){
      vararray=&variables;
      meanarray=&means;
    }
    else{
      meanarray=&variables;
      vararray=&means;
    }
    xmax = dynamic_cast<const RooRealVar*>(vararray->at(intVarpr))->getMax();
    xmin = dynamic_cast<const RooRealVar*>(vararray->at(intVarpr))->getMin();
    //y = dynamic_cast<const RooRealVar*>(vararray->at((intVarpr+1)%3))->getVal();
    //z = dynamic_cast<const RooRealVar*>(vararray->at((intVarpr+2)%3))->getVal();
    xb = dynamic_cast<const RooRealVar*>(meanarray->at(intVarpr))->getVal();
    //yb = dynamic_cast<const RooRealVar*>(meanarray->at((intVarpr+1)%3))->getVal();
    //zb = dynamic_cast<const RooRealVar*>(meanarray->at((intVarpr+2)%3))->getVal();
    dxmax = xmax-xb;
    dxmin = xmin-xb;
    //dy = y-yb;
    //dz = z-zb;
    dy=0; dz=0;

    if (verbosity==RooGaussianMomConstraint::kVerbose){
      cout << "xb, xmax, xmin = " << xb << ", " << xmax << ", " << xmin << endl;
      //cout << "y, yb = " << y << ", " << yb << endl;
      //cout << "z, zb = " << z << ", " << zb << endl;
    }

    Double_t invCovMat[3][3]={ { 0 } };
    for (int ix=0; ix<nDims; ix++){
      Int_t rx = ((intVarpr+ix)%3);
      /*
      if (
      (((fixCode%prime_var1)==0 || (fixCode%prime_mean1)==0) && rx==0) ||
      (((fixCode%prime_var2)==0 || (fixCode%prime_mean2)==0) && rx==1) ||
      (((fixCode%prime_var3)==0 || (fixCode%prime_mean3)==0) && rx==2)
      ) continue;
      */
      for (int iy=0; iy<nDims; iy++){
        Int_t ry = ((intVarpr+iy)%3);
        /*
        if (
          (((fixCode%prime_var1)==0 || (fixCode%prime_mean1)==0) && ry==0) ||
          (((fixCode%prime_var2)==0 || (fixCode%prime_mean2)==0) && ry==1) ||
          (((fixCode%prime_var3)==0 || (fixCode%prime_mean3)==0) && ry==2)
          ) continue;
        */
        invCovMat[rx][ry] = dynamic_cast<const RooRealVar*>(matrixElement.at(3*rx+ry))->getVal();
        if (verbosity==RooGaussianMomConstraint::kVerbose)
          cout << "invCovMat(" << rx << ", " << ry << ") = " << invCovMat[rx][ry] << endl;
      }
    }

    Double_t offxdiff = ((invCovMat[0][1]+invCovMat[1][0])*dy + (invCovMat[0][2]+invCovMat[2][0])*dz)/2.;
    Double_t diffothers = invCovMat[1][1]*pow(dy, 2) + (invCovMat[1][2]+invCovMat[2][1])*dy*dz + invCovMat[2][2]*pow(dz, 2);
    const Double_t tolerance = epsilon*0.01;
    if (fabs(invCovMat[0][0])>tolerance){
      Double_t exponent = 0.5*(
        pow(offxdiff, 2)/invCovMat[0][0]
        - diffothers
        );

      value = 1./sqrt(invCovMat[0][0]) * sqrt(Pi()/2.) * exp(exponent) * (
        Erf((invCovMat[0][0]*dxmax + offxdiff)/sqrt(2.*invCovMat[0][0]))
        - Erf((invCovMat[0][0]*dxmin + offxdiff)/sqrt(2.*invCovMat[0][0]))
        );

      if (verbosity==RooGaussianMomConstraint::kVerbose)
        cout << "RooGaussianMomConstraint::computeCaseXYZ case 1 exponent=" << exponent << " value=" << value << endl;
    }
    else if (fabs(offxdiff)>tolerance){
      value =
        exp(
        -0.5*diffothers
        )*(
        exp(-offxdiff*dxmin) - exp(-offxdiff*dxmax)
        )/offxdiff;

      if (verbosity==RooGaussianMomConstraint::kVerbose)
        cout << "RooGaussianMomConstraint::computeCaseXYZ case 2 value=" << value << endl;
    }
    else{
      value =
        exp(
        -0.5*diffothers
        )*(xmax-xmin);

      if (verbosity==RooGaussianMomConstraint::kVerbose)
        cout << "RooGaussianMomConstraint::computeCaseXYZ case 3 value=" << value << endl;
    }

    if (verbosity==RooGaussianMomConstraint::kVerbose)
      cout << "RooGaussianMomConstraint::computeCaseXYZ code " << code << " with result " << value << endl;

    if (!(value==value) || value<epsilon) return epsilon;
  }

  if (verbosity==RooGaussianMomConstraint::kVerbose)
    cout << "End RooGaussianMomConstraint::computeCaseXYZ with result " << value << endl;

  return value;
}
Double_t RooGaussianMomConstraint::computeCaseRhoLambdaPhi(const Int_t code) const{ // The code can be either for var1, (not)2 or 3, but it can only be for one of them.
  Int_t intVar=-1;
  const int nDims=1;
  if (code>0){
    if ((code%prime_var1)==0)intVar=0;
    //else if ((code%prime_var2)==0)intVar=1; // No idea how to implement integration over the Lambda coordinate.
    else if ((code%prime_var2)==0){ cerr << "RooGaussianMomConstraint::computeCaseRhoLambdaPhi: Something went terribly wrong! code%prime_var2==0 should not have happened!" << endl; assert(0); }
    else if ((code%prime_var3)==0)intVar=2;
    else if ((code%prime_mean1)==0)intVar=3;
    //else if ((code%prime_mean2)==0)intVar=4; // No idea how to implement integration over the Lambda coordinate.
    else if ((code%prime_mean2)==0){ cerr << "RooGaussianMomConstraint::computeCaseRhoLambdaPhi: Something went terribly wrong! code%prime_mean2==0 should not have happened!" << endl; assert(0); }
    else if ((code%prime_mean3)==0)intVar=5;
  }

  if (verbosity==RooGaussianMomConstraint::kVerbose)
    cout << "Begin RooGaussianMomConstraint::computeCaseRhoLambdaPhi with code = " << code << " and intVar = " << intVar << endl;

  Double_t value=0;
  const Double_t epsilon = 1e-10;
  if (intVar<0 || (code%prime_var2)==0 || (code%prime_mean2)==0){ // Only extra term is (Rho Sec[Lambda])^2, so numerical value with either no integration or integration over Phi is unchanged other than this "constant" value.
    if (verbosity==RooGaussianMomConstraint::kVerbose)
      cout << "RooGaussianMomConstraint::computeCaseRhoLambdaPhi case 0" << endl;

    value = computeCaseXYZ(code); // Just do regular computation

    const RooListProxy* vararray;
    if (intVar<3) vararray=&variables;
    else vararray=&means;

    Double_t x = dynamic_cast<const RooRealVar*>(vararray->at(0))->getVal();
    Double_t y;
    //y = dynamic_cast<const RooRealVar*>(vararray->at(1))->getVal();
    y = 0;
    Double_t jacobian = pow(x, 2);
    if (fabs(cos(y))<sqrt(epsilon)) jacobian /= epsilon;
    else jacobian /= pow(cos(y), 2);
    value *= jacobian;

    if (verbosity==RooGaussianMomConstraint::kVerbose){
      cout << "RooGaussianMomConstraint::computeCaseRhoLambdaPhi x, y, jacobian " << x << ", " << y << ", " << jacobian << endl;
      cout << "RooGaussianMomConstraint::computeCaseRhoLambdaPhi code " << code << " with result " << value << endl;
    }

    if (!(value==value) || value<epsilon) value = epsilon;
    return value;
  }
  // return is called, so else-if does not make sense!

  Int_t intVarpr=intVar%3;
  if (intVarpr==0){ // Integration over cylindrical Rho is different wrt. how it would have been in the XYZ coordinates.
    if (verbosity==RooGaussianMomConstraint::kVerbose)
      cout << "RooGaussianMomConstraint::computeCaseRhoLambdaPhi case 1" << endl;

    //Double_t xmin, xmax, y, z, xb, yb, zb, dxmax, dxmin, dy, dz;
    Double_t xmin, xmax, y, xb, dxmax, dxmin, dy, dz;
    const RooListProxy* vararray;
    const RooListProxy* meanarray;
    if (intVar<3){
      vararray=&variables;
      meanarray=&means;
    }
    else{
      meanarray=&variables;
      vararray=&means;
    }
    xmax = dynamic_cast<const RooRealVar*>(vararray->at(intVarpr))->getMax();
    xmin = dynamic_cast<const RooRealVar*>(vararray->at(intVarpr))->getMin();
    //y = dynamic_cast<const RooRealVar*>(vararray->at((intVarpr+1)%3))->getVal();
    //z = dynamic_cast<const RooRealVar*>(vararray->at((intVarpr+2)%3))->getVal();
    xb = dynamic_cast<const RooRealVar*>(meanarray->at(intVarpr))->getVal();
    //yb = dynamic_cast<const RooRealVar*>(meanarray->at((intVarpr+1)%3))->getVal();
    //zb = dynamic_cast<const RooRealVar*>(meanarray->at((intVarpr+2)%3))->getVal();
    dxmax = xmax-xb;
    dxmin = xmin-xb;
    //dy = y-yb;
    //dz = z-zb;
    y=0; dy=0; dz=0;

    if (verbosity==RooGaussianMomConstraint::kVerbose){
      cout << "xb, xmax, xmin = " << xb << ", " << xmax << ", " << xmin << endl;
      //cout << "y, yb = " << y << ", " << yb << endl;
      //cout << "z, zb = " << z << ", " << zb << endl;
    }

    Double_t invCovMat[3][3]={ { 0 } };
    for (int ix=0; ix<nDims; ix++){
      Int_t rx = ((intVarpr+ix)%3);
      /*
      if (
        (((fixCode%prime_var1)==0 || (fixCode%prime_mean1)==0) && rx==0) ||
        (((fixCode%prime_var2)==0 || (fixCode%prime_mean2)==0) && rx==1) ||
        (((fixCode%prime_var3)==0 || (fixCode%prime_mean3)==0) && rx==2)
        ) continue;
      */
      for (int iy=0; iy<nDims; iy++){
        Int_t ry = ((intVarpr+iy)%3);
        /*
        if (
          (((fixCode%prime_var1)==0 || (fixCode%prime_mean1)==0) && ry==0) ||
          (((fixCode%prime_var2)==0 || (fixCode%prime_mean2)==0) && ry==1) ||
          (((fixCode%prime_var3)==0 || (fixCode%prime_mean3)==0) && ry==2)
          ) continue;
        */
        invCovMat[rx][ry] = dynamic_cast<const RooRealVar*>(matrixElement.at(3*rx+ry))->getVal();
      }
    }

    Double_t offxdiff = ((invCovMat[0][1]+invCovMat[1][0])*dy + (invCovMat[0][2]+invCovMat[2][0])*dz)/2.;
    Double_t diffothers = invCovMat[1][1]*pow(dy, 2) + (invCovMat[1][2]+invCovMat[2][1])*dy*dz + invCovMat[2][2]*pow(dz, 2);
    const Double_t tolerance = epsilon*0.01;
    if (fabs(invCovMat[0][0])>tolerance){
      Double_t overallfac
        = exp(-0.5*diffothers) / sqrt(2.*pow(invCovMat[0][0], 5));
      Double_t term1
        = sqrt(2.*invCovMat[0][0])
        * (
        exp(-0.5*dxmax*(2.*offxdiff + invCovMat[0][0]*dxmax))*(offxdiff - invCovMat[0][0]*(xb + xmax))
        - exp(-0.5*dxmin*(2.*offxdiff + invCovMat[0][0]*dxmin))*(offxdiff - invCovMat[0][0]*(xb + xmin))
        );
      Double_t term2
        = exp(0.5*pow(offxdiff, 2)/invCovMat[0][0])
        * sqrt(TMath::Pi())
        * (invCovMat[0][0] + pow(offxdiff, 2) - 2.*invCovMat[0][0]*offxdiff*xb + pow(invCovMat[0][0]*xb, 2))
        *(
        TMath::Erf((offxdiff + invCovMat[0][0]*dxmax)/sqrt(2.*invCovMat[0][0]))
        - TMath::Erf((offxdiff + invCovMat[0][0]*dxmin)/sqrt(2.*invCovMat[0][0]))
        );

      value = (term1+term2)*overallfac;
    }
    else if (fabs(offxdiff)>tolerance){
      value =
        exp(
        -0.5*diffothers
        )*(
        exp(-offxdiff*dxmin)*(2.+2.*offxdiff*xmin+pow(offxdiff*xmin, 2))
        - exp(-offxdiff*dxmax)*(2.+2.*offxdiff*xmax+pow(offxdiff*xmax, 2))
        )/pow(offxdiff, 3);
    }
    else{
      value =
        exp(
        -0.5*diffothers
        )*(pow(xmax, 3)-pow(xmin, 3))/3.;
    }

    // intVarpr==0, so y==lambda
    if (fabs(cos(y))<sqrt(epsilon)) value /= epsilon;
    else value /= pow(cos(y), 2);
  }

  if (verbosity==RooGaussianMomConstraint::kVerbose)
    cout << "RooGaussianMomConstraint::computeCaseRhoLambdaPhi code " << code << " with result " << value << endl;

  if (!(value==value) || value<epsilon) value = epsilon;

  if (verbosity==RooGaussianMomConstraint::kVerbose)
    cout << "End RooGaussianMomConstraint::computeCaseRhoLambdaPhi with result " << value << endl;

  return value;
}



