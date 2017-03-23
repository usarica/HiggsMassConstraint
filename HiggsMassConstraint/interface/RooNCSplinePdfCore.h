#ifndef ROONCSPLINEPDFCORE
#define ROONCSPLINEPDFCORE  

#include <vector>
#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooRealProxy.h"
#include "RooConstVar.h"
#include "RooArgList.h"

class RooNCSplinePdfCore : public RooAbsPdf{
public:
  typedef Double_t T;

  enum VerbosityLevel{
    kSilent,
    kVerbose
  };

  RooNCSplinePdfCore();
  RooNCSplinePdfCore(
    const char* name,
    const char* title
    );
  RooNCSplinePdfCore(
    const char* name,
    const char* title,
    RooAbsReal& inXVar,
    std::vector<T>& inXList
    );
  RooNCSplinePdfCore(const RooNCSplinePdfCore& other, const char* name=0);
  virtual TObject* clone(const char* newname)const=0;
  inline virtual ~RooNCSplinePdfCore(){}

  virtual void setVerbosity(VerbosityLevel flag);

protected:
  VerbosityLevel verbosity;

  RooRealProxy theXVar;
  std::vector<T> XList;
  
  unsigned int npointsX()const{ return XList.size()-1; }

  virtual Int_t getWhichBin(const T& val, const Int_t whichDirection)const = 0;
  virtual void getKappa(std::vector<T>& kappas, const Int_t whichDirection)const = 0;
  virtual T getTVar(const std::vector<T>& kappas, const T& val, const Int_t& bin, const Int_t whichDirection)const = 0;

  virtual Double_t interpolateFcn(Int_t code, const char* rangeName=0)const = 0;
  virtual Double_t evaluate()const = 0;
  virtual Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0)const = 0;
  virtual Double_t analyticalIntegral(Int_t code, const char* rangeName=0)const = 0;

  virtual void getBArray(const std::vector<T>& kappas, const std::vector<T>& fcnList, std::vector<T>& BArray)const;
  virtual void getAArray(const std::vector<T>& kappas, std::vector<std::vector<T>>& AArray)const;
  virtual std::vector<std::vector<T>> getCoefficientsAlongDirection(const std::vector<T>& kappas, const TMatrixD& Ainv, const std::vector<T>& fcnList, const Int_t pickBin)const;
  virtual std::vector<T> getCoefficients(const TVectorD& S, const std::vector<T>& kappas, const std::vector<T>& fcnList, const Int_t& bin)const;

  virtual Double_t evalSplineSegment(const std::vector<T>& coefs, const T kappa, T tup, T tdn=0, Bool_t doIntegrate=false)const;

protected:
  ClassDef(RooNCSplinePdfCore, 0)
};
 
#endif
