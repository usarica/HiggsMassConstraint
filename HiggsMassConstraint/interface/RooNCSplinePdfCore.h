#ifndef ROONCSPLINEPDFCORE
#define ROONCSPLINEPDFCORE  

#include <vector>
#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooRealVar.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include "TH3F.h"
#include "TH1.h"
#include "RooDataHist.h"
#include "RooHistFunc.h"

class RooNCSplinePdfCore : public RooAbsPdf{
public:
  RooNCSplinePdfCore();
  RooNCSplinePdfCore(
    const char* name,
    const char* title
    );
  RooNCSplinePdfCore(const RooNCSplinePdfCore& other, const char* name=0);
  virtual TObject* clone(const char* newname)const=0;
  inline virtual ~RooNCSplinePdfCore(){}

protected:
  virtual Int_t getWhichBin(const Double_t& val, const Int_t whichDirection)const = 0;
  virtual void getKappa(std::vector<Double_t>& kappas, const Int_t whichDirection)const = 0;
  virtual Double_t getTVar(const std::vector<Double_t>& kappas, const Double_t& val, const Int_t& bin, const Int_t whichDirection)const = 0;

  virtual Double_t interpolateFcn(Int_t code, const char* rangeName=0)const = 0;
  virtual Double_t evaluate()const = 0;
  virtual Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0)const = 0;
  virtual Double_t analyticalIntegral(Int_t code, const char* rangeName=0)const = 0;

  virtual void getBArray(const std::vector<Double_t>& kappas, const std::vector<Double_t>& fcnList, std::vector<Double_t>& BArray)const;
  virtual void getAArray(const std::vector<Double_t>& kappas, std::vector<std::vector<Double_t>>& AArray)const;
  virtual std::vector<std::vector<Double_t>> getCoefficientsAlongDirection(const std::vector<Double_t>& kappas, const TMatrixD& Ainv, const std::vector<Double_t>& fcnList, const Int_t pickBin)const;
  virtual std::vector<Double_t> getCoefficients(const TVectorD& S, const std::vector<Double_t>& kappas, const std::vector<Double_t>& fcnList, const Int_t& bin)const;

  virtual Double_t evalSplineSegment(const std::vector<Double_t>& coefs, Double_t tup, Double_t tdn=0, Bool_t doIntegrate=false)const;

protected:
  ClassDef(RooNCSplinePdfCore, 0)
};
 
#endif
