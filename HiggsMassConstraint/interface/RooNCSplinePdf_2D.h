#ifndef ROONCSPLINEPDF_2D
#define ROONCSPLINEPDF_2D

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

class RooNCSplinePdf_2D : public RooAbsPdf{
protected:

  RooRealProxy theXVar;
  RooRealProxy theYVar;
  RooListProxy XList; // List of X values
  RooListProxy YList; // List of Y values
  RooListProxy ZList; // List of Z values
  int npointsX;
  int npointsY;

public:
  RooNCSplinePdf_2D();
  RooNCSplinePdf_2D(
    const char* name,
    const char* title,
    RooAbsReal& inXVar,
    RooAbsReal& inYVar,
    const RooArgList& inXList, // X and Y define the grid
    const RooArgList& inYList,
    const RooArgList& inZList // Z has dimension dim(X)*dim(Y) with Z[i][j] corresponding to X[i], Y[j]
    );
  RooNCSplinePdf_2D(const RooNCSplinePdf_2D& other, const char* name=0);
	virtual TObject* clone(const char* newname)const { return new RooNCSplinePdf_2D(*this, newname); }
	inline virtual ~RooNCSplinePdf_2D(){}

protected:
  virtual Int_t getWhichBin(const Double_t& val, const Int_t whichDirection)const;
  virtual Double_t getTVar(const std::vector<Double_t>& kappas, const Double_t& val, const Int_t& bin, const Int_t whichDirection)const;

  virtual void getKappa(std::vector<Double_t>& kappas, const Int_t whichDirection)const;

  virtual void getBArray(const std::vector<Double_t>& kappas, const std::vector<Double_t>& fcnList, std::vector<Double_t>& BArray)const;
  virtual void getAArray(const std::vector<Double_t>& kappas, std::vector<std::vector<Double_t>>& AArray)const;
  virtual std::vector<std::vector<Double_t>> getCoefficientsAlongDirection(const std::vector<Double_t>& kappas, const TMatrixD& Ainv, const std::vector<Double_t>& fcnList, const Int_t pickBin)const;

  virtual std::vector<Double_t> getCoefficients(const TVectorD& S, const std::vector<Double_t>& kappas, const std::vector<Double_t>& fcnList, const Int_t& bin)const;

  virtual std::vector<std::vector<Double_t>> getCoefficientsPerY(const std::vector<Double_t>& kappaX, const TMatrixD& xAinv, const Int_t& ybin, const Int_t xbin)const; // xbin can be -1, which means push all of them

  virtual Double_t interpolateFcn(Int_t code, const char* rangeName=0)const;
  virtual Double_t evaluate()const;
  virtual Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const;
  virtual Double_t analyticalIntegral(Int_t code, const char* rangeName=0)const;

private:
  ClassDef(RooNCSplinePdf_2D, 1)
};
 
#endif
