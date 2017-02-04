#ifndef ROONCSPLINEPDF
#define ROONCSPLINEPDF

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


class RooNCSplinePdf : public RooAbsPdf {
protected:

	RooRealProxy theVar;
  RooListProxy XList; // List of X values
  RooListProxy YList; // List of Y values
  Int_t npoints;

public:
  RooNCSplinePdf();
  RooNCSplinePdf(
    const char *name,
    const char *title,
    RooAbsReal& inVar,
    const RooArgList& inXList,
    const RooArgList& inYList
    );
  RooNCSplinePdf(const RooNCSplinePdf& other, const char* name=0);
	virtual TObject* clone(const char* newname)const { return new RooNCSplinePdf(*this, newname); }
	inline virtual ~RooNCSplinePdf(){}

protected:
  virtual void getKappa(std::vector<Double_t>& kappas)const;
  virtual void getBArray(const std::vector<Double_t>& kappas, std::vector<Double_t>& BArray)const;
  virtual void getAArray(const std::vector<Double_t>& kappas, std::vector<std::vector<Double_t>>& AArray)const;
  virtual Int_t getWhichBin(const Double_t& xval)const;
  virtual Double_t getTVar(const std::vector<Double_t>& kappas, const Double_t& xval, const Int_t& bin)const;
  virtual std::vector<Double_t> getCoefficients(const TVectorD& S, const std::vector<Double_t>& kappas, const Int_t& bin)const;

  virtual Double_t interpolateFcn(Bool_t doIntegrate, const char* rangeName=0)const;
  virtual Double_t evaluate()const;
  virtual Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const;
  virtual Double_t analyticalIntegral(Int_t code, const char* rangeName=0)const;

};
 
#endif
