#ifndef ROOSLICEPDF
#define ROOSLICEPDF

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


class RooSlicePdf : public RooAbsPdf {
protected:

	RooRealProxy sliceVar;
  RooListProxy depList; // List of pdf dependents
  RooListProxy sliceList; // List of slice coordinates
  RooListProxy funcList; // List of slice pdfs
  RooListProxy coefList; // List of slice interpolation coefficients

  Int_t nSlices;

  std::vector<Double_t> sliceIntegral;
  Int_t smoothAlgo; // 0: Quadratic interpolation null at zero and continuous at boundaries but not smooth at boundaries (default); 1: Quadratic interpolation that is everywhere differentiable but not null at zero
  Double_t smoothRegion; // Usually 1.

  virtual Double_t evaluate() const;

public:
  RooSlicePdf();
	RooSlicePdf(
		const char* name,
		const char* title,
		RooAbsReal& inSliceVar,
    const RooArgList& inDepList,
    const RooArgList& inSliceList,
    const RooArgList& inFuncList,
    const RooArgList& inCoefList,

    Int_t inSmoothAlgo=0,
    Double_t inSmoothRegion=1.
		);

	RooSlicePdf(const RooSlicePdf& other, const char* name = 0);
	virtual TObject* clone(const char* newname) const { return new RooSlicePdf(*this, newname); }
	inline virtual ~RooSlicePdf(){}

	virtual Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName = 0) const;
  virtual Double_t analyticalIntegral(Int_t code, const char* rangeName = 0) const;

protected:
  virtual Int_t findNeighborBins() const; // Returns index_low for sliceVar pdfs
  virtual Double_t findBinWidth(Int_t bin) const;
  virtual Double_t interpolateBin(Int_t intCode) const;
  virtual Double_t interpolateVariation(Double_t theta, Double_t valueCenter, Double_t valueHigh, Double_t valueLow) const;

  virtual void setListProxyFromArgList(const RooArgList& argList, RooListProxy& proxyList);

private:
  ClassDef(RooSlicePdf, 1)
};
 
#endif
