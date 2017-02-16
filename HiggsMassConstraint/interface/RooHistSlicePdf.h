#ifndef ROOHISTSLICEPDF
#define ROOHISTSLICEPDF

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


class RooHistSlicePdf : public RooAbsPdf {
protected:

	RooRealProxy sliceVar;
  RooListProxy depList; // List of pdf dependents
  RooListProxy funcList; // List of histogram pdfs
  RooListProxy coefList; // List of histogram interpolation coefficients

  Int_t nSlices;
  // Valid slice min/max does not necessarily have to be the same as sliceVar.getMin/Max()
  Double_t sliceVar_validitymin;
  Double_t sliceVar_validitymax;
  std::vector<Double_t> sliceIntegral;
  Int_t smoothAlgo; // 0: Quadratic interpolation null at zero and continuous at boundaries but not smooth at boundaries (default); 1: Quadratic interpolation that is everywhere differentiable but not null at zero
  Double_t smoothRegion; // Usually 1.

  virtual Double_t evaluate() const;

public:
  RooHistSlicePdf();
	RooHistSlicePdf(
		const char* name,
		const char* title,
		RooAbsReal& inSliceVar,
    const RooArgList& inDepList,
    const RooArgList& inFuncList,
    const RooArgList& inCoefList,

		Double_t inSliceVar_validitymin,
		Double_t inSliceVar_validitymax,

    Int_t inSmoothAlgo=0,
    Double_t inSmoothRegion=1.
		);

	RooHistSlicePdf(const RooHistSlicePdf& other, const char* name = 0);
	virtual TObject* clone(const char* newname) const { return new RooHistSlicePdf(*this, newname); }
	inline virtual ~RooHistSlicePdf(){}

	virtual Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName = 0) const;
  virtual Double_t analyticalIntegral(Int_t code, const char* rangeName = 0) const;

protected:
  virtual Int_t findNeighborBins() const; // Returns index_low for sliceVar pdfs
  virtual Double_t interpolateBin(Bool_t doIntegrate) const;
  virtual Double_t interpolateVariation(Double_t theta, Double_t valueCenter, Double_t valueHigh, Double_t valueLow) const;

private:
  ClassDef(RooHistSlicePdf, 1)
};
 
#endif
