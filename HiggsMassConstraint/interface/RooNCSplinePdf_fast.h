#ifndef ROONCSPLINEPDF_FAST
#define ROONCSPLINEPDF_FAST

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
#include "RooNCSplinePdf.h"


class RooNCSplinePdf_fast : public RooNCSplinePdf{
protected:
  std::vector<Double_t> kappas;
  std::vector<std::vector<Double_t>> coefficients;

public:
  RooNCSplinePdf_fast();
  RooNCSplinePdf_fast(
    const char* name,
    const char* title,
    RooAbsReal& inVar,
    const RooArgList& inXList,
    const RooArgList& inYList
    );
  RooNCSplinePdf_fast(const RooNCSplinePdf_fast& other, const char* name=0);
	virtual TObject* clone(const char* newname)const { return new RooNCSplinePdf_fast(*this, newname); }
	inline virtual ~RooNCSplinePdf_fast(){}

protected:
  virtual Double_t interpolateFcn(Int_t code, const char* rangeName=0)const;
  virtual Double_t evaluate()const;
  //virtual Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const;
  virtual Double_t analyticalIntegral(Int_t code, const char* rangeName=0)const;

private:
  ClassDef(RooNCSplinePdf_fast, 1)
};
 
#endif
