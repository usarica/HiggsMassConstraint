#ifndef ROONCSPLINEPDF_FAST
#define ROONCSPLINEPDF_FAST

#include <vector>
#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooRealVar.h"
#include "RooAbsReal.h"
#include "RooNCSplinePdf_1D.h"


class RooNCSplinePdf_1D_fast : public RooNCSplinePdf_1D{
protected:
  std::vector<Double_t> kappaX;
  std::vector<std::vector<Double_t>> coefficients;

public:
  RooNCSplinePdf_1D_fast();
  RooNCSplinePdf_1D_fast(
    const char* name,
    const char* title,
    RooAbsReal* inXVar,
    const RooArgList* inXList,
    const RooArgList* inFcnList
    );
  RooNCSplinePdf_1D_fast(const RooNCSplinePdf_1D_fast& other, const char* name=0);
	virtual TObject* clone(const char* newname)const { return new RooNCSplinePdf_1D_fast(*this, newname); }
	inline virtual ~RooNCSplinePdf_1D_fast(){}

protected:
  virtual Double_t interpolateFcn(Int_t code, const char* rangeName=0)const;

  virtual Double_t evaluate()const;
  virtual Double_t analyticalIntegral(Int_t code, const char* rangeName=0)const;

private:
  ClassDef(RooNCSplinePdf_1D_fast, 1)
};
 
#endif
