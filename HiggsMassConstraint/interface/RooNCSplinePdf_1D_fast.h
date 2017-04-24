#ifndef ROONCSPLINEPDF_FAST
#define ROONCSPLINEPDF_FAST

#include <vector>
#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooRealVar.h"
#include "RooAbsReal.h"
#include "RooNCSplinePdfCore.h"


class RooNCSplinePdf_1D_fast : public RooNCSplinePdfCore{
protected:
  std::vector<T> FcnList; // List of function values

  std::vector<T> kappaX;
  std::vector<std::vector<T>> coefficients;

public:
  RooNCSplinePdf_1D_fast();
  RooNCSplinePdf_1D_fast(
    const char* name,
    const char* title,
    RooAbsReal& inXVar,
    std::vector<T>& inXList,
    std::vector<T>& inFcnList
    );
  RooNCSplinePdf_1D_fast(const RooNCSplinePdf_1D_fast& other, const char* name=0);
	virtual TObject* clone(const char* newname)const { return new RooNCSplinePdf_1D_fast(*this, newname); }
	inline virtual ~RooNCSplinePdf_1D_fast(){}

protected:
  virtual void emptyFcnList(){ std::vector<T> tmp; FcnList.swap(tmp); }

  virtual Int_t getWhichBin(const T& val, const Int_t whichDirection)const;
  virtual void getKappas(std::vector<T>& kappas, const Int_t whichDirection)const;
  virtual T getTVar(const std::vector<T>& kappas, const T& val, const Int_t& bin, const Int_t whichDirection)const;

  virtual T interpolateFcn(Int_t code, const char* rangeName=0)const;

  virtual Double_t evaluate()const;
  virtual Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0)const;
  virtual Double_t analyticalIntegral(Int_t code, const char* rangeName=0)const;

private:

  ClassDef(RooNCSplinePdf_1D_fast, 1)

};
 
#endif
