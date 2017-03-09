#ifndef ROONCSPLINEPDF_1D
#define ROONCSPLINEPDF_1D

#include <vector>
#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooRealVar.h"
#include "RooAbsReal.h"
#include "RooNCSplinePdfCore.h"

class RooNCSplinePdf_1D : public RooNCSplinePdfCore{
protected:
  RooListProxy FcnList; // List of function values

public:
  RooNCSplinePdf_1D();
  RooNCSplinePdf_1D(
    const char* name,
    const char* title,
    RooAbsReal* inXVar,
    const RooArgList* inXList,
    const RooArgList* inFcnList,
    bool inUseConst=false
    );
  RooNCSplinePdf_1D(const RooNCSplinePdf_1D& other, const char* name=0);
	virtual TObject* clone(const char* newname)const { return new RooNCSplinePdf_1D(*this, newname); }
	inline virtual ~RooNCSplinePdf_1D(){}

protected:
  virtual Int_t getWhichBin(const Double_t& val, const Int_t whichDirection)const;
  virtual void getKappa(std::vector<Double_t>& kappas, const Int_t whichDirection)const;
  virtual Double_t getTVar(const std::vector<Double_t>& kappas, const Double_t& val, const Int_t& bin, const Int_t whichDirection)const;

  virtual Double_t interpolateFcn(Int_t code, const char* rangeName=0)const;
  virtual Double_t evaluate()const;
  virtual Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0)const;
  virtual Double_t analyticalIntegral(Int_t code, const char* rangeName=0)const;

private:
  ClassDef(RooNCSplinePdf_1D, 1)
};
 
#endif
