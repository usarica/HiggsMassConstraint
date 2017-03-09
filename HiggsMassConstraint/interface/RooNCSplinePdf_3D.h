#ifndef ROONCSPLINEPDF_3D
#define ROONCSPLINEPDF_3D

#include <vector>
#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooRealVar.h"
#include "RooAbsReal.h"
#include "RooNCSplinePdfCore.h"

class RooNCSplinePdf_3D : public RooNCSplinePdfCore{
protected:
  RooRealProxy theYVar;
  RooRealProxy theZVar;
  RooListProxy YList;
  RooListProxy ZList;
  Int_t npointsY;
  Int_t npointsZ;

  std::vector<std::vector<RooListProxy>> FcnList;

public:
  RooNCSplinePdf_3D();
  RooNCSplinePdf_3D(
    const char* name,
    const char* title,
    RooAbsReal* inXVar,
    RooAbsReal* inYVar,
    RooAbsReal* inZVar,
    const RooArgList* inXList,
    const RooArgList* inYList,
    const RooArgList* inZList,
    std::vector<std::vector<const RooArgList*>>& inFcnList,
    bool inUseConst=false
    );
  RooNCSplinePdf_3D(const RooNCSplinePdf_3D& other, const char* name=0);
	virtual TObject* clone(const char* newname)const { return new RooNCSplinePdf_3D(*this, newname); }
	inline virtual ~RooNCSplinePdf_3D(){}

protected:
  virtual Int_t getWhichBin(const Double_t& val, const Int_t whichDirection)const;
  virtual Double_t getTVar(const std::vector<Double_t>& kappas, const Double_t& val, const Int_t& bin, const Int_t whichDirection)const;
  virtual void getKappa(std::vector<Double_t>& kappas, const Int_t whichDirection)const;

  virtual std::vector<std::vector<Double_t>> getCoefficientsPerYPerZ(
    const std::vector<Double_t>& kappaX, const TMatrixD& xAinv,
    const Int_t& ybin, const Int_t& zbin,
    const Int_t xbin
    )const; // xbin can be -1, which means push all of them

  virtual Double_t interpolateFcn(Int_t code, const char* rangeName=0)const;
  virtual Double_t evaluate()const;
  virtual Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const;
  virtual Double_t analyticalIntegral(Int_t code, const char* rangeName=0)const;

private:
  ClassDef(RooNCSplinePdf_3D, 1)
};
 
#endif
