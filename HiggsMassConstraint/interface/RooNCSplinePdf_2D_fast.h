#ifndef ROONCSPLINEPDF_2D_FAST
#define ROONCSPLINEPDF_2D_FAST

#include <vector>
#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooRealVar.h"
#include "RooAbsReal.h"
#include "RooNCSplinePdf_2D.h"

class RooNCSplinePdf_2D_fast : public RooNCSplinePdf_2D{
protected:
  std::vector<Double_t> kappaX;
  std::vector<Double_t> kappaY;
  std::vector<std::vector<std::vector<std::vector<Double_t>>>> coefficients; // [ix][A_x,B_x,C_x,D_x][iy][A_x_y,B_x_y,C_x_y,D_x_y]

public:
  RooNCSplinePdf_2D_fast();
  RooNCSplinePdf_2D_fast(
    const char* name,
    const char* title,
    RooAbsReal* inXVar,
    RooAbsReal* inYVar,
    const RooArgList* inXList,
    const RooArgList* inYList,
    std::vector<const RooArgList*>& inFcnList
    );
  RooNCSplinePdf_2D_fast(const RooNCSplinePdf_2D_fast& other, const char* name=0);
	virtual TObject* clone(const char* newname)const { return new RooNCSplinePdf_2D_fast(*this, newname); }
	inline virtual ~RooNCSplinePdf_2D_fast(){}

protected:
  virtual Double_t interpolateFcn(Int_t code, const char* rangeName=0)const;
  virtual Double_t evaluate()const;
  virtual Double_t analyticalIntegral(Int_t code, const char* rangeName=0)const;

private:
  ClassDef(RooNCSplinePdf_2D_fast, 1)
};
 
#endif
