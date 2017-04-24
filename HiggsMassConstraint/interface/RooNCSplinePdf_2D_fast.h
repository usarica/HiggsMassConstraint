#ifndef ROONCSPLINEPDF_2D_FAST
#define ROONCSPLINEPDF_2D_FAST

#include <vector>
#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooRealVar.h"
#include "RooAbsReal.h"
#include "RooNCSplinePdfCore.h"

class RooNCSplinePdf_2D_fast : public RooNCSplinePdfCore{
protected:
  RooRealProxy theYVar;
  std::vector<T> YList;

  std::vector<std::vector<T>> FcnList;

  std::vector<T> kappaX;
  std::vector<T> kappaY;
  std::vector<std::vector<std::vector<std::vector<T>>>> coefficients; // [ix][A_x,B_x,C_x,D_x][iy][A_x_y,B_x_y,C_x_y,D_x_y]

public:
  RooNCSplinePdf_2D_fast();
  RooNCSplinePdf_2D_fast(
    const char* name,
    const char* title,
    RooAbsReal& inXVar,
    RooAbsReal& inYVar,
    std::vector<T>& inXList,
    std::vector<T>& inYList,
    std::vector<std::vector<T>>& inFcnList
    );
  RooNCSplinePdf_2D_fast(const RooNCSplinePdf_2D_fast& other, const char* name=0);
	virtual TObject* clone(const char* newname)const { return new RooNCSplinePdf_2D_fast(*this, newname); }
	inline virtual ~RooNCSplinePdf_2D_fast(){}

protected:
  virtual void emptyFcnList(){ std::vector<std::vector<T>> tmp; FcnList.swap(tmp); }

  unsigned int npointsY()const{ return YList.size(); }

  virtual Int_t getWhichBin(const T& val, const Int_t whichDirection)const;
  virtual T getTVar(const std::vector<T>& kappas, const T& val, const Int_t& bin, const Int_t whichDirection)const;
  virtual void getKappas(std::vector<T>& kappas, const Int_t whichDirection)const;

  virtual std::vector<std::vector<T>> getCoefficientsPerY(const std::vector<T>& kappaX, const TMatrix_t& xAinv, const Int_t& ybin, const Int_t xbin)const; // xbin can be -1, which means push all of them

  virtual T interpolateFcn(Int_t code, const char* rangeName=0)const;

  virtual Double_t evaluate()const;
  virtual Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const;
  virtual Double_t analyticalIntegral(Int_t code, const char* rangeName=0)const;

private:

  ClassDef(RooNCSplinePdf_2D_fast, 1)

};
 
#endif
