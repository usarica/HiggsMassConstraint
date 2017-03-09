#ifndef NCSPLINEPDFFACTORY_3D
#define NCSPLINEPDFFACTORY_3D

#include <vector>
#include <utility>
#include <algorithm>
#include "RooConstVar.h"
#include "RooNCSplinePdf_3D_fast.h"

template<typename T> struct quadruplet{
  T value[4];
  quadruplet(T i1, T i2, T i3, T i4){
    value[0]=i1;
    value[1]=i2;
    value[2]=i3;
    value[3]=i4;
  }
  quadruplet(T i1){ for (unsigned int idim=0; idim<4; idim++) value[idim] = i1; }
  quadruplet(){}
  T& operator[](std::size_t ipos){ return value[ipos]; } // Return by reference
  const T& operator[](std::size_t ipos)const{ return value[ipos]; } // Return by const reference
};
typedef quadruplet<Double_t> doubleQuadruplet_t;

class NCSplinePdfFactory_3D{
protected:
  TString appendName;

  RooAbsReal* XVar;
  RooAbsReal* YVar;
  RooAbsReal* ZVar;
  RooNCSplinePdf_3D_fast* PDF;

  std::vector<RooConstVar*> Xcoord;
  std::vector<RooConstVar*> Ycoord;
  std::vector<RooConstVar*> Zcoord;
  std::vector<std::vector<std::vector<RooConstVar*>>> FcnVal;

public:
  NCSplinePdfFactory_3D(RooAbsReal* XVar_, RooAbsReal* YVar_, RooAbsReal* ZVar_, TString appendName_="");
  ~NCSplinePdfFactory_3D();

  RooNCSplinePdf_3D_fast* getPDF();

  void setPoints(const std::vector<doubleQuadruplet_t>& pList);

protected:
  void destroyPoints();
  void initPoints(const std::vector<doubleQuadruplet_t>& pList);

  void destroyPDF();
  void initPDF();

  void addUnique(std::vector<Double_t>& list, Double_t val);

};


#endif



