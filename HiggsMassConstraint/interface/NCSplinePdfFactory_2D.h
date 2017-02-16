#ifndef NCSPLINEPDFFACTORY_2D
#define NCSPLINEPDFFACTORY_2D

#include <vector>
#include <utility>
#include <algorithm>
#include "RooConstVar.h"
#include "RooNCSplinePdf_2D_fast.h"

template<typename T> struct triplet{
  T value[3];
  triplet(T i1, T i2, T i3){
    value[0]=i1;
    value[1]=i2;
    value[2]=i3;
  }
  triplet(T i1){
    value[0]=i1;
    value[1]=i1;
    value[2]=i1;
  }
  triplet(){}
  T& operator[](std::size_t ipos){ return value[ipos]; } // Return by reference
  const T& operator[](std::size_t ipos)const{ return value[ipos]; } // Return by const reference
};
typedef triplet<Double_t> doubleTriplet_t;

class NCSplinePdfFactory_2D{
protected:
  TString appendName;

  RooAbsReal* XVar;
  RooAbsReal* YVar;
  RooNCSplinePdf_2D_fast* PDF;

  std::vector<RooConstVar*> Xcoord;
  std::vector<RooConstVar*> Ycoord;
  std::vector<RooConstVar*> FcnVal;

public:
  NCSplinePdfFactory_2D(RooAbsReal* XVar_, RooAbsReal* YVar_, TString appendName_="");
  ~NCSplinePdfFactory_2D();

  RooNCSplinePdf_2D_fast* getPDF();

  void setPoints(const std::vector<doubleTriplet_t>& pList);

protected:
  void destroyPoints();
  void initPoints(const std::vector<doubleTriplet_t>& pList);

  void destroyPDF();
  void initPDF();

  void addUnique(std::vector<Double_t>& list, Double_t val);

};


#endif



