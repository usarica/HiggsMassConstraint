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
typedef triplet<RooNCSplinePdfCore::T> doubleTriplet_t;

class NCSplinePdfFactory_2D{
protected:
  TString appendName;

  RooAbsReal* XVar;
  RooAbsReal* YVar;
  RooNCSplinePdf_2D_fast* PDF;

  const std::vector<doubleTriplet_t> getPoints(const std::vector<RooNCSplinePdfCore::T>& XList, const std::vector<RooNCSplinePdfCore::T>& YList, const std::vector<RooNCSplinePdfCore::T>& FcnList);

  void destroyPDF();
  void initPDF(const std::vector<doubleTriplet_t>& pList);

  void addUnique(std::vector<RooNCSplinePdfCore::T>& list, RooNCSplinePdfCore::T val);

public:
  NCSplinePdfFactory_2D(RooAbsReal& XVar_, RooAbsReal& YVar_, TString appendName_="");
  ~NCSplinePdfFactory_2D();

  RooNCSplinePdf_2D_fast* getPDF();

  void setPoints(const std::vector<doubleTriplet_t>& pList){ initPDF(pList); }
  template<typename inType> void setPoints(const std::vector<inType>& XList, const std::vector<inType>& YList, const std::vector<inType>& FcnList){
    std::vector<RooNCSplinePdfCore::T> transXList;
    std::vector<RooNCSplinePdfCore::T> transYList;
    std::vector<RooNCSplinePdfCore::T> transFcnList;
    for (unsigned int ip=0; ip<XList.size(); ip++) transXList.push_back((RooNCSplinePdfCore::T)XList.at(ip));
    for (unsigned int ip=0; ip<YList.size(); ip++) transYList.push_back((RooNCSplinePdfCore::T)YList.at(ip));
    for (unsigned int ip=0; ip<FcnList.size(); ip++) transFcnList.push_back((RooNCSplinePdfCore::T)FcnList.at(ip));
    const std::vector<doubleTriplet_t> pList = getPoints(transXList, transYList, transFcnList);
    setPoints(pList);
  }

};

template void NCSplinePdfFactory_2D::setPoints<Float_t>(const std::vector<Float_t>& XList, const std::vector<Float_t>& YList, const std::vector<Float_t>& FcnList);
template void NCSplinePdfFactory_2D::setPoints<Double_t>(const std::vector<Double_t>& XList, const std::vector<Double_t>& YList, const std::vector<Double_t>& FcnList);

#endif



