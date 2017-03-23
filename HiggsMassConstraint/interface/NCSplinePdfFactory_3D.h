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
typedef quadruplet<RooNCSplinePdfCore::T> doubleQuadruplet_t;

class NCSplinePdfFactory_3D{
protected:
  TString appendName;

  RooAbsReal* XVar;
  RooAbsReal* YVar;
  RooAbsReal* ZVar;
  RooNCSplinePdf_3D_fast* PDF;

  const std::vector<doubleQuadruplet_t> getPoints(const std::vector<RooNCSplinePdfCore::T>& XList, const std::vector<RooNCSplinePdfCore::T>& YList, const std::vector<RooNCSplinePdfCore::T>& ZList, const std::vector<RooNCSplinePdfCore::T>& FcnList);

  void destroyPDF();
  void initPDF(const std::vector<doubleQuadruplet_t>& pList);

  void addUnique(std::vector<RooNCSplinePdfCore::T>& list, RooNCSplinePdfCore::T val);

public:
  NCSplinePdfFactory_3D(RooAbsReal& XVar_, RooAbsReal& YVar_, RooAbsReal& ZVar_, TString appendName_="");
  ~NCSplinePdfFactory_3D();

  RooNCSplinePdf_3D_fast* getPDF();

  void setPoints(const std::vector<doubleQuadruplet_t>& pList){ initPDF(pList); }
  template<typename inType> void setPoints(const std::vector<inType>& XList, const std::vector<inType>& YList, const std::vector<inType>& ZList, const std::vector<inType>& FcnList){
    std::vector<RooNCSplinePdfCore::T> transXList;
    std::vector<RooNCSplinePdfCore::T> transYList;
    std::vector<RooNCSplinePdfCore::T> transZList;
    std::vector<RooNCSplinePdfCore::T> transFcnList;
    for (unsigned int ip=0; ip<XList.size(); ip++) transXList.push_back((RooNCSplinePdfCore::T)XList.at(ip));
    for (unsigned int ip=0; ip<YList.size(); ip++) transYList.push_back((RooNCSplinePdfCore::T)YList.at(ip));
    for (unsigned int ip=0; ip<ZList.size(); ip++) transZList.push_back((RooNCSplinePdfCore::T)ZList.at(ip));
    for (unsigned int ip=0; ip<FcnList.size(); ip++) transFcnList.push_back((RooNCSplinePdfCore::T)FcnList.at(ip));
    const std::vector<doubleQuadruplet_t> pList = getPoints(transXList, transYList, transZList, transFcnList);
    setPoints(pList);
  }

};

template void NCSplinePdfFactory_3D::setPoints<Float_t>(const std::vector<Float_t>& XList, const std::vector<Float_t>& YList, const std::vector<Float_t>& ZList, const std::vector<Float_t>& FcnList);
template void NCSplinePdfFactory_3D::setPoints<Double_t>(const std::vector<Double_t>& XList, const std::vector<Double_t>& YList, const std::vector<Double_t>& ZList, const std::vector<Double_t>& FcnList);


#endif



