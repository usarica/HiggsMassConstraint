#ifndef NCSPLINEPDFFACTORY_1D
#define NCSPLINEPDFFACTORY_1D

#include <vector>
#include <utility>
#include <algorithm>
#include "RooNCSplinePdf_1D_fast.h"
#include "TGraph.h"

class NCSplinePdfFactory_1D{
protected:
  TString appendName;

  RooAbsReal* splineVar;
  RooNCSplinePdf_1D_fast* PDF;

  const std::vector<std::pair<RooNCSplinePdfCore::T, RooNCSplinePdfCore::T>> getPoints(const std::vector<RooNCSplinePdfCore::T>& XList, const std::vector<RooNCSplinePdfCore::T>& FcnList);

  void destroyPDF();
  void initPDF(const std::vector<std::pair<RooNCSplinePdfCore::T, RooNCSplinePdfCore::T>>& pList);

public:
  NCSplinePdfFactory_1D(RooAbsReal& splineVar_, TString appendName_="");
  ~NCSplinePdfFactory_1D();

  RooNCSplinePdf_1D_fast* getPDF();
  void setPoints(TGraph* tg);
  void setPoints(const std::vector<std::pair<RooNCSplinePdfCore::T, RooNCSplinePdfCore::T>>& pList){ initPDF(pList); }
  template<typename inType> void setPoints(const std::vector<inType>& XList, const std::vector<inType>& FcnList){
    std::vector<RooNCSplinePdfCore::T> transXList;
    std::vector<RooNCSplinePdfCore::T> transFcnList;
    for (unsigned int ip=0; ip<XList.size(); ip++) transXList.push_back((RooNCSplinePdfCore::T)XList.at(ip));
    for (unsigned int ip=0; ip<FcnList.size(); ip++) transFcnList.push_back((RooNCSplinePdfCore::T)FcnList.at(ip));
    const std::vector<std::pair<RooNCSplinePdfCore::T, RooNCSplinePdfCore::T>> pList = getPoints(transXList, transFcnList);
    initPDF(pList);
  }

};

template void NCSplinePdfFactory_1D::setPoints<Float_t>(const std::vector<Float_t>& XList, const std::vector<Float_t>& FcnList);
template void NCSplinePdfFactory_1D::setPoints<Double_t>(const std::vector<Double_t>& XList, const std::vector<Double_t>& FcnList);

#endif



