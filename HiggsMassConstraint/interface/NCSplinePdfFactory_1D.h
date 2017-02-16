#ifndef NCSPLINEPDFFACTORY_1D
#define NCSPLINEPDFFACTORY_1D

#include <vector>
#include <utility>
#include <algorithm>
#include "RooConstVar.h"
#include "RooNCSplinePdf_1D_fast.h"
#include "TGraph.h"

class NCSplinePdfFactory_1D{
protected:
  TString appendName;

  RooAbsReal* splineVar;
  RooNCSplinePdf_1D_fast* PDF;

  std::vector<std::pair<RooConstVar*, RooConstVar*>> points;

public:
  NCSplinePdfFactory_1D(RooAbsReal* splineVar_, TString appendName_="");
  ~NCSplinePdfFactory_1D();

  RooNCSplinePdf_1D_fast* getPDF();
  void setPoints(TGraph* tg);
  void setPoints(const std::vector<Double_t>& XList, const std::vector<Double_t>& FcnList);

protected:
  void destroyPoints();
  void initPoints(const std::vector<Double_t>& XList, const std::vector<Double_t>& FcnList);

  void destroyPDF();
  void initPDF();

};


#endif



