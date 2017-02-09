#ifndef SLICEPDFFACTORY_H
#define SLICEPDFFACTORY_H

#include <vector>
#include <utility>
#include <algorithm>
#include "RooConstVar.h"
#include "RooSlicePdf.h"
#include "TGraph.h"

class SlicePdfFactory{
protected:
  TString appendName;

  RooAbsReal* sliceVar;
  const RooArgList depList;
  const RooArgList funcList;
  const RooArgList coefList;
  RooSlicePdf* PDF;

  std::vector<RooConstVar*> points;

public:
  SlicePdfFactory(RooAbsReal* sliceVar_, const RooArgList& depList_, const RooArgList& funcList_, const RooArgList& coefList_, TString appendName_="");
  ~SlicePdfFactory();

  RooSlicePdf* getPDF();
  void setXCoordinates(std::vector<Double_t> coordx);

protected:
  void destroyX();
  void initX(std::vector<Double_t> coordx);

  void destroyPDF();
  void initPDF();

};


#endif



