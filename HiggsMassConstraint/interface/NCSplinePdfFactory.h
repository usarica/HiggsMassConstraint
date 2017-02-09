#ifndef NCSPLINEPDFFACTORY_H
#define NCSPLINEPDFFACTORY_H

#include <vector>
#include <utility>
#include <algorithm>
#include "RooConstVar.h"
#include "RooNCSplinePdf.h"
#include "RooNCSplinePdf_fast.h"
#include "TGraph.h"

class NCSplinePdfFactory{
protected:
  TString appendName;

  RooAbsReal* splineVar;
  RooNCSplinePdf_fast* PDF;

  std::vector<std::pair<RooConstVar*, RooConstVar*>> points;

public:
  NCSplinePdfFactory(RooAbsReal* splineVar_, TString appendName_="");
  ~NCSplinePdfFactory();

  RooNCSplinePdf_fast* getPDF();
  void setGraph(TGraph* tg);

protected:
  void destroyXY();
  void initXY(TGraph* tg);

  void destroyPDF();
  void initPDF();

};


#endif



