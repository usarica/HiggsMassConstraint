#include "NCSplinePdfFactory.h"

using namespace std;


NCSplinePdfFactory::NCSplinePdfFactory(RooAbsReal* splineVar_, TString appendName_) :
appendName(appendName_),
splineVar(splineVar_),
PDF(0)
{}
NCSplinePdfFactory::~NCSplinePdfFactory(){
  destroyPDF();
  destroyXY();
}
void NCSplinePdfFactory::destroyXY(){
  for (unsigned int ip=0; ip<points.size(); ip++){
    delete points.at(ip).first; delete points.at(ip).second;
  }
  points.clear();
}
void NCSplinePdfFactory::setGraph(TGraph* tg){
  initXY(tg);
  initPDF();
}
void NCSplinePdfFactory::initXY(TGraph* tg){
  // Flush the points array
  destroyXY();

  double* xx = tg->GetX();
  double* yy = tg->GetY();
  int n = tg->GetN();
  for (int ip=0; ip<n; ip++){
    TString xname = Form("point_%i_x", ip);
    TString yname = Form("point_%i_y", ip);
    if (appendName!=""){
      xname = Form("%s_%s", xname.Data(), appendName.Data());
      yname = Form("%s_%s", yname.Data(), appendName.Data());
    }
    RooConstVar* xvar = new RooConstVar(xname, xname, xx[ip]);
    RooConstVar* yvar = new RooConstVar(yname, yname, yy[ip]);

    points.push_back(pair<RooConstVar*, RooConstVar*>(xvar, yvar));
  }
}

void NCSplinePdfFactory::destroyPDF(){ delete PDF; PDF=0; }
void NCSplinePdfFactory::initPDF(){
  destroyPDF();

  RooArgList XList, YList;
  for (unsigned int ip=0; ip<points.size(); ip++){
    XList.add(*(points.at(ip).first));
    YList.add(*(points.at(ip).second));
  }

  TString name = "PDF";
  if (appendName!="") name = Form("%s_%s", name.Data(), appendName.Data());
  TString title=name;
  PDF = new RooNCSplinePdf_fast(
    name.Data(),
    title.Data(),
    *splineVar,
    XList,
    YList
    );
}
RooNCSplinePdf_fast* NCSplinePdfFactory::getPDF(){ return PDF; }
