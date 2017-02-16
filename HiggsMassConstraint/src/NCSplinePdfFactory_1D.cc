#include "NCSplinePdfFactory_1D.h"

using namespace std;


NCSplinePdfFactory_1D::NCSplinePdfFactory_1D(RooAbsReal* splineVar_, TString appendName_) :
appendName(appendName_),
splineVar(splineVar_),
PDF(0)
{}
NCSplinePdfFactory_1D::~NCSplinePdfFactory_1D(){
  destroyPDF();
  destroyPoints();
}
void NCSplinePdfFactory_1D::destroyPoints(){
  for (unsigned int ip=0; ip<points.size(); ip++){
    delete points.at(ip).first; delete points.at(ip).second;
  }
  points.clear();
}
void NCSplinePdfFactory_1D::setPoints(TGraph* tg){
  vector<Double_t> XList, FcnList;
  double* xx = tg->GetX();
  double* yy = tg->GetY();
  int n = tg->GetN();
  for (int ip=0; ip<n; ip++){ XList.push_back(xx[ip]); FcnList.push_back(yy[ip]); }
  setPoints(XList, FcnList);
}
void NCSplinePdfFactory_1D::setPoints(const std::vector<Double_t>& XList, const std::vector<Double_t>& FcnList){
  initPoints(XList, FcnList);
  initPDF();
}
void NCSplinePdfFactory_1D::initPoints(const std::vector<Double_t>& XList, const std::vector<Double_t>& FcnList){
  // Flush the points array
  destroyPoints();

  unsigned int n = XList.size();
  for (unsigned int ip=0; ip<n; ip++){
    TString xname = Form("point_%i_x", ip);
    TString yname = Form("point_%i_y", ip);
    if (appendName!=""){
      xname = Form("%s_%s", xname.Data(), appendName.Data());
      yname = Form("%s_%s", yname.Data(), appendName.Data());
    }
    RooConstVar* xvar = new RooConstVar(xname, xname, XList.at(ip));
    RooConstVar* yvar = new RooConstVar(yname, yname, FcnList.at(ip));
    points.push_back(pair<RooConstVar*, RooConstVar*>(xvar, yvar));
  }
}

void NCSplinePdfFactory_1D::destroyPDF(){ delete PDF; PDF=0; }
void NCSplinePdfFactory_1D::initPDF(){
  destroyPDF();

  RooArgList XList, YList;
  for (unsigned int ip=0; ip<points.size(); ip++){
    XList.add(*(points.at(ip).first));
    YList.add(*(points.at(ip).second));
  }

  TString name = "PDF";
  if (appendName!="") name = Form("%s_%s", name.Data(), appendName.Data());
  TString title=name;
  PDF = new RooNCSplinePdf_1D_fast(
    name.Data(),
    title.Data(),
    *splineVar,
    XList,
    YList
    );
}
RooNCSplinePdf_1D_fast* NCSplinePdfFactory_1D::getPDF(){ return PDF; }
