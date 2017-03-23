#include "NCSplinePdfFactory_1D.h"
#include <cassert>

using namespace std;


NCSplinePdfFactory_1D::NCSplinePdfFactory_1D(RooAbsReal& splineVar_, TString appendName_) :
appendName(appendName_),
splineVar(&splineVar_),
PDF(0)
{}
NCSplinePdfFactory_1D::~NCSplinePdfFactory_1D(){
  destroyPDF();
}
void NCSplinePdfFactory_1D::setPoints(TGraph* tg){
  vector<Double_t> XList, FcnList;
  double* xx = tg->GetX();
  double* yy = tg->GetY();
  int n = tg->GetN();
  for (int ip=0; ip<n; ip++){ XList.push_back(xx[ip]); FcnList.push_back(yy[ip]); }
  setPoints(XList, FcnList);
}
const std::vector<std::pair<RooNCSplinePdfCore::T, RooNCSplinePdfCore::T>> NCSplinePdfFactory_1D::getPoints(const std::vector<RooNCSplinePdfCore::T>& XList, const std::vector<RooNCSplinePdfCore::T>& FcnList){
  unsigned int nX = XList.size();
  unsigned int n = FcnList.size();
  if (nX!=n){
    cerr << "NCSplinePdfFactory_1D::getPoints: nX=" << nX << " != nFcn=" << n << endl;
    assert(0);
  }
  std::vector<std::pair<RooNCSplinePdfCore::T, RooNCSplinePdfCore::T>> pList; pList.reserve(n);
  for (unsigned int ip=0; ip<n; ip++) pList.push_back(pair<RooNCSplinePdfCore::T, RooNCSplinePdfCore::T>(XList.at(ip), FcnList.at(ip)));
  return pList;
}

void NCSplinePdfFactory_1D::destroyPDF(){ delete PDF; PDF=0; }
void NCSplinePdfFactory_1D::initPDF(const std::vector<std::pair<RooNCSplinePdfCore::T, RooNCSplinePdfCore::T>>& pList){
  destroyPDF();

  std::vector<RooNCSplinePdfCore::T> XList;
  std::vector<RooNCSplinePdfCore::T> FcnList;
  for (unsigned int ip=0; ip<pList.size(); ip++){
    XList.push_back(pList.at(ip).first);
    FcnList.push_back(pList.at(ip).second);
  }

  TString name = "PDF";
  if (appendName!="") name = Form("%s_%s", name.Data(), appendName.Data());
  TString title=name;
  PDF = new RooNCSplinePdf_1D_fast(
    name.Data(),
    title.Data(),
    *splineVar,
    XList, FcnList
    );
}
RooNCSplinePdf_1D_fast* NCSplinePdfFactory_1D::getPDF(){ return PDF; }
