#include "NCSplinePdfFactory_2D.h"

using namespace std;


NCSplinePdfFactory_2D::NCSplinePdfFactory_2D(RooAbsReal* XVar_, RooAbsReal* YVar_, TString appendName_) :
appendName(appendName_),
XVar(XVar_), YVar(YVar_),
PDF(0)
{}
NCSplinePdfFactory_2D::~NCSplinePdfFactory_2D(){
  destroyPDF();
  destroyPoints();
}
void NCSplinePdfFactory_2D::destroyPoints(){
  for (unsigned int ip=0; ip<Xcoord.size(); ip++) delete Xcoord.at(ip);
  for (unsigned int ip=0; ip<Ycoord.size(); ip++) delete Ycoord.at(ip);
  for (unsigned int ip=0; ip<FcnVal.size(); ip++) delete FcnVal.at(ip);
  Xcoord.clear(); Ycoord.clear(); FcnVal.clear();
}
void NCSplinePdfFactory_2D::setPoints(const std::vector<doubleTriplet_t>& pList){
  initPoints(pList);
  initPDF();
}
void NCSplinePdfFactory_2D::addUnique(std::vector<Double_t>& list, Double_t val){
  for (unsigned int ip=0; ip<list.size(); ip++){ if (list.at(ip)==val) return; }
  list.push_back(val);
}
void NCSplinePdfFactory_2D::initPoints(const std::vector<doubleTriplet_t>& pList){
  // Flush the points array
  destroyPoints();

  unsigned int n = pList.size();
  vector<Double_t> Xval;
  vector<Double_t> Yval;
  vector<Double_t> Fcnval;
  for (unsigned int ip=0; ip<n; ip++){
    addUnique(Xval, (pList.at(ip))[0]);
    addUnique(Yval, (pList.at(ip))[1]);
    Fcnval.push_back((pList.at(ip))[2]); // Do not use unique here
  }
  for (unsigned int ip=0; ip<Xval.size(); ip++){
    TString name = Form("point_x_%i", ip);
    if (appendName!="") name = Form("%s_%s", name.Data(), appendName.Data());
    RooConstVar* var = new RooConstVar(name, name, Xval.at(ip));
    Xcoord.push_back(var);
  }
  for (unsigned int ip=0; ip<Yval.size(); ip++){
    TString name = Form("point_y_%i", ip);
    if (appendName!="") name = Form("%s_%s", name.Data(), appendName.Data());
    RooConstVar* var = new RooConstVar(name, name, Yval.at(ip));
    Ycoord.push_back(var);
  }
  for (unsigned int ip=0; ip<Fcnval.size(); ip++){
    TString name = Form("point_fcn_%i", ip);
    if (appendName!="") name = Form("%s_%s", name.Data(), appendName.Data());
    RooConstVar* var = new RooConstVar(name, name, Fcnval.at(ip));
    FcnVal.push_back(var);
  }
}

void NCSplinePdfFactory_2D::destroyPDF(){ delete PDF; PDF=0; }
void NCSplinePdfFactory_2D::initPDF(){
  destroyPDF();

  RooArgList XList;
  RooArgList YList;
  RooArgList FcnList;
  for (unsigned int ip=0; ip<Xcoord.size(); ip++) XList.add(*(Xcoord.at(ip)));
  for (unsigned int ip=0; ip<Ycoord.size(); ip++) YList.add(*(Ycoord.at(ip)));
  for (unsigned int ip=0; ip<FcnVal.size(); ip++) FcnList.add(*(FcnVal.at(ip)));

  TString name = "PDF";
  if (appendName!="") name = Form("%s_%s", name.Data(), appendName.Data());
  TString title=name;
  PDF = new RooNCSplinePdf_2D_fast(
    name.Data(),
    title.Data(),
    *XVar, *YVar,
    XList,
    YList,
    FcnList
    );
}
RooNCSplinePdf_2D_fast* NCSplinePdfFactory_2D::getPDF(){ return PDF; }
