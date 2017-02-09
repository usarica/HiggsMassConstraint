#include "SlicePdfFactory.h"

using namespace std;


SlicePdfFactory::SlicePdfFactory(RooAbsReal* sliceVar_, const RooArgList& depList_, const RooArgList& funcList_, const RooArgList& coefList_, TString appendName_) :
appendName(appendName_),
sliceVar(sliceVar_),
depList(depList_),
funcList(funcList_),
coefList(coefList_),
PDF(0)
{}
SlicePdfFactory::~SlicePdfFactory(){
  destroyPDF();
  destroyX();
}
void SlicePdfFactory::destroyX(){
  for (unsigned int ip=0; ip<points.size(); ip++) delete points.at(ip);
  points.clear();
}
void SlicePdfFactory::setXCoordinates(vector<Double_t> coordx){
  initX(coordx);
  initPDF();
}
void SlicePdfFactory::initX(vector<Double_t> coordx){
  // Flush the points array
  destroyX();

  for (unsigned int ip=0; ip<coordx.size(); ip++){
    TString xname = Form("point_%i_x", ip);
    if (appendName!="") xname = Form("%s_%s", xname.Data(), appendName.Data());
    RooConstVar* xvar = new RooConstVar(xname, xname, coordx[ip]);
    points.push_back(xvar);
  }
}

void SlicePdfFactory::destroyPDF(){ delete PDF; PDF=0; }
void SlicePdfFactory::initPDF(){
  destroyPDF();

  RooArgList XList;
  for (unsigned int ip=0; ip<points.size(); ip++) XList.add(*(points.at(ip)));

  TString name = "PDF";
  if (appendName!="") name = Form("%s_%s", name.Data(), appendName.Data());
  TString title=name;
  PDF = new RooSlicePdf(
    name.Data(),
    title.Data(),
    *sliceVar,
    depList,
    XList,
    funcList,
    coefList
    );
}
RooSlicePdf* SlicePdfFactory::getPDF(){ return PDF; }
