#include "NCSplinePdfFactory_3D.h"

using namespace std;


NCSplinePdfFactory_3D::NCSplinePdfFactory_3D(RooAbsReal* XVar_, RooAbsReal* YVar_, RooAbsReal* ZVar_, TString appendName_) :
appendName(appendName_),
XVar(XVar_), YVar(YVar_), ZVar(ZVar_),
PDF(0)
{}
NCSplinePdfFactory_3D::~NCSplinePdfFactory_3D(){
  destroyPDF();
  destroyPoints();
}
void NCSplinePdfFactory_3D::destroyPoints(){
  for (unsigned int ip=0; ip<Xcoord.size(); ip++) delete Xcoord.at(ip);
  for (unsigned int ip=0; ip<Ycoord.size(); ip++) delete Ycoord.at(ip);
  for (unsigned int ip=0; ip<Zcoord.size(); ip++) delete Zcoord.at(ip);
  for (unsigned int iz=0; iz<FcnVal.size(); iz++){
    for (unsigned int iy=0; iy<FcnVal.at(iz).size(); iy++){
      for (unsigned int ix=0; ix<FcnVal.at(iz).at(iy).size(); ix++) delete FcnVal.at(iz).at(iy).at(ix);
      FcnVal.at(iz).at(iy).clear();
    }
    FcnVal.at(iz).clear();
  }
  Xcoord.clear(); Ycoord.clear(); Zcoord.clear(); FcnVal.clear();
}
void NCSplinePdfFactory_3D::setPoints(const std::vector<doubleQuadruplet_t>& pList){
  initPoints(pList);
  initPDF();
}
void NCSplinePdfFactory_3D::addUnique(std::vector<Double_t>& list, Double_t val){
  for (unsigned int ip=0; ip<list.size(); ip++){ if (list.at(ip)==val) return; }
  list.push_back(val);
}
void NCSplinePdfFactory_3D::initPoints(const std::vector<doubleQuadruplet_t>& pList){
  // Flush the points array
  destroyPoints();

  unsigned int n = pList.size();

  //cout << "NCSplinePdfFactory_3D::initPoints: Received " << n << " points" << endl;

  vector<Double_t> Xval;
  vector<Double_t> Yval;
  vector<Double_t> Zval;
  vector<vector<vector<Double_t>>> Fcnval;
  for (unsigned int ip=0; ip<n; ip++){
    addUnique(Xval, (pList.at(ip))[0]);
    addUnique(Yval, (pList.at(ip))[1]);
    addUnique(Zval, (pList.at(ip))[2]);
  }
  for (unsigned int iz=0; iz<Zval.size(); iz++){
    vector<vector<Double_t>> dumz;
    Fcnval.push_back(dumz);
    for (unsigned int iy=0; iy<Yval.size(); iy++){
      vector<Double_t> dumy;
      Fcnval.at(iz).push_back(dumy);
      for (unsigned int ix=0; ix<Xval.size(); ix++){
        unsigned int ip = Zval.size()*(Yval.size()*ix + iy) + iz;
        Fcnval.at(iz).at(iy).push_back((pList.at(ip))[3]);
      }
    }
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
  for (unsigned int ip=0; ip<Zval.size(); ip++){
    TString name = Form("point_z_%i", ip);
    if (appendName!="") name = Form("%s_%s", name.Data(), appendName.Data());
    RooConstVar* var = new RooConstVar(name, name, Zval.at(ip));
    Zcoord.push_back(var);
  }

  for (unsigned int iz=0; iz<Fcnval.size(); iz++){
    vector<vector<RooConstVar*>> dumz;
    FcnVal.push_back(dumz);
    for (unsigned int iy=0; iy<Fcnval.at(iz).size(); iy++){
      vector<RooConstVar*> dumy;
      FcnVal.at(iz).push_back(dumy);
      for (unsigned int ix=0; ix<Fcnval.at(iz).at(iy).size(); ix++){
        TString name = Form("point_fcn_Z%iY%iX%i", iz, iy, ix);
        if (appendName!="") name = Form("%s_%s", name.Data(), appendName.Data());
        RooConstVar* var = new RooConstVar(name, name, Fcnval.at(iz).at(iy).at(ix));
        FcnVal.at(iz).at(iy).push_back(var);
      }
    }
  }
}

void NCSplinePdfFactory_3D::destroyPDF(){ delete PDF; PDF=0; }
void NCSplinePdfFactory_3D::initPDF(){
  destroyPDF();

  RooArgList XList;
  RooArgList YList;
  RooArgList ZList;
  vector<vector<const RooArgList*>> FcnList;
  for (unsigned int ip=0; ip<Xcoord.size(); ip++) XList.add(*(Xcoord.at(ip)));
  for (unsigned int ip=0; ip<Ycoord.size(); ip++) YList.add(*(Ycoord.at(ip)));
  for (unsigned int ip=0; ip<Zcoord.size(); ip++) ZList.add(*(Zcoord.at(ip)));
  for (unsigned int iz=0; iz<FcnVal.size(); iz++){
    vector<const RooArgList*> dumz;
    FcnList.push_back(dumz);
    for (unsigned int iy=0; iy<FcnVal.at(iz).size(); iy++){
      RooArgList* list = new RooArgList();
      for (unsigned int ix=0; ix<FcnVal.at(iz).at(iy).size(); ix++) list->add(*(FcnVal.at(iz).at(iy).at(ix)));
      FcnList.at(iz).push_back(list);
    }
  }

  //cout << "NCSplinePdfFactory_3D::initPDF: Initializing an " << Xcoord.size() << " x " << Ycoord.size() << " x " << Zcoord.size() << " pdf" << endl;

  TString name = "PDF";
  if (appendName!="") name = Form("%s_%s", name.Data(), appendName.Data());
  TString title=name;
  PDF = new RooNCSplinePdf_3D_fast(
    name.Data(),
    title.Data(),
    XVar, YVar, ZVar,
    &XList, &YList, &ZList,
    FcnList
    );

  for (unsigned int iz=0; iz<FcnList.size(); iz++){
    for (unsigned int iy=0; iy<FcnList.at(iz).size(); iy++) delete FcnList.at(iz).at(iy);
    FcnList.at(iz).clear();
  }
}
RooNCSplinePdf_3D_fast* NCSplinePdfFactory_3D::getPDF(){ return PDF; }
