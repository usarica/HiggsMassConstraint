#include "NCSplinePdfFactory_3D.h"

using namespace std;


NCSplinePdfFactory_3D::NCSplinePdfFactory_3D(RooAbsReal& XVar_, RooAbsReal& YVar_, RooAbsReal& ZVar_, TString appendName_) :
appendName(appendName_),
XVar(&XVar_), YVar(&YVar_), ZVar(&ZVar_),
PDF(0)
{}
NCSplinePdfFactory_3D::~NCSplinePdfFactory_3D(){
  destroyPDF();
}

void NCSplinePdfFactory_3D::addUnique(std::vector<RooNCSplinePdfCore::T>& list, RooNCSplinePdfCore::T val){
  for (unsigned int ip=0; ip<list.size(); ip++){ if (list.at(ip)==val) return; }
  list.push_back(val);
}
const std::vector<doubleQuadruplet_t> NCSplinePdfFactory_3D::getPoints(
  const std::vector<RooNCSplinePdfCore::T>& XList,
  const std::vector<RooNCSplinePdfCore::T>& YList,
  const std::vector<RooNCSplinePdfCore::T>& ZList,
  const std::vector<RooNCSplinePdfCore::T>& FcnList
  ){
  unsigned int nX = XList.size();
  unsigned int nY = YList.size();
  unsigned int nZ = ZList.size();
  unsigned int n = FcnList.size();
  if (nX*nY*nZ!=n){
    cerr << "NCSplinePdfFactory_3D::getPoints: nX=" << nX << " x nY=" << nY << " x nZ=" << nZ << " != nFcn=" << n << endl;
    assert(0);
  }

  std::vector<doubleQuadruplet_t> pList; pList.reserve(n);
  for (unsigned int ix=0; ix<nX; ix++){
    RooNCSplinePdfCore::T xval = XList.at(ix);
    for (unsigned int iy=0; iy<nY; iy++){
      RooNCSplinePdfCore::T yval = YList.at(iy);
      for (unsigned int iz=0; iz<nZ; iz++){
        RooNCSplinePdfCore::T zval = ZList.at(iz);
        unsigned int ip = nZ*(nY*ix + iy) + iz;
        pList.push_back(doubleQuadruplet_t(xval, yval, zval, FcnList.at(ip)));
      }
    }
  }
  return pList;
}

void NCSplinePdfFactory_3D::destroyPDF(){ delete PDF; PDF=0; }
void NCSplinePdfFactory_3D::initPDF(const std::vector<doubleQuadruplet_t>& pList){
  destroyPDF();

  unsigned int n = pList.size();
  vector<RooNCSplinePdfCore::T> XList;
  vector<RooNCSplinePdfCore::T> YList;
  vector<RooNCSplinePdfCore::T> ZList;
  vector<vector<vector<RooNCSplinePdfCore::T>>> FcnList;
  for (unsigned int ip=0; ip<n; ip++){
    addUnique(XList, (pList.at(ip))[0]);
    addUnique(YList, (pList.at(ip))[1]);
    addUnique(ZList, (pList.at(ip))[2]);
  }
  FcnList.reserve(ZList.size());
  for (unsigned int iz=0; iz<ZList.size(); iz++){
    vector<vector<RooNCSplinePdfCore::T>> dumz;
    dumz.reserve(YList.size());
    for (unsigned int iy=0; iy<YList.size(); iy++){
      vector<RooNCSplinePdfCore::T> dumy;
      dumy.reserve(XList.size());
      for (unsigned int ix=0; ix<XList.size(); ix++){
        unsigned int ip = ZList.size()*(YList.size()*ix + iy) + iz;
        dumy.push_back((pList.at(ip))[3]); // Do not use unique here
      }
      dumz.push_back(dumy);
    }
    FcnList.push_back(dumz);
  }

  TString name = "PDF";
  if (appendName!="") name = Form("%s_%s", name.Data(), appendName.Data());
  TString title=name;
  PDF = new RooNCSplinePdf_3D_fast(
    name.Data(),
    title.Data(),
    *XVar, *YVar, *ZVar,
    XList, YList, ZList, FcnList
    );
}
RooNCSplinePdf_3D_fast* NCSplinePdfFactory_3D::getPDF(){ return PDF; }
