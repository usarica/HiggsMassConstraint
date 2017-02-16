#include <iostream>
#include <cmath>
#include <string>
#include <algorithm>
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussModel.h"
#include "RooRealIntegral.h"
#include "RooDecay.h"
#include "RooBMixDecay.h"
#include "RooCategory.h"
#include "RooBinning.h"
#include "RooPlot.h"
#include "TSystem.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"
#include "TGaxis.h"
#include "TString.h"
#include "TChain.h"
#include <ZZMatrixElement/MELA/interface/Mela.h>
#include <ZZMatrixElement/MELA/interface/ScalarPdfFactory_HVV.h>

using namespace RooFit;
using namespace std;


// Global variables
const double m1m2_BinWidth = 1.;
const double m1Range[2]={ 40., 120. };
const double m2Range[2]={ 12., 120. };
const Int_t intCodeStart = RooSpin::prime_h1*RooSpin::prime_h2*RooSpin::prime_Phi*RooSpin::prime_hs*RooSpin::prime_Phi1;
const bool recordPDFHisto=false;

void splitOption(const string rawoption, string& wish, string& value, char delimiter);
void splitOptionRecursive(const string rawoption, vector<string>& splitoptions, char delimiter);
Bool_t checkListVariable(const vector<string>& list, const string& var);
string set_gHVV(RooSpinZero::modelCouplings& couplings, string gSet);
void getIntegral_single_HVV(string gSet, string decaytype);

void getIntegral(){
  getIntegral_single_HVV("ghz1=1,0", "ZZ4l");
}


void getIntegral_single_HVV(string gSet, string decaytype){
  RooRealVar* m12 = new RooRealVar("m12", "M_{H} (GeV)", 70., 20000.);
  RooRealVar* m1 = new RooRealVar("m1", "m_{1} (GeV)", m1Range[0], m1Range[0], m1Range[1]);
  RooRealVar* m2 = new RooRealVar("m2", "m_{2} (GeV)", m2Range[0], m2Range[0], m2Range[1]);
  m1->setBins((m1->getMax()-m1->getMin())/m1m2_BinWidth);
  m2->setBins((m2->getMax()-m2->getMin())/m1m2_BinWidth);
  RooRealVar* hs = new RooRealVar("hs", "cos#theta^{*}", -1, 1);
  RooRealVar* h1 = new RooRealVar("h1", "cos#theta_{1}", -1, 1);
  RooRealVar* h2 = new RooRealVar("h2", "cos#theta_{2}", -1, 1);
  RooRealVar* Phi = new RooRealVar("phi", "#Phi", -TMath::Pi(), TMath::Pi());
  RooRealVar* Phi1 = new RooRealVar("phi1", "#Phi_{1}", -TMath::Pi(), TMath::Pi());
  RooRealVar* Y = new RooRealVar("Y", "Y", 0);

  RooSpin::modelMeasurables measurables;
  measurables.m1 = m1;
  measurables.m2 = m2;
  measurables.m12 = m12;
  if (intCodeStart%RooSpin::prime_h1 != 0) measurables.h1 = h1;
  else measurables.h1 = 0;
  if (intCodeStart%RooSpin::prime_h2 != 0) measurables.h2 = h2;
  else measurables.h2 = 0;
  if (intCodeStart%RooSpin::prime_Phi != 0) measurables.Phi = Phi;
  else measurables.Phi = 0;
  if (intCodeStart%RooSpin::prime_hs != 0) measurables.hs = hs;
  else measurables.hs = 0;
  if (intCodeStart%RooSpin::prime_Phi1 != 0) measurables.Phi1 = Phi1;
  else measurables.Phi1 = 0;
  measurables.Y = Y;

  RooSpin::VdecayType Vdecay1=(decaytype.find("WW")!=string::npos ? RooSpin::kVdecayType_Wany : RooSpin::kVdecayType_Zll);
  RooSpin::VdecayType Vdecay2=(decaytype.find("WW")!=string::npos ? RooSpin::kVdecayType_Wany : RooSpin::kVdecayType_Zll);
  if (decaytype.find("ZG")!=string::npos){
    Vdecay2=RooSpin::kVdecayType_GammaOnshell;
    if (decaytype.find("2l")!=string::npos) Vdecay1=RooSpin::kVdecayType_Zll;
    else if (decaytype.find("2nu")!=string::npos) Vdecay1=RooSpin::kVdecayType_Znn;
    else if (decaytype.find("2q")!=string::npos) Vdecay1=RooSpin::kVdecayType_Zud;
    else{
      cerr << "Could not find the Z decays! Exiting." << endl;
      return;
    }
  }
  else if (decaytype.find("GG")!=string::npos && decaytype.find("GGto")==string::npos){ Vdecay1=RooSpin::kVdecayType_GammaOnshell; Vdecay2=RooSpin::kVdecayType_GammaOnshell; }
  else if (decaytype.find("ZZ")!=string::npos){
    if (decaytype.find("4l")!=string::npos){ Vdecay1=RooSpin::kVdecayType_Zll; Vdecay2=RooSpin::kVdecayType_Zll; }
    else if (decaytype.find("4nu")!=string::npos){ Vdecay1=RooSpin::kVdecayType_Znn; Vdecay2=RooSpin::kVdecayType_Znn; }
    else if (decaytype.find("4q")!=string::npos){ Vdecay1=RooSpin::kVdecayType_Zud; Vdecay2=RooSpin::kVdecayType_Zud; }
    else if (decaytype.find("2l2q")!=string::npos || decaytype.find("2q2l")!=string::npos){ Vdecay1=RooSpin::kVdecayType_Zll; Vdecay2=RooSpin::kVdecayType_Zud; }
    else if (decaytype.find("2l2nu")!=string::npos || decaytype.find("2nu2l")!=string::npos){ Vdecay1=RooSpin::kVdecayType_Zll; Vdecay2=RooSpin::kVdecayType_Znn; }
    else if (decaytype.find("2q2nu")!=string::npos || decaytype.find("2nu2q")!=string::npos){ Vdecay1=RooSpin::kVdecayType_Zll; Vdecay2=RooSpin::kVdecayType_Zud; }
    else{
      cerr << "Could not find the Z decays! Exiting." << endl;
      return;
    }
  }
  else if (decaytype.find("WW")==string::npos){
    cerr << "Could not find the V decays! Exiting." << endl;
    return;
  }
  if (Vdecay1==RooSpin::kVdecayType_GammaOnshell){ m1->setVal(0); m1->setConstant(true); }
  if (Vdecay2==RooSpin::kVdecayType_GammaOnshell){ m2->setVal(0); m2->setConstant(true); }

  cout << "Decay modes found are " << Vdecay1 << '\t' << Vdecay2 << endl;
  ScalarPdfFactory_HVV* someHiggs = new ScalarPdfFactory_HVV(measurables, false, Vdecay1, Vdecay2, true);

  // Set coouplings from string
  someHiggs->setZZ4fOrdering(false);
  someHiggs->makeParamsConst(false);
  string strcoupl = set_gHVV(someHiggs->couplings, gSet);
  someHiggs->makeParamsConst(true);

  RooSpinZero_7DComplex_withAccep_HVV* pdf = (RooSpinZero_7DComplex_withAccep_HVV*)someHiggs->getPDF();
  pdf->alwaysIntegrate(intCodeStart);
  RooRealIntegral pdf_int("pdf_int", "", *pdf, RooArgSet(*(measurables.m1), *(measurables.m2)));

  //RooAbsReal* pdf_int = pdf->createIntegral(RooArgSet(*(measurables.m1), *(measurables.m2)));
  //cout << "pdf integral: " << pdf_int->getVal() << endl;
  //pdf_int->Print("v");
  //delete pdf_int;

  TFile* foutput = TFile::Open(Form("H%sDecay_mHstarShape_NoInterf_%s%s", decaytype.c_str(), strcoupl.c_str(), ".root"), "recreate");
  cout << "Will compute mH=\n";
  vector<double> masses;
  double massmin=(int)(m1Range[0]+m2Range[0]+m1m2_BinWidth+0.5);
  double mass=massmin;
  while (mass<=15000){
    cout << mass << "\n";
    masses.push_back(mass);
    double massinc;
    if (mass<200.) massinc=1;
    else if (mass<600.) massinc=2.;
    else if (mass<1500.) massinc=5.;
    else if (mass<3000.) massinc=10.;
    else if (mass<10000.) massinc=50.;
    else massinc=100.;
    mass += massinc;
  }
  const unsigned int nbins=masses.size();
  cout << "Total number of mass points to evaluate: " << nbins << endl;
  cout << endl;
  double* xyarray[2]={ new double[nbins], new double[nbins] };
  for (unsigned int ibin=0; ibin<nbins; ibin++){
    double mPOLE = masses.at(ibin);
    m12->setConstant(false);
    m12->setVal(mPOLE);
    m12->setConstant(true);
    xyarray[0][ibin]=mPOLE;
    xyarray[1][ibin]=pdf_int.getVal();
    cout << "PDF[mH=" << xyarray[0][ibin] << "] = " << xyarray[1][ibin] << endl;
    if (recordPDFHisto){
      TH2F* hpdf = (TH2F*)pdf->createHistogram(Form("h_pdf_mH%.0f", mPOLE), *m1, YVar(*m2));
      hpdf->SetOption("colz");
      for (int ix=1; ix<=hpdf->GetNbinsX(); ix++){
        for (int iy=1; iy<=hpdf->GetNbinsY(); iy++){
          if ((hpdf->GetXaxis()->GetBinCenter(ix)+hpdf->GetYaxis()->GetBinCenter(iy))>xyarray[0][ibin]) hpdf->SetBinContent(ix, iy, 0.);
        }
      }
      hpdf->Scale(xyarray[1][ibin]/hpdf->Integral("width"));
      //cout << "Histogram integral: " << hpdf->Integral("width") << " or " << hpdf->Integral() << " (no width)" << endl;
      foutput->WriteTObject(hpdf);
      delete hpdf;
    }
  }
  TGraph* tgint = new TGraph(nbins, xyarray[0], xyarray[1]);
  tgint->SetName("tg_pdf_mH");
  tgint->SetTitle("Integral at different mH");
  tgint->GetXaxis()->SetTitle(m12->GetTitle());
  tgint->GetYaxis()->SetTitle(Form("H%s amplitude",decaytype.c_str()));
  foutput->WriteTObject(tgint);
  delete tgint;
  foutput->Close();

  for (unsigned int ixy=0; ixy<2; ixy++) delete[] xyarray[ixy];
  delete someHiggs;
  delete m12;
  delete m1;
  delete m2;
  delete hs;
  delete h1;
  delete h2;
  delete Phi;
  delete Phi1;
  delete Y;
}

string set_gHVV(RooSpinZero::modelCouplings& couplings, string gSet){ // gSet in the form ghz1=1,0;ghz2=0,1 etc.
  string strcoupl="";
  vector<string> vtmp;
  splitOptionRecursive(gSet, vtmp, ';');
  for (unsigned int iopt=0; iopt<vtmp.size(); iopt++){
    string opt=vtmp.at(iopt);
    string wish, strVal, strValRe, strValIm;
    // Use double precision for couplings
    Double_t valRe=0;
    Double_t valIm=0;
    splitOption(opt, wish, strVal, '=');
    // Lambda and cz/cw couplings have no imaginary components, so do not expect to parse them with ','.
    if (
      wish.find("Lambda")==string::npos
      &&
      wish.find("q1sq")==string::npos
      &&
      wish.find("q2sq")==string::npos
      &&
      wish.find("q12sq")==string::npos
      &&
      wish!="separateWWZZcouplings"
      ){
      splitOption(strVal, strValRe, strValIm, ',');
      valRe = atof(strValRe.c_str());
      valIm = atof(strValIm.c_str());
    }
    else valRe = atof(strVal.c_str());

    strcoupl += wish;
    cout << "Coupling " << wish << " is set to " << valRe;
    if (valIm!=0.) cout << "," << valIm;
    cout << endl;

    int ig=-1;
    int ilam=0;

    if (wish.find("prime2")!=string::npos) ilam=2;
    else if (wish.find("prime3")!=string::npos) ilam=3;
    else if (wish.find("prime4")!=string::npos) ilam=4;
    else if (wish.find("prime5")!=string::npos) ilam=5;
    else if (wish.find("prime6")!=string::npos) ilam=6;
    else if (wish.find("prime7")!=string::npos) ilam=7;
    else if (wish.find("prime")!=string::npos) ilam=1;

    if (wish.find("ghz1")!=string::npos) ig=0;
    else if (wish.find("ghz2")!=string::npos) ig=1;
    else if (wish.find("ghz3")!=string::npos) ig=2;
    else if (wish.find("ghz4")!=string::npos) ig=3;
    else if (wish.find("ghzgs1")!=string::npos){ ig=4; ilam=2; }
    else if (wish.find("ghzgs2")!=string::npos){ ig=5; ilam=0; }
    else if (wish.find("ghzgs3")!=string::npos){ ig=6; ilam=0; }
    else if (wish.find("ghzgs4")!=string::npos){ ig=7; ilam=0; }
    else if (wish.find("ghgsgs2")!=string::npos){ ig=8; ilam=0; }
    else if (wish.find("ghgsgs3")!=string::npos){ ig=9; ilam=0; }
    else if (wish.find("ghgsgs4")!=string::npos){ ig=10; ilam=0; }


    if (wish=="cz_q1sq"){ ((RooRealVar*)couplings.cLambda_qsq[cLambdaHIGGS_VV_QSQ1])->setVal((int)valRe); }
    else if (wish=="cz_q2sq"){ ((RooRealVar*)couplings.cLambda_qsq[cLambdaHIGGS_VV_QSQ2])->setVal((int)valRe); }
    else if (wish=="cz_q12sq"){ ((RooRealVar*)couplings.cLambda_qsq[cLambdaHIGGS_VV_QSQ12])->setVal((int)valRe); }

    else if (wish=="Lambda_z10"){ ((RooRealVar*)couplings.Lambda_z1qsq[cLambdaHIGGS_VV_QSQ12])->setVal((int)valRe); }
    else if (wish=="Lambda_z11"){ ((RooRealVar*)couplings.Lambda_z1qsq[cLambdaHIGGS_VV_QSQ1])->setVal((int)valRe); }
    else if (wish=="Lambda_z12"){ ((RooRealVar*)couplings.Lambda_z1qsq[cLambdaHIGGS_VV_QSQ2])->setVal((int)valRe); }
    else if (wish=="Lambda_z20"){ ((RooRealVar*)couplings.Lambda_z2qsq[cLambdaHIGGS_VV_QSQ12])->setVal((int)valRe); }
    else if (wish=="Lambda_z21"){ ((RooRealVar*)couplings.Lambda_z2qsq[cLambdaHIGGS_VV_QSQ1])->setVal((int)valRe); }
    else if (wish=="Lambda_z22"){ ((RooRealVar*)couplings.Lambda_z2qsq[cLambdaHIGGS_VV_QSQ2])->setVal((int)valRe); }
    else if (wish=="Lambda_z30"){ ((RooRealVar*)couplings.Lambda_z3qsq[cLambdaHIGGS_VV_QSQ12])->setVal((int)valRe); }
    else if (wish=="Lambda_z31"){ ((RooRealVar*)couplings.Lambda_z3qsq[cLambdaHIGGS_VV_QSQ1])->setVal((int)valRe); }
    else if (wish=="Lambda_z32"){ ((RooRealVar*)couplings.Lambda_z3qsq[cLambdaHIGGS_VV_QSQ2])->setVal((int)valRe); }
    else if (wish=="Lambda_z40"){ ((RooRealVar*)couplings.Lambda_z4qsq[cLambdaHIGGS_VV_QSQ12])->setVal((int)valRe); }
    else if (wish=="Lambda_z41"){ ((RooRealVar*)couplings.Lambda_z4qsq[cLambdaHIGGS_VV_QSQ1])->setVal((int)valRe); }
    else if (wish=="Lambda_z42"){ ((RooRealVar*)couplings.Lambda_z4qsq[cLambdaHIGGS_VV_QSQ2])->setVal((int)valRe); }

    else{
      if (ig==0){
        ((RooRealVar*)couplings.g1List[ilam][0])->setVal(valRe);
        ((RooRealVar*)couplings.g1List[ilam][1])->setVal(valIm);
      }
      else if (ig==1){
        ((RooRealVar*)couplings.g2List[ilam][0])->setVal(valRe);
        ((RooRealVar*)couplings.g2List[ilam][1])->setVal(valIm);
      }
      else if (ig==2){
        ((RooRealVar*)couplings.g3List[ilam][0])->setVal(valRe);
        ((RooRealVar*)couplings.g3List[ilam][1])->setVal(valIm);
      }
      else if (ig==3){
        ((RooRealVar*)couplings.g4List[ilam][0])->setVal(valRe);
        ((RooRealVar*)couplings.g4List[ilam][1])->setVal(valIm);
      }
      else if (ig==4){
        ((RooRealVar*)couplings.gzgs1List[ilam-2][0])->setVal(valRe);
        ((RooRealVar*)couplings.gzgs1List[ilam-2][1])->setVal(valIm);
      }
      else if (ig==5){
        ((RooRealVar*)couplings.gzgs2List[ilam][0])->setVal(valRe);
        ((RooRealVar*)couplings.gzgs2List[ilam][1])->setVal(valIm);
      }
      else if (ig==6){
        ((RooRealVar*)couplings.gzgs3List[ilam][0])->setVal(valRe);
        ((RooRealVar*)couplings.gzgs3List[ilam][1])->setVal(valIm);
      }
      else if (ig==7){
        ((RooRealVar*)couplings.gzgs4List[ilam][0])->setVal(valRe);
        ((RooRealVar*)couplings.gzgs4List[ilam][1])->setVal(valIm);
      }
      else if (ig==8){
        ((RooRealVar*)couplings.ggsgs2List[ilam][0])->setVal(valRe);
        ((RooRealVar*)couplings.ggsgs2List[ilam][1])->setVal(valIm);
      }
      else if (ig==9){
        ((RooRealVar*)couplings.ggsgs3List[ilam][0])->setVal(valRe);
        ((RooRealVar*)couplings.ggsgs3List[ilam][1])->setVal(valIm);
      }
      else if (ig==10){
        ((RooRealVar*)couplings.ggsgs4List[ilam][0])->setVal(valRe);
        ((RooRealVar*)couplings.ggsgs4List[ilam][1])->setVal(valIm);
      }
    }

  }
  return strcoupl;
}


void splitOption(const string rawoption, string& wish, string& value, char delimiter){
  size_t posEq = rawoption.find(delimiter);
  if (posEq!=string::npos){
    wish=rawoption;
    value=rawoption.substr(posEq+1);
    wish.erase(wish.begin()+posEq, wish.end());
  }
  else{
    wish="";
    value=rawoption;
  }
}
void splitOptionRecursive(const string rawoption, vector<string>& splitoptions, char delimiter){
  string suboption=rawoption, result=rawoption;
  string remnant;
  while (result!=""){
    splitOption(suboption, result, remnant, delimiter);
    if (result!="" && !checkListVariable(splitoptions, result)) splitoptions.push_back(result);
    suboption = remnant;
  }
  if (remnant!="") splitoptions.push_back(remnant);
}
Bool_t checkListVariable(const vector<string>& list, const string& var){
  for (unsigned int v=0; v<list.size(); v++){
    if (list.at(v)==var) return true; // Look for exact match
  }
  return false;
}
