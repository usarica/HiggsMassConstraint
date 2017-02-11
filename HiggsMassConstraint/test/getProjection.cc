#include <iostream>
#include <cmath>
#include <string>
#include <algorithm>
#include "TSystem.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"
#include "TGaxis.h"
#include "TString.h"
#include "TChain.h"
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
#include "RooNumIntConfig.h"
#include "RooWorkspace.h"
#include <ZZMatrixElement/MELA/interface/Mela.h>
#include <ZZMatrixElement/MELA/interface/ScalarPdfFactory_HVV.h>
#include "NCSplinePdfFactory.h"
#include "RooSlicePdf.h"

using namespace RooFit;
using namespace std;

#ifndef varm12range
#define varm12range 1
#endif
#ifndef varm1m2range
#define varm1m2range 1
#endif

// Global variables
const double m1m2_BinWidth = 1.;
const double m1Range[2]={ 40., 120. };
const double m2Range[2]={ 12., 120. };
const Int_t intCodeStart = RooSpin::prime_h1*RooSpin::prime_h2*RooSpin::prime_Phi*RooSpin::prime_hs*RooSpin::prime_Phi1;
const bool recordPDFHisto=false;

TString savePlot(TString name, TCanvas* c);
void splitOption(const string rawoption, string& wish, string& value, char delimiter);
void splitOptionRecursive(const string rawoption, vector<string>& splitoptions, char delimiter);
Bool_t checkListVariable(const vector<string>& list, const string& var);
void addByLowest(double val, std::vector<double>& valArray);
string set_gHVV(RooSpinZero::modelCouplings& couplings, string gSet);
void getProjection_single_HVV(string gSet, string decaytype, string strprojvar);

void getProjection(){
  gSystem->Exec("mkdir -p ./plots");
  getProjection_single_HVV("ghz1=1,0", "ZZ4l", "m1");
  getProjection_single_HVV("ghz1=1,0", "ZZ4l", "m2");
}

void getProjection_single_HVV(string gSet, string decaytype, string strprojvar){
  RooRealVar* projvar=0;
  RooRealVar* m12 = new RooRealVar("m12", "m_{H} (GeV)", m1Range[0]+m2Range[0], 20000.); if (strprojvar==m12->GetName()) projvar=m12;
  RooRealVar* m1 = new RooRealVar("m1", "m_{1} (GeV)", m1Range[0], m1Range[0], m1Range[1]); if (strprojvar==m1->GetName()) projvar=m1;
  RooRealVar* m2 = new RooRealVar("m2", "m_{2} (GeV)", m2Range[0], m2Range[0], m2Range[1]); if (strprojvar==m2->GetName()) projvar=m2;
  m1->setBins((m1->getMax()-m1->getMin())/m1m2_BinWidth);
  m2->setBins((m2->getMax()-m2->getMin())/m1m2_BinWidth);
  RooRealVar* hs = new RooRealVar("hs", "cos#theta^{*}", -1, 1); if (strprojvar==hs->GetName()) projvar=hs;
  RooRealVar* h1 = new RooRealVar("h1", "cos#theta_{1}", -1, 1); if (strprojvar==h1->GetName()) projvar=h1;
  RooRealVar* h2 = new RooRealVar("h2", "cos#theta_{2}", -1, 1); if (strprojvar==h2->GetName()) projvar=h2;
  RooRealVar* Phi = new RooRealVar("phi", "#Phi", -TMath::Pi(), TMath::Pi()); if (strprojvar==Phi->GetName()) projvar=Phi;
  RooRealVar* Phi1 = new RooRealVar("phi1", "#Phi_{1}", -TMath::Pi(), TMath::Pi()); if (strprojvar==Phi1->GetName()) projvar=Phi1;
  RooRealVar* Y = new RooRealVar("Y", "Y", 0); if (strprojvar==Y->GetName()) projvar=Y;
  if (projvar!=0){
    gSystem->Exec(Form("mkdir -p ./plots/%s", projvar->GetName()));

#if varm1m2range==1
    const unsigned int nevals=101;
#else
    unsigned int nevals=101;
    if (projvar==m1 || projvar==m2) nevals=201;
#endif
    const double projmin=projvar->getMin();
    const double projmax=projvar->getMax();
    const double projwidth=(projmax-projmin)/((const double)(nevals-1));

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
    measurables.Y = 0;

    RooArgSet intSet;
    if (projvar!=m1) intSet.add(*(measurables.m1));
    if (projvar!=m2) intSet.add(*(measurables.m2));
    if (projvar!=h1 && measurables.h1!=0) intSet.add(*(measurables.h1));
    if (projvar!=h2 && measurables.h2!=0) intSet.add(*(measurables.h2));
    if (projvar!=hs && measurables.hs!=0) intSet.add(*(measurables.hs));
    if (projvar!=Phi && measurables.Phi!=0) intSet.add(*(measurables.Phi));
    if (projvar!=Phi1 && measurables.Phi1!=0) intSet.add(*(measurables.Phi1));
    if (projvar!=Y && measurables.Y!=0) intSet.add(*(measurables.Y));

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
    //someHiggs->setZZ4fOrdering(false);
    someHiggs->makeParamsConst(false);
    string strcoupl = set_gHVV(someHiggs->couplings, gSet);
    someHiggs->makeParamsConst(true);

    RooSpinZero_7DComplex_withAccep_HVV* pdf = (RooSpinZero_7DComplex_withAccep_HVV*)someHiggs->getPDF();
    pdf->alwaysIntegrate(intCodeStart);
    pdf->defaultIntegratorConfig()->method1D().setLabel("RooAdaptiveGaussKronrodIntegrator1D");
    pdf->defaultIntegratorConfig()->getConfigSection("RooAdaptiveGaussKronrodIntegrator1D").setRealValue("maxSeg", 100);;
    pdf->defaultIntegratorConfig()->setEpsAbs(1e-5);
    pdf->defaultIntegratorConfig()->setEpsRel(1e-5);
    RooRealIntegral* pdf_int=0;
    if (intSet.getSize()>0) pdf_int = new RooRealIntegral(Form("pdf_proj_%s",projvar->GetName()), "", *pdf, intSet);
    if (pdf_int!=0){
      delete pdf_int; // Delete dummy integration. Need to create at each mPOLE to optimize integration.

      TFile* foutput = TFile::Open(Form("H%sDecay_%sProjection_NoInterf_%s%s", decaytype.c_str(), strprojvar.c_str(), strcoupl.c_str(), ".root"), "recreate");
      cout << "Will compute mH=\n";
      vector<double> masses;
#if varm12range==1
      double massmin=(int)(m1Range[0]+m2Range[0]+m1m2_BinWidth+0.5);
#else
      double massmin=65;
#endif
      double mass=massmin;
      TTree* masses_tree = new TTree("masses", "");
      masses_tree->Branch("mH", &mass);
#if varm12range==1
      while (mass<=15000){
        cout << mass << "\n";
        masses.push_back(mass); masses_tree->Fill();
        double massinc;
        if (mass<200.) massinc=1;
        else if (mass<600.) massinc=10.;
        else if (mass<1500.) massinc=50.;
        else if (mass<3000.) massinc=100.;
        else if (mass<10000.) massinc=500.;
        else massinc=1000.;
        mass += massinc;
      }
#else
      while (mass<=500){
        cout << mass << "\n";
        masses.push_back(mass); masses_tree->Fill();
        double massinc=2;
        mass += massinc;
      }
#endif
      const unsigned int nbins=masses.size();
      cout << "Total number of mass points to evaluate: " << nbins << "==" << masses_tree->GetEntries() << endl;
      cout << endl;
      foutput->WriteTObject(masses_tree);
      delete masses_tree;

      vector<double> projarray;
      for (unsigned int ip=0; ip<nevals; ip++){
        double projval=projmin+projwidth*ip;
        projarray.push_back(projval);
      }
#if varm1m2range==1
      if (projvar==m1 || projvar==m2){
        double mv, gamv;
        pdf->getMVGamV(&mv, &gamv);
        cout << "Found mV, gamV = " << mv << " " << gamv << endl;
        double projnewmin=mv-2.*gamv;
        double projnewmax=mv+2.*gamv;
        unsigned int nextrabins=41;
        if (projvar==m2){
          projnewmin=mv-1.*gamv;
          projnewmax=mv+1.*gamv;
          nextrabins=81;
        }
        double projnewwidth=(projnewmax-projnewmin)/((const double)(nextrabins-1));
        cout << "Adding extra " << nextrabins << " points from " << projnewmin << " to " << projnewmax  << " in steps of " << projnewwidth << endl;
        for (unsigned int ip=0; ip<nextrabins; ip++){
          double projnewval=projnewmin+projnewwidth*ip;
          addByLowest(projnewval, projarray);
        }
      }
#endif
      unsigned int nprojbins = projarray.size();

      for (unsigned int ip=0; ip<nprojbins; ip++) cout << projarray.at(ip) << '\n';
      cout << "Total number of projection points to evaluate: " << nprojbins << endl;
      cout << endl;

      double* xyarray[2]={ new double[nprojbins], new double[nprojbins] };

      TString gifCmd;
      for (unsigned int ibin=0; ibin<nbins; ibin++){
        double mPOLE = masses.at(ibin);
        cout << "Evaluating mH=" << mPOLE << endl;
        m12->setConstant(false);
        m12->setVal(mPOLE);
        m12->setConstant(true);
        pdf_int = new RooRealIntegral(Form("pdf_proj_%s", projvar->GetName()), "", *pdf, intSet);
        for (unsigned int ip=0; ip<nprojbins; ip++){
          xyarray[0][ip]=projarray.at(ip);
          projvar->setConstant(false);
          projvar->setVal(xyarray[0][ip]);
          projvar->setConstant(true);
          xyarray[1][ip]=pdf_int->getVal();
        }
        delete pdf_int;
        TGraph* tgint = new TGraph(nprojbins, xyarray[0], xyarray[1]);
        tgint->SetName(Form("tg_pdf_mH%.0f", mPOLE));
        tgint->SetTitle(Form("Projection at mH=%.0f GeV", mPOLE));
        tgint->GetXaxis()->SetTitle(projvar->GetTitle());
        tgint->GetYaxis()->SetTitle(Form("H%s amplitude", decaytype.c_str()));
        tgint->SetMarkerStyle(20);
        tgint->SetMarkerSize(0.8);

        if (projvar==m1 || projvar==m2){
          NCSplinePdfFactory* spFactory = new NCSplinePdfFactory(projvar, Form("mH%.0f", mPOLE));
          spFactory->setGraph(tgint);
          RooNCSplinePdf_1D_fast* spPDF = spFactory->getPDF();

          RooRealIntegral spPDFint("spPDFint", "", *spPDF, RooArgSet(*projvar));
          double spint = spPDFint.getVal();
          double* yy = tgint->GetY();
          for (int ip=0; ip<tgint->GetN(); ip++) yy[ip]/=spint;
          RooPlot* plot = projvar->frame((projvar->getMax()-projvar->getMin())/m1m2_BinWidth);
          plot->GetXaxis()->CenterTitle();
          plot->GetYaxis()->SetTitleOffset(1.2);
          plot->GetYaxis()->CenterTitle();
          plot->GetYaxis()->SetTitle(tgint->GetYaxis()->GetTitle());
          plot->GetXaxis()->SetTitle(tgint->GetXaxis()->GetTitle());
          plot->GetXaxis()->SetNdivisions(-505);
          plot->SetTitle(tgint->GetTitle());
          spPDF->plotOn(plot, LineColor(kRed), LineWidth(2));
          TCanvas* canvas = new TCanvas(Form("c_%s", tgint->GetName()), "", 600, 600);
          plot->Draw();
          tgint->Draw("psame");
          canvas->Modified();
          canvas->Update();
          foutput->WriteTObject(canvas);
          gifCmd += savePlot(projvar->GetName(), canvas);
          gifCmd += " ";
          canvas->Close();
        }

        foutput->WriteTObject(tgint);
        delete tgint;
      }
      foutput->Close();

      for (unsigned int ixy=0; ixy<2; ixy++) delete[] xyarray[ixy];

      if (gifCmd!=""){
        gifCmd.Prepend("convert ");
        gifCmd.Append(Form("./plots/%s/all.gif", projvar->GetName()));
        gSystem->Exec(gifCmd);
      }
    }

    delete someHiggs;
  }
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
void addByLowest(double val, std::vector<double>& valArray){
  bool inserted = false;
  for (std::vector<double>::iterator it = valArray.begin(); it<valArray.end(); it++){
    if (*it>val){
      inserted=true;
      valArray.insert(it, val);
      break;
    }
  }
  if (!inserted) valArray.push_back(val);
}

TString savePlot(TString name, TCanvas* c){
  TString plotname = Form("./plots/%s/%s.png", name.Data(), c->GetName());
  c->SaveAs(plotname);
  return plotname;
}