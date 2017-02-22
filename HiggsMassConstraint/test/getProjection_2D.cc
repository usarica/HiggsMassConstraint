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
#include "NCSplinePdfFactory_1D.h"
#include "NCSplinePdfFactory_2D.h"

using namespace RooFit;
using namespace std;

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
Bool_t checkProjVar(const vector<string>& strprojvars, RooRealVar*& var, vector<RooRealVar*>* addArray=0);
void addByLowest(double val, std::vector<double>& valArray);
string set_gHVV(RooSpinZero::modelCouplings& couplings, string gSet);
void getProjection_2D_single_HVV(string gSet, string decaytype, string strprojvar);

void getProjection_2D(){
  gSystem->Exec("mkdir -p ./plots");
  getProjection_2D_single_HVV("ghz1=1,0", "ZZ4l", "m1:m12");
}

void getProjection_2D_single_HVV(string gSet, string decaytype, string strprojvar){
  std::vector<RooRealVar*> projvars;
  vector<string> strprojvars;
  splitOptionRecursive(strprojvar, strprojvars, ':');

  TString outdir="";
  for (unsigned int is=0; is<strprojvars.size(); is++){
    if (is!=strprojvars.size()-1) outdir.Append(Form("%s_", strprojvars.at(is).c_str()));
    else outdir.Append(Form("%s", strprojvars.at(is).c_str()));
  }

  RooRealVar* m12 = new RooRealVar("m12", "m_{H} (GeV)", m1Range[0]+m2Range[0], 20000.);
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
  checkProjVar(strprojvars, m12, &projvars);
  checkProjVar(strprojvars, m1, &projvars);
  checkProjVar(strprojvars, m2, &projvars);
  int whichIsMH=-1;
  for (unsigned int ip=0; ip<projvars.size(); ip++){
    if (projvars.at(ip)==m12) whichIsMH = ip;
  }

  if (projvars.size()==2){
    gSystem->Exec(Form("mkdir -p ./plots/%s", outdir.Data()));

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

    cout << "Will compute mH=\n";
    vector<double> masses;
    double massmin=(int)(m1Range[0]+m2Range[0]+m1m2_BinWidth+0.5);
    double mass=massmin;
    while (mass<=15000){
      cout << mass << "\n";
      masses.push_back(mass);
      double massinc;
      if (mass<200.) massinc=1;
      else if (mass<600.) massinc=10.;
      else if (mass<1500.) massinc=50.;
      else if (mass<3000.) massinc=100.;
      else if (mass<10000.) massinc=500.;
      else massinc=1000.;
      mass += massinc;
    }
    const unsigned int nbins=masses.size();
    cout << "Total number of mass points to evaluate: " << nbins << endl;
    cout << endl;

    vector<vector<double>> projvals;
    for (unsigned int ip=0; ip<projvars.size(); ip++){
      unsigned int nevals = 101;
      double projmin = projvars.at(ip)->getMin();
      double projmax = projvars.at(ip)->getMax();
      double projwidth = (projmax-projmin)/((double)(nevals-1));
      if ((int)ip==whichIsMH){
        nevals = nbins;
        projmin = masses.at(0);
        projmax = masses.at(masses.size()-1);
        projwidth = 0;
      }

      if ((int)ip==whichIsMH) projvals.push_back(masses);
      else{
        vector<double> vals;
        for (unsigned int ie=0; ie<nevals; ie++){
          double val = projmin + ((double)ie)*projwidth;
          vals.push_back(val);
        }
        if (projvars.at(ip)==m1 || projvars.at(ip)==m2){
          double mv, gamv;
          pdf->getMVGamV(&mv, &gamv);
          cout << "Found mV, gamV = " << mv << " " << gamv << endl;
          double projnewmin=mv-2.*gamv;
          double projnewmax=mv+2.*gamv;
          unsigned int nextrabins=41;
          if (projvars.at(ip)==m2){
            projnewmin=mv-1.*gamv;
            projnewmax=mv+1.*gamv;
            nextrabins=81;
          }
          double projnewwidth=(projnewmax-projnewmin)/((const double)(nextrabins-1));
          cout << "Adding extra " << nextrabins << " points from " << projnewmin << " to " << projnewmax  << " in steps of " << projnewwidth << " to the variable " << projvars.at(ip)->GetName() << endl;
          for (unsigned int ip=0; ip<nextrabins; ip++){
            double projnewval=projnewmin+projnewwidth*ip;
            addByLowest(projnewval, vals);
          }
        }
        projvals.push_back(vals);
      }
      cout << "Projection variable " << projvars.at(ip)->GetName() << " [" << projvals.at(ip).at(0) << ", " << projvals.at(ip).at(projvals.at(ip).size()-1) << " : " << projvals.at(ip).size() << "]" << endl;
    }

    RooArgSet intSet;
    for (unsigned int ip=0; ip<projvars.size(); ip++){
      if (!checkProjVar(strprojvars, m1, 0)) intSet.add(*(measurables.m1));
      if (!checkProjVar(strprojvars, m2, 0)) intSet.add(*(measurables.m2));
      if (!checkProjVar(strprojvars, h1, 0) && measurables.h1!=0) intSet.add(*(measurables.h1));
      if (!checkProjVar(strprojvars, h2, 0) && measurables.h2!=0) intSet.add(*(measurables.h2));
      if (!checkProjVar(strprojvars, hs, 0) && measurables.hs!=0) intSet.add(*(measurables.hs));
      if (!checkProjVar(strprojvars, Phi, 0) && measurables.Phi!=0) intSet.add(*(measurables.Phi));
      if (!checkProjVar(strprojvars, Phi1, 0) && measurables.Phi1!=0) intSet.add(*(measurables.Phi1));
      if (!checkProjVar(strprojvars, Y, 0) && measurables.Y!=0) intSet.add(*(measurables.Y));
    }
    cout << "Integration variables:\n";
    intSet.Print("v");
    cout << endl;

    RooRealIntegral* pdf_int=0;
    if (intSet.getSize()>0) pdf_int = new RooRealIntegral(Form("pdf_proj_%s", outdir.Data()), "", *pdf, intSet);
    if (pdf_int!=0){
      delete pdf_int; // Delete dummy integration. Need to create at each mPOLE to optimize integration.

      double mPOLE = 125.;
      vector<doubleTriplet_t> points;
      vector<unsigned int> ndim; for (unsigned int iv=0; iv<projvals.size(); iv++) ndim.push_back(projvals.at(iv).size());

      double* dim = new double[projvars.size()];
      double fcn = 0;
      double spfcn0=0;
      double spfcn=0;

      TFile* foutput = TFile::Open(Form("H%sDecay_%sProjection_NoInterf_%s%s", decaytype.c_str(), outdir.Data(), strcoupl.c_str(), ".root"), "recreate");
      TTree* pointsTree = new TTree("points", "");
      for (unsigned int iv=0; iv<projvars.size(); iv++) pointsTree->Branch(Form("d%i", (int)iv+1), &(dim[iv]));
      pointsTree->Branch("fcn", &fcn);
      pointsTree->Branch("spfcn", &spfcn);
      pointsTree->Branch("spfcn0", &spfcn0);

      m12->setConstant(false);
      m12->setVal(mPOLE);
      m12->setConstant(true);

      pdf_int = new RooRealIntegral(Form("pdf_proj_%s", outdir.Data()), "", *pdf, intSet);
      for (unsigned int ix=0; ix<ndim[0]; ix++){
        dim[0]=projvals.at(0).at(ix);
        projvars.at(0)->setConstant(false);
        projvars.at(0)->setVal(dim[0]);
        projvars.at(0)->setConstant(true);
        for (unsigned int iy=0; iy<ndim[1]; iy++){
          dim[1]=projvals.at(1).at(iy);
          projvars.at(1)->setConstant(false);
          projvars.at(1)->setVal(dim[1]);
          projvars.at(1)->setConstant(true);
          fcn = pdf_int->getVal();
          points.push_back(doubleTriplet_t(dim[0], dim[1], fcn));
        }

        // Calculate Riemann integral for xcheck
        Double_t manualIntegral=0;
        for (unsigned int iy=0; iy<ndim[1]-1; iy++){
          unsigned int ip = ix*ndim[1]+iy;
          unsigned int ip_ypo = ip+1;
          manualIntegral += 0.5*(points.at(ip_ypo)[2]+points.at(ip)[2])*(points.at(ip_ypo)[1]-points.at(ip)[1]);
        }
        cout << "Original Riemann integral of PDF at " << projvars.at(0)->GetName() << "=" << projvals.at(0).at(ix) << ": " << manualIntegral << endl;
      }

      TString gifCmd;

      NCSplinePdfFactory_2D* spFactory = new NCSplinePdfFactory_2D(projvars.at(0), projvars.at(1));
      spFactory->setPoints(points);
      RooNCSplinePdf_2D_fast* spPDF = spFactory->getPDF();

      projvars.at(0)->setConstant(false);
      projvars.at(1)->setConstant(false);
      RooRealIntegral* spPDFint = 0;
      RooArgSet splineIntvars;
      if (whichIsMH>=0){
        for (unsigned int ip=0; ip<projvars.size(); ip++){ if ((int)ip!=whichIsMH) splineIntvars.add(*(projvars.at(ip))); }
        cout << "Spline integration vars except m4l:" << endl;
        splineIntvars.Print("v");
        spPDFint = new RooRealIntegral("spPDFint", "", *spPDF, splineIntvars);
      }

      for (unsigned int ix=0; ix<ndim[0]; ix++){
        Double_t manualIntegral=0;
        for (unsigned int iy=0; iy<ndim[1]; iy++){
          unsigned int ip = ix*ndim[1]+iy;
          unsigned int ip_ypo = ip; if (iy<(ndim[1]-1)) ip_ypo++;

          dim[0] = points.at(ip)[0];
          dim[1] = points.at(ip)[1];
          fcn = points.at(ip)[2];

          if (dim[0]==155.) cout << "Constructing the renormalized spline for point (x,y,fcn)=(" << dim[0] << ", " << dim[1] << ", " << fcn << ")" << endl;

          projvars.at(0)->setVal(dim[0]);
          projvars.at(1)->setVal(dim[1]);

          spfcn0 = spPDF->getVal();
          spfcn = spfcn0;
          if (spPDFint!=0){
            double integral = spPDFint->getVal();
            spfcn /= integral;
            points.at(ip)[2] /= integral;
            if (dim[0]==155.) cout << "\tIntegral: " << integral << endl;
          }
          pointsTree->Fill();
        }

        for (unsigned int iy=0; iy<ndim[1]-1; iy++){
          unsigned int ip = ix*ndim[1]+iy;
          unsigned int ip_ypo = ip+1;
          manualIntegral += 0.5*(points.at(ip_ypo)[2]+points.at(ip)[2])*(points.at(ip_ypo)[1]-points.at(ip)[1]);
        }
        cout << "Renormalized Riemann integral of (spline) PDF at " << projvars.at(0)->GetName() << "=" << projvals.at(0).at(ix) << ": " << manualIntegral << endl;
      }
      if (spPDFint!=0){
        delete spPDFint;
        spFactory->setPoints(points); // Reset spline points such that they are normalized along the m4l slices
        spPDF = spFactory->getPDF();
      }
      foutput->WriteTObject(pointsTree);
      delete pointsTree;

      if (whichIsMH>=0){
        mPOLE = 155.36;
        //mPOLE = 155;
        cout << "Plotting PDFs for central value " << mPOLE << endl;
        m12->setVal(mPOLE);
        m12->setConstant(true);

        spPDFint = new RooRealIntegral("spPDFint", "", *spPDF, splineIntvars);
        double spint = spPDFint->getVal();
        cout << "Spline integral at mH = " << mPOLE << ": " << spint << endl;
        delete spPDFint;

        RooPlot* plot = projvars.at(1-whichIsMH)->frame(projvals.at(1-whichIsMH).size()-1);
        plot->GetXaxis()->CenterTitle();
        plot->GetYaxis()->SetTitleOffset(1.2);
        plot->GetYaxis()->CenterTitle();
        plot->GetXaxis()->SetTitle(projvars.at(1-whichIsMH)->GetTitle());
        plot->GetYaxis()->SetTitle(Form("H%s amplitude", decaytype.c_str()));
        plot->GetXaxis()->SetNdivisions(-505);
        plot->SetTitle(Form("Projection at mH=%.2f GeV", mPOLE));
        //spPDF->setVerbosity(RooNCSplinePdfCore::kVerbose);
        spPDF->plotOn(plot, LineColor(kRed), LineWidth(2));
        pdf->plotOn(plot, LineColor(kBlack), LineWidth(2), LineStyle(2), Project(intSet));

        m12->setConstant(false);
        mPOLE = (int)mPOLE;
        m12->setVal(mPOLE);
        m12->setConstant(true);
        cout << "Plotting PDF for edge value " << mPOLE << endl;
        pdf->plotOn(plot, LineColor(kBlue), LineWidth(2), LineStyle(2), Project(intSet));
        m12->setConstant(false);
        mPOLE = mPOLE+1.;
        m12->setVal(mPOLE);
        m12->setConstant(true);
        cout << "Plotting PDF for edge value " << mPOLE << endl;
        pdf->plotOn(plot, LineColor(kGreen+2), LineWidth(2), LineStyle(2), Project(intSet));

        TCanvas* canvas = new TCanvas(Form("c_%s", outdir.Data()), "", 800, 800);
        plot->Draw();
        canvas->Modified();
        canvas->Update();
        foutput->WriteTObject(canvas);
        gifCmd += savePlot(outdir, canvas);
        gifCmd += " ";
        canvas->Close();
        delete plot;
      }
      delete spFactory;
      delete pdf_int;
      foutput->Close();
      /*
      if (gifCmd!=""){
        gifCmd.Prepend("convert ");
        gifCmd.Append(Form("./plots/%s/all.gif", projvar->GetName()));
        gSystem->Exec(gifCmd);
      }
      */
      delete[] dim;
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

Bool_t checkProjVar(const vector<string>& strprojvars, RooRealVar*& var, vector<RooRealVar*>* addArray){
  if (var==0) return false;
  Bool_t found = checkListVariable(strprojvars, var->GetName());
  if (found && addArray!=0){
    for (unsigned int ia=0; ia<addArray->size(); ia++){
      if (string(var->GetName())==string(addArray->at(ia)->GetName())) return found;
    }
    addArray->push_back(var);
  }
  return found;
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