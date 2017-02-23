#include <cassert>
#include "HiggsMassConstraint.h"
#include "RooFit.h"
#include "RooNLLVar.h"
#include "RooMinimizer.h"
#include "RooMinuit.h"
#include "RooNameReg.h"
#include "RooCmdConfig.h"
#include "RooGlobalFunc.h"
#include "RooAddition.h"
#include "RooNumIntConfig.h"
#include "RooProjectedPdf.h"
#include "RooInt.h"
#include "RooCustomizer.h"
#include "RooConstraintSum.h"
#include "RooCachedReal.h"
#include "RooXYChi2Var.h"
#include "RooChi2Var.h"
#include "RooRealIntegral.h"
#include "TIterator.h"
#include "Math/Minimizer.h"


using namespace std;
using namespace reco;
using namespace pat;
using namespace HMCtoolkit;


HiggsMassConstraint::HiggsMassConstraint(
  const Double_t sqrts_,
  RooSpin::VdecayType Vdecay1_,
  RooSpin::VdecayType Vdecay2_,
  const Int_t X_spin_,
  const Int_t intCodeStart_
  ) :
  sqrts(sqrts_*1000.),
  Vdecay1(Vdecay1_),
  Vdecay2(Vdecay2_),
  X_spin(X_spin_),
  intCodeStart(intCodeStart_),
  fastPDF(0)
{
  if (X_spin!=0 && X_spin!=2){
    cerr << "HiggsMassConstraint::HiggsMassConstraint: Only spin-0 or -2 are supported!" << endl;
    assert(0);
  }

  resetInitialMassErrors(); // Reset the double arrays

  setFastPDF(); // By default, this is set to false. If true, it will attempt to use the product of two BWs.

  setFitMomentumStrategy(); // Set default momentum strategy
  setFitVVStrategy(); // Set default VV strategy
  setJECUserFloatString(""); // Set JEC uncertainty default string
  setMuonKalmanCorrectedPtErrorString(""); // Set string to obtain a userFloat for the corrected pT error on muons after Kalman fit (+ smearing in MC)

  setPtEtaCuts(); // Set default cuts on pT, eta and phi of leptons, jets and FSR. Note the the cut targets are numbers!
  constructVariables();
  setM1M2Cuts(); // Set default minimum cuts on m1, m2, mA, mB. Note that the cut targets are RooRealVars, so constructVariables needs to be called first!
  constructPdfFactory();
  constructSplinePDF();
  constructConstraintPdfs();
  constructCompoundPdf();
}

HiggsMassConstraint::~HiggsMassConstraint(){
  destroyCompoundPdf();
  destroyConstraintPdfs();
  destroySplinePDF();
  destroyPdfFactory();
  destroyVariables();
}

RooAbsPdf* HiggsMassConstraint::getPDF(){ return PDF; }
SpinPdfFactory* HiggsMassConstraint::getPDFFactory(){ return pdfFactory; }

void HiggsMassConstraint::constructVariables(){
  fitResult=0;
  varZero = new RooRealVar("varZero", "", 0.);
  varOne = new RooRealVar("varOne", "", 1.);

  m1lowcut = new RooRealVar("m1lowcut", "", 0., 0., sqrts);
  m2lowcut = new RooRealVar("m2lowcut", "", 0., 0., sqrts);
  m1highcut = new RooRealVar("m1highcut", "", sqrts, 0., sqrts);
  m2highcut = new RooRealVar("m2highcut", "", sqrts, 0., sqrts);
  mFFOScut = new RooRealVar("mFFOScut", "", 0., 0., sqrts);
  mFFSScut = new RooRealVar("mFFSScut", "", 0., 0., sqrts);

  for (int iZ=0; iZ<2; iZ++) mManip[iZ] = new RooRealVar(Form("m%i_manip", iZ+1), "", 0., 0., sqrts); // Special variables for spline manipulation
  mManip[2] = new RooRealVar("m12_manip", "", 0., 0., sqrts);

  RooArgList m12_args;
  for (int iZ=0; iZ<2; iZ++){
    RooArgList mHdaughter_args;
    for (int iferm=0; iferm<2; iferm++){
      pT_ferm[iZ][iferm] = new RooRealVar(Form("pTRefit_Z%iFermion%i", iZ+1, iferm+1), "", 0., 0., sqrts);
      lambda_ferm[iZ][iferm] = new RooRealVar(Form("lambdaRefit_Z%iFermion%i", iZ+1, iferm+1), "", -piovertwo_val, piovertwo_val);
      phi_ferm[iZ][iferm] = new RooRealVar(Form("phiRefit_Z%iFermion%i", iZ+1, iferm+1), "", -pi_val, pi_val);
      massbar_ferm[iZ][iferm] = new RooRealVar(Form("massInit_Z%iFermion%i", iZ+1, iferm+1), "", 0., -sqrts, sqrts);

      pT_fsr[iZ][iferm] = new RooRealVar(Form("pTRefit_Z%iFermion%iFSR", iZ+1, iferm+1), "", 0., 0., sqrts);
      lambda_fsr[iZ][iferm] = new RooRealVar(Form("lambdaRefit_Z%iFermion%iFSR", iZ+1, iferm+1), "", -piovertwo_val, piovertwo_val);
      phi_fsr[iZ][iferm] = new RooRealVar(Form("phiRefit_Z%iFermion%iFSR", iZ+1, iferm+1), "", -pi_val, pi_val);

      pTobs_ferm[iZ][iferm] = new RooRealVar(Form("pTObs_Z%iFermion%i", iZ+1, iferm+1), "", 0., 0., sqrts);
      lambdaobs_ferm[iZ][iferm] = new RooRealVar(Form("lambdaObs_Z%iFermion%i", iZ+1, iferm+1), "", -piovertwo_val, piovertwo_val);
      phiobs_ferm[iZ][iferm] = new RooRealVar(Form("phiObs_Z%iFermion%i", iZ+1, iferm+1), "", -pi_val, pi_val);
      pTobs_fsr[iZ][iferm] = new RooRealVar(Form("pTObs_Z%iFermion%iFSR", iZ+1, iferm+1), "", 0., 0., sqrts);
      lambdaobs_fsr[iZ][iferm] = new RooRealVar(Form("lambdaObs_Z%iFermion%iFSR", iZ+1, iferm+1), "", -piovertwo_val, piovertwo_val);
      phiobs_fsr[iZ][iferm] = new RooRealVar(Form("phiObs_Z%iFermion%iFSR", iZ+1, iferm+1), "", -pi_val, pi_val);

      massbar_ferm[iZ][iferm]->removeMin();
      massbar_ferm[iZ][iferm]->removeMax();

      // Covariance elements
      for (int ix=0; ix<3; ix++){
        TString title_ix="";
        if (ix==0) title_ix="pT";
        else if (ix==1) title_ix="lambda";
        else title_ix="phi";
        for (int iy=0; iy<3; iy++){
          TString title_iy="";
          if (iy==0) title_iy="pT";
          else if (iy==1) title_iy="lambda";
          else title_iy="phi";
          Double_t minval=-1e9, maxval=1e9;
          if (ix==iy) minval=0;
          invcov_ferm[iZ][iferm][3*ix+iy] = new RooRealVar(Form("%s_vs_%s_Z%iFermion%i", title_ix.Data(), title_iy.Data(), iZ+1, iferm+1), "", 0., minval, maxval);
          if (ix!=iy) invcov_ferm[iZ][iferm][3*ix+iy]->removeMin();
          invcov_ferm[iZ][iferm][3*ix+iy]->removeMax();
          invcov_fsr[iZ][iferm][3*ix+iy] = new RooRealVar(Form("%s_vs_%s_Z%iFermion%iFSR", title_ix.Data(), title_iy.Data(), iZ+1, iferm+1), "", 0., minval, maxval);
          if (ix!=iy) invcov_fsr[iZ][iferm][3*ix+iy]->removeMin();
          invcov_fsr[iZ][iferm][3*ix+iy]->removeMax();
        }
      }

      E_ferm[iZ][iferm] = new RooFormulaVar(Form("ERefit_Z%iFermion%i", iZ+1, iferm+1), "sqrt( pow(@0,2)*TMath::Sign(1.,pow(@0,2)) + pow(@1/cos(@2),2) )", RooArgList(*(massbar_ferm[iZ][iferm]), *(pT_ferm[iZ][iferm]), *(lambda_ferm[iZ][iferm])));
      px_ferm[iZ][iferm] = new RooFormulaVar(Form("pxRefit_Z%iFermion%i", iZ+1, iferm+1), "(@0*cos(@1))", RooArgList(*(pT_ferm[iZ][iferm]), *(phi_ferm[iZ][iferm])));
      py_ferm[iZ][iferm] = new RooFormulaVar(Form("pyRefit_Z%iFermion%i", iZ+1, iferm+1), "(@0*sin(@1))", RooArgList(*(pT_ferm[iZ][iferm]), *(phi_ferm[iZ][iferm])));
      pz_ferm[iZ][iferm] = new RooFormulaVar(Form("pzRefit_Z%iFermion%i", iZ+1, iferm+1), "(@0*tan(@1))", RooArgList(*(pT_ferm[iZ][iferm]), *(lambda_ferm[iZ][iferm])));
      E_fsr[iZ][iferm] = new RooFormulaVar(Form("ERefit_Z%iFermion%iFSR", iZ+1, iferm+1), "(@0/cos(@1))", RooArgList(*(pT_fsr[iZ][iferm]), *(lambda_fsr[iZ][iferm])));
      px_fsr[iZ][iferm] = new RooFormulaVar(Form("pxRefit_Z%iFermion%iFSR", iZ+1, iferm+1), "(@0*cos(@1))", RooArgList(*(pT_fsr[iZ][iferm]), *(phi_fsr[iZ][iferm])));
      py_fsr[iZ][iferm] = new RooFormulaVar(Form("pyRefit_Z%iFermion%iFSR", iZ+1, iferm+1), "(@0*sin(@1))", RooArgList(*(pT_fsr[iZ][iferm]), *(phi_fsr[iZ][iferm])));
      pz_fsr[iZ][iferm] = new RooFormulaVar(Form("pzRefit_Z%iFermion%iFSR", iZ+1, iferm+1), "(@0*tan(@1))", RooArgList(*(pT_fsr[iZ][iferm]), *(lambda_fsr[iZ][iferm])));

      E_Hdaughter[iZ][iferm] = new RooFormulaVar(Form("ERefit_Z%iDau%i", iZ+1, iferm+1), "(@0+@1)", RooArgList(*(E_ferm[iZ][iferm]), *(E_fsr[iZ][iferm])));
      px_Hdaughter[iZ][iferm] = new RooFormulaVar(Form("pxRefit_Z%iDau%i", iZ+1, iferm+1), "(@0+@1)", RooArgList(*(px_ferm[iZ][iferm]), *(px_fsr[iZ][iferm])));
      py_Hdaughter[iZ][iferm] = new RooFormulaVar(Form("pyRefit_Z%iDau%i", iZ+1, iferm+1), "(@0+@1)", RooArgList(*(py_ferm[iZ][iferm]), *(py_fsr[iZ][iferm])));
      pz_Hdaughter[iZ][iferm] = new RooFormulaVar(Form("pzRefit_Z%iDau%i", iZ+1, iferm+1), "(@0+@1)", RooArgList(*(pz_ferm[iZ][iferm]), *(pz_fsr[iZ][iferm])));
      // There is no Sign needed for m_Hdaughter since it is only used in beta_Vdaughter, which does not care about the sign.
      m_Hdaughter[iZ][iferm] = new RooFormulaVar(Form("m_Z%iDau%iRefit", iZ+1, iferm+1), "sqrt( abs(pow(@0,2)-pow(@1,2)-pow(@2,2)-pow(@3,2)) )", RooArgList(*(E_Hdaughter[iZ][iferm]), *(px_Hdaughter[iZ][iferm]), *(py_Hdaughter[iZ][iferm]), *(pz_Hdaughter[iZ][iferm])));

      mHdaughter_args.add(*(E_Hdaughter[iZ][iferm]));
      mHdaughter_args.add(*(px_Hdaughter[iZ][iferm]));
      mHdaughter_args.add(*(py_Hdaughter[iZ][iferm]));
      mHdaughter_args.add(*(pz_Hdaughter[iZ][iferm]));
    }

    // Add mHdaughter arguments into m12 arguments
    m12_args.add(mHdaughter_args);

    // Construct m1/m2
    m[iZ] = new RooFormulaVar(Form("m%iRefit", iZ+1), "sqrt( TMath::Max(1e-15, pow(@0+@4,2)-pow(@1+@5,2)-pow(@2+@6,2)-pow(@3+@7,2)) )", mHdaughter_args);

    // This beta should multiply the spinPDF bc. having massive fermions have additional scale factors. These factors are even more relevant when FSR is present!
    beta_Vdaughter[iZ] = new RooFormulaVar(Form("betaV%iRefit", iZ+1), "(@0>0. ? sqrt( TMath::Max(1e-15, ( 1.-pow((@1+@2)/@0,2) )*( 1.-pow((@1-@2)/@0,2) ) ) ) : 1.)", RooArgList(*(m[iZ]), *(m_Hdaughter[iZ][0]), *(m_Hdaughter[iZ][1])));
  }

  // Construct m12
  m[2] = new RooFormulaVar("m12Refit", "sqrt( pow(@0+@4+@8+@12,2)-pow(@1+@5+@9+@13,2)-pow(@2+@6+@10+@14,2)-pow(@3+@7+@11+@15,2) )", m12_args);

  // Construct mA, mB
  RooArgList mA_args;
  RooArgList mB_args;
  RooArgList mAp_args;
  RooArgList mBp_args;
  for (int iZ=0; iZ<2; iZ++){
    int iferm = 0;
    if (iZ==0){
      mA_args.add(*(E_Hdaughter[iZ][iferm]));
      mA_args.add(*(px_Hdaughter[iZ][iferm]));
      mA_args.add(*(py_Hdaughter[iZ][iferm]));
      mA_args.add(*(pz_Hdaughter[iZ][iferm]));
    }
    else{
      mB_args.add(*(E_Hdaughter[iZ][iferm]));
      mB_args.add(*(px_Hdaughter[iZ][iferm]));
      mB_args.add(*(py_Hdaughter[iZ][iferm]));
      mB_args.add(*(pz_Hdaughter[iZ][iferm]));
    }
    mAp_args.add(*(E_Hdaughter[iZ][iferm]));
    mAp_args.add(*(px_Hdaughter[iZ][iferm]));
    mAp_args.add(*(py_Hdaughter[iZ][iferm]));
    mAp_args.add(*(pz_Hdaughter[iZ][iferm]));

    iferm = 1;
    if (iZ==0){
      mA_args.add(*(E_Hdaughter[1-iZ][iferm]));
      mA_args.add(*(px_Hdaughter[1-iZ][iferm]));
      mA_args.add(*(py_Hdaughter[1-iZ][iferm]));
      mA_args.add(*(pz_Hdaughter[1-iZ][iferm]));
    }
    else{
      mB_args.add(*(E_Hdaughter[1-iZ][iferm]));
      mB_args.add(*(px_Hdaughter[1-iZ][iferm]));
      mB_args.add(*(py_Hdaughter[1-iZ][iferm]));
      mB_args.add(*(pz_Hdaughter[1-iZ][iferm]));
    }
    mBp_args.add(*(E_Hdaughter[1-iZ][iferm]));
    mBp_args.add(*(px_Hdaughter[1-iZ][iferm]));
    mBp_args.add(*(py_Hdaughter[1-iZ][iferm]));
    mBp_args.add(*(pz_Hdaughter[1-iZ][iferm]));
  }
  mAB[0] = new RooFormulaVar("mAOSRefit", "sqrt( TMath::Max(1e-15, pow(@0+@4,2)-pow(@1+@5,2)-pow(@2+@6,2)-pow(@3+@7,2)) )", mA_args);
  mAB[1] = new RooFormulaVar("mBOSRefit", "sqrt( TMath::Max(1e-15, pow(@0+@4,2)-pow(@1+@5,2)-pow(@2+@6,2)-pow(@3+@7,2)) )", mB_args);
  mAB[2] = new RooFormulaVar("mASSRefit", "sqrt( TMath::Max(1e-15, pow(@0+@4,2)-pow(@1+@5,2)-pow(@2+@6,2)-pow(@3+@7,2)) )", mAp_args);
  mAB[3] = new RooFormulaVar("mBSSRefit", "sqrt( TMath::Max(1e-15, pow(@0+@4,2)-pow(@1+@5,2)-pow(@2+@6,2)-pow(@3+@7,2)) )", mBp_args);

  // Variables integrated over
  // Should be re-written in terms of pT, lambda and phi at some point...
  hs = new RooRealVar("Gencosthetastar", "cos#theta^{*}", -1, 1);
  h1 = new RooRealVar("GenhelcosthetaZ1", "cos#theta_{1}", -1, 1);
  h2 = new RooRealVar("GenhelcosthetaZ2", "cos#theta_{2}", -1, 1);
  Phi = new RooRealVar("Genhelphi", "#Phi", -pi_val, pi_val);
  Phi1 = new RooRealVar("GenphistarZ1", "#Phi_{1}", -pi_val, pi_val);
  Y = new RooRealVar("GenY", "Y", 0);

  // Initialize the meaurables: Set those always integrated over to 0
  measurables.m1 = m[0];
  measurables.m2 = m[1];
  measurables.m12 = m[2];
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
}
void HiggsMassConstraint::constructPdfFactory(){
  hvvFactory=0;
  hvvFastFactory=0;
  xvvFactory=0;
  if (X_spin==0){
    hvvFactory = new ScalarPdfFactory_HVV(measurables, false, Vdecay1, Vdecay2, true); // First false for acceptance, then true for always-on-shell H
    hvvFactory->makeParamsConst(false); // So that we can play with couplings
    hvvFactory->setZZ4fOrdering(false); // Disable m1/m2 ordering
    //spinPDF = hvvFactory->getPDF();
    //pdfFactory = hvvFactory;

#ifdef _hmcpkgpathstr_
    const string HMCPKGPATH = _hmcpkgpathstr_;
    setFastIntegrationGraph(Form("%s/data/HZZ4lDecay_mHstarShape_NoInterf_ghz1.root", HMCPKGPATH.c_str()), "tg_pdf_mH");
#else
    cout << "HiggsMassConstraint::constructPdfFactory: Package path is undefined! Please modify BuildFile.xml." << endl;
    assert(0);
#endif
    hvvFastFactory = new ScalarPdfFactory_HVV_fast(measurables, false, Vdecay1, Vdecay2, true); // First false for acceptance, then true for always-on-shell H
    hvvFastFactory->makeParamsConst(false); // So that we can play with couplings
    hvvFastFactory->setZZ4fOrdering(false); // Disable m1/m2 ordering
    hvvFastFactory->setFastPDFOverallIntegration(tgint);
    spinPDF = hvvFastFactory->getPDF();
    pdfFactory = hvvFastFactory;
  }
  else if (X_spin==2){
    xvvFactory = new TensorPdfFactory_ppHVV(measurables, Vdecay1, Vdecay2, true); // true for always-on-shell X
    xvvFactory->makeParamsConst(false); // So that we can play with couplings
    xvvFactory->setZZ4fOrdering(false); // Disable m1/m2 ordering
    spinPDF = xvvFactory->getPDF();
    pdfFactory = xvvFactory;
  }
  spinPDF->alwaysIntegrate(intCodeStart);

  bwProdPDF = new RooRelBWProduct(
    "RelBWProductPDF", "RelBWProductPDF",
    pdfFactory->measurables,
    pdfFactory->parameters,
    Vdecay1, Vdecay2
    );

  RooGenericPdf* simpleBWPDF1 = new RooGenericPdf(
    "SimpleBreitWignerPDF1", "2*@0/( pow(pow(@0, 2)-pow(@1, 2), 2)+pow(@1*@2, 2) )",
    RooArgSet(*(pdfFactory->measurables.m1), *(pdfFactory->parameters.mZ), *(pdfFactory->parameters.gamZ))
    );
  RooGenericPdf* simpleBWPDF2 = new RooGenericPdf(
    "SimpleBreitWignerPDF2", "2*@0/( pow(pow(@0, 2)-pow(@1, 2), 2)+pow(@1*@2, 2) )",
    RooArgSet(*(pdfFactory->measurables.m2), *(pdfFactory->parameters.mZ), *(pdfFactory->parameters.gamZ))
    );
  simpleBWPDF.push_back(simpleBWPDF1);
  simpleBWPDF.push_back(simpleBWPDF2);

#if hmc_debug==1
  cout << "spinPDF: " << endl;
  spinPDF->Print("v");
  RooArgSet* pdfPars = spinPDF->getParameters((RooArgSet*)0, true);
  TIterator* parIter = pdfPars->createIterator();
  RooAbsArg* thePar;
  while ((thePar = (RooAbsArg*)parIter->Next())){
    cout << thePar->GetName() << endl;
  }
  cout << endl;
  delete parIter;
  delete pdfPars;

  cout << "bwProdPDF: " << endl;
  bwProdPDF->Print("v");
  pdfPars = bwProdPDF->getParameters((RooArgSet*)0, true);
  parIter = pdfPars->createIterator();
  while ((thePar = (RooAbsArg*)parIter->Next())){
    cout << thePar->GetName() << endl;
  }
  cout << endl;
  delete parIter;
  delete pdfPars;

  cout << "SimpleBWPDF: " << endl;
  for (unsigned int ip=0; ip<simpleBWPDF.size(); ip++){
    simpleBWPDF.at(ip)->Print("v");
    pdfPars = simpleBWPDF.at(ip)->getParameters((RooArgSet*)0, true);
    parIter = pdfPars->createIterator();
    while ((thePar = (RooAbsArg*)parIter->Next())){
      cout << thePar->GetName() << endl;
    }
    cout << endl;
    delete parIter;
    delete pdfPars;
  }
#endif
}
void HiggsMassConstraint::constructSplinePDF(){
#ifdef _hmcpkgpathstr_
  const string HMCPKGPATH = _hmcpkgpathstr_;
  TFile* fin = TFile::Open(Form("%s/data/HZZ4lDecay_m12_m1Projection_NoInterf_ghz1.root", HMCPKGPATH.c_str()), "read");
  TTree* tin = (TTree*)fin->Get("points");
  const unsigned  int ndim=2;
  vector<doubleTriplet_t> points;
  double fcn;
  double dim[ndim]={ 0 };
  tin->SetBranchAddress("spfcn", &fcn);
  for (unsigned int idim=0; idim<ndim; idim++) tin->SetBranchAddress(Form("d%i", idim+1), dim+idim);
  double massmax=0, massmin=1e5;
  for (int ev=0; ev<tin->GetEntries(); ev++){
    tin->GetEntry(ev);
    points.push_back(doubleTriplet_t(dim[0], dim[1], fcn));
    massmin = min(dim[0], massmin);
    massmax = max(dim[0], massmax);
    for (unsigned int ic=0; ic<2; ic++) addUnique(dim[ic], splineCoord[ic]);
  }

  // Construct 2D spline factory and its PDF
  spline2DFactory = new NCSplinePdfFactory_2D(mManip[2], mManip[0], "2D");
  spline2DFactory->setPoints(points);
  //((RooNCSplinePdfCore*)spline2DFactory->getPDF())->setVerbosity(RooNCSplinePdfCore::kVerbose);
  sqrts=min(massmax, sqrts);
  mManip[2]->setRange(massmin, sqrts);
  mManip[0]->setRange(splineCoord[1].at(0), splineCoord[1].at(splineCoord[1].size()-1));
  fin->Close();

  // Construct 1D spline factory but not its PDF
  spline1DFactory = new NCSplinePdfFactory_1D(m[0], "1D");
#else
  cout << "HiggsMassConstraint::constructSplinePDFs: Package path is undefined! Please modify BuildFile.xml." << endl;
  assert(0);
#endif
}
void HiggsMassConstraint::constructConstraintPdfs(){
  RooArgList constraints;
  for (int iZ=0; iZ<2; iZ++){
    for (int iferm=0; iferm<2; iferm++){
      RooArgList vars_ferm, means_ferm, me_ferm;
      vars_ferm.add(*(pT_ferm[iZ][iferm]));
      vars_ferm.add(*(lambda_ferm[iZ][iferm]));
      vars_ferm.add(*(phi_ferm[iZ][iferm]));
      means_ferm.add(*(pTobs_ferm[iZ][iferm]));
      means_ferm.add(*(lambdaobs_ferm[iZ][iferm]));
      means_ferm.add(*(phiobs_ferm[iZ][iferm]));
      for (int im=0; im<9; im++) me_ferm.add(*(invcov_ferm[iZ][iferm][im]));
      gausConstraintsPDF[iZ][iferm][0] = new RooGaussianMomConstraint(Form("gausConstraintsPDF_Z%iFermion%i", iZ+1, iferm+1), Form("gausConstraintsPDF_Z%iFermion%i", iZ+1, iferm+1), vars_ferm, means_ferm, me_ferm, RooGaussianMomConstraint::kRhoLambdaPhi);
      if (iZ==0) constraints.add(*(gausConstraintsPDF[iZ][iferm][0]));
      //gausConstraintsPDF[iZ][iferm][0]->setVerbosity(RooGaussianMomConstraint::kVerbose);

      RooArgList vars_fsr, means_fsr, me_fsr;
      vars_fsr.add(*(pT_fsr[iZ][iferm]));
      vars_fsr.add(*(lambda_fsr[iZ][iferm]));
      vars_fsr.add(*(phi_fsr[iZ][iferm]));
      means_fsr.add(*(pTobs_fsr[iZ][iferm]));
      means_fsr.add(*(lambdaobs_fsr[iZ][iferm]));
      means_fsr.add(*(phiobs_fsr[iZ][iferm]));
      for (int im=0; im<9; im++) me_fsr.add(*(invcov_fsr[iZ][iferm][im]));
      gausConstraintsPDF[iZ][iferm][1] = new RooGaussianMomConstraint(Form("gausConstraintsPDF_Z%iFermion%iFSR", iZ+1, iferm+1), Form("gausConstraintsPDF_Z%iFermion%iFSR", iZ+1, iferm+1), vars_fsr, means_fsr, me_fsr, RooGaussianMomConstraint::kRhoLambdaPhi);
      if (iZ==0) constraints.add(*(gausConstraintsPDF[iZ][iferm][1]));
      //gausConstraintsPDF[iZ][iferm][1]->setVerbosity(RooGaussianMomConstraint::kVerbose);
    }
  }

  RooArgList massCuts_args;
  massCuts_args.add(*(mAB[0])); // @0
  massCuts_args.add(*(mAB[1])); // @1
  massCuts_args.add(*mFFOScut); // @2
  massCuts_args.add(*(mAB[2])); // @3
  massCuts_args.add(*(mAB[3])); // @4
  massCuts_args.add(*mFFSScut); // @5
  massCuts_args.add(*(m[0])); // @6
  massCuts_args.add(*m1lowcut); // @7
  massCuts_args.add(*m1highcut); // @8
  massCuts_args.add(*(m[1])); // @9
  massCuts_args.add(*m2lowcut); // @10
  massCuts_args.add(*m2highcut); // @11
  massCuts = new RooFormulaVar("mABCutParameterization", "( (@0>=@2 && @1>=@2  &&  @3>=@5 && @4>=@5  &&  @6>=@7 && @6<=@8  &&  @9>=@10 && @9<=@11) ? 1. : 1.e-15)", massCuts_args);
  auxilliaryConstraintsPDF = new RooGenericPdf("auxilliaryConstraintsPDF", "@0*@1*@2", RooArgList(*(beta_Vdaughter[0]), *(beta_Vdaughter[1]), *massCuts));

  //constraints.add(*(auxilliaryConstraintsPDF));
  //constraints.add(*(DiracDeltaPDF));
  constraintsPDF = new RooProdPdf("constraintsPDF", "constraintsPDF", constraints);
#if hmc_debug==1
  RooArgSet* pdfPars;
  TIterator* parIter;
  RooAbsArg* thePar;

  constraintsPDF->Print("v");
  cout << "HiggsMassConstraint::constructConstraintPdfs: constraintsPDF parameters:" << endl;
  pdfPars = constraintsPDF->getParameters((RooArgSet*)0, true);
  parIter = pdfPars->createIterator();
  while ((thePar = (RooAbsArg*)parIter->Next())) cout << thePar->GetName() << endl;
  cout << endl;
  delete parIter;
  delete pdfPars;

  cout << "HiggsMassConstraint::constructConstraintPdfs: constraintsPDF observables:" << endl;
  pdfPars = constraintsPDF->getObservables((RooArgSet*)0, true);
  parIter = pdfPars->createIterator();
  while ((thePar = (RooAbsArg*)parIter->Next())) cout << thePar->GetName() << endl;
  cout << endl;
  delete parIter;
  delete pdfPars;

  cout << "HiggsMassConstraint::constructConstraintPdfs: constraintsPDF dependents:" << endl;
  pdfPars = constraintsPDF->getDependents((RooArgSet*)0);
  parIter = pdfPars->createIterator();
  while ((thePar = (RooAbsArg*)parIter->Next())) cout << thePar->GetName() << endl;
  cout << endl;
  delete parIter;
  delete pdfPars;
#endif
}
void HiggsMassConstraint::constructCompoundPdf(){
  RooArgList pdfList(*spinPDF, *constraintsPDF);
  PDF = new RooProdPdf("HiggsMassConstraint_PDF", "HiggsMassConstraint_PDF", pdfList);

#if hmc_debug==1
  RooArgSet* pdfPars;
  TIterator* parIter;
  RooAbsArg* thePar;

  PDF->Print("v");
  cout << "HiggsMassConstraint::constructCompoundPdf: PDF parameters:" << endl;
  pdfPars = PDF->getParameters((RooArgSet*)0, true);
  parIter = pdfPars->createIterator();
  while ((thePar = (RooAbsArg*)parIter->Next())) cout << thePar->GetName() << endl;
  cout << endl;
  delete parIter;
  delete pdfPars;

  cout << "HiggsMassConstraint::constructCompoundPdf: PDF observables:" << endl;
  pdfPars = PDF->getObservables((RooArgSet*)0, true);
  parIter = pdfPars->createIterator();
  while ((thePar = (RooAbsArg*)parIter->Next())) cout << thePar->GetName() << endl;
  cout << endl;
  delete parIter;
  delete pdfPars;

  cout << "HiggsMassConstraint::constructCompoundPdf: PDF dependents:" << endl;
  pdfPars = PDF->getDependents((RooArgSet*)0);
  parIter = pdfPars->createIterator();
  while ((thePar = (RooAbsArg*)parIter->Next())) cout << thePar->GetName() << endl;
  cout << endl;
  delete parIter;
  delete pdfPars;
#endif
}
void HiggsMassConstraint::constructCompoundFastPdf(){
  vector<Double_t> fcnList;
  for (unsigned int ip=0; ip<splineCoord[1].size(); ip++){
    mManip[0]->setVal(splineCoord[1].at(ip));
    Double_t fcn = spline2DFactory->getPDF()->getVal();
#if hmc_debug==1
    cout << "spline1DPDF(m1=" << mManip[0]->getVal() << ") = " << fcn << endl;
#endif
    fcnList.push_back(fcn);
  }
  spline1DFactory->setPoints(splineCoord[1], fcnList);
  RooNCSplinePdf_1D_fast* splinePDF = (RooNCSplinePdf_1D_fast*)spline1DFactory->getPDF();
  if (splinePDF!=0){
    RooArgList fastpdfList;
    fastpdfList.add(*splinePDF);
    fastpdfList.add(*constraintsPDF);
    fastPDF = new RooProdPdf("HiggsMassConstraint_FastPDF", "HiggsMassConstraint_FastPDF", fastpdfList);
  }
  else{
    cerr << "HiggsMassConstraint::constructCompoundFastPdf: ERROR: 1D spline could not be constructed. Aborting..." << endl;
    assert(0);
  }
#if hmc_debug==1
  RooArgSet* pdfPars;
  TIterator* parIter;
  RooAbsArg* thePar;

  fastPDF->Print("v");
  cout << "HiggsMassConstraint::constructCompoundFastPdf: fastPDF parameters:" << endl;
  pdfPars = fastPDF->getParameters((RooArgSet*)0, true);
  parIter = pdfPars->createIterator();
  while ((thePar = (RooAbsArg*)parIter->Next())) cout << thePar->GetName() << endl;
  cout << endl;
  delete parIter;
  delete pdfPars;

  cout << "HiggsMassConstraint::constructCompoundFastPdf: fastPDF observables:" << endl;
  pdfPars = fastPDF->getObservables((RooArgSet*)0, true);
  parIter = pdfPars->createIterator();
  while ((thePar = (RooAbsArg*)parIter->Next())) cout << thePar->GetName() << endl;
  cout << endl;
  delete parIter;
  delete pdfPars;

  cout << "HiggsMassConstraint::constructCompoundFastPdf: fastPDF dependents:" << endl;
  pdfPars = fastPDF->getDependents((RooArgSet*)0);
  parIter = pdfPars->createIterator();
  while ((thePar = (RooAbsArg*)parIter->Next())) cout << thePar->GetName() << endl;
  cout << endl;
  delete parIter;
  delete pdfPars;
#endif
}

void HiggsMassConstraint::destroyVariables(){
  deletePtr(fitResult);

  // Destroy in reverse order of creation
  deletePtr(h1);
  deletePtr(h2);
  deletePtr(hs);
  deletePtr(Phi);
  deletePtr(Phi1);
  deletePtr(Y);

  for (int iZ=3; iZ>=0; iZ--) deletePtr(mAB[iZ]);

  deletePtr(m[2]);
  for (int iZ=1; iZ>=0; iZ--){
    deletePtr(beta_Vdaughter[iZ]);
    deletePtr(m[iZ]);

    for (int iferm=1; iferm>=0; iferm--){
      deletePtr(E_Hdaughter[iZ][iferm]);
      deletePtr(px_Hdaughter[iZ][iferm]);
      deletePtr(py_Hdaughter[iZ][iferm]);
      deletePtr(pz_Hdaughter[iZ][iferm]);
      deletePtr(m_Hdaughter[iZ][iferm]);

      deletePtr(E_fsr[iZ][iferm]);
      deletePtr(px_fsr[iZ][iferm]);
      deletePtr(py_fsr[iZ][iferm]);
      deletePtr(pz_fsr[iZ][iferm]);
      deletePtr(E_ferm[iZ][iferm]);
      deletePtr(px_ferm[iZ][iferm]);
      deletePtr(py_ferm[iZ][iferm]);
      deletePtr(pz_ferm[iZ][iferm]);

      for (int ix=0; ix<3; ix++){
        for (int iy=0; iy<3; iy++){
          deletePtr(invcov_ferm[iZ][iferm][3*ix+iy]);
          deletePtr(invcov_fsr[iZ][iferm][3*ix+iy]);
        }
      }
      
      deletePtr(pTobs_ferm[iZ][iferm]);
      deletePtr(lambdaobs_ferm[iZ][iferm]);
      deletePtr(phiobs_ferm[iZ][iferm]);
      deletePtr(pTobs_fsr[iZ][iferm]);
      deletePtr(lambdaobs_fsr[iZ][iferm]);
      deletePtr(phiobs_fsr[iZ][iferm]);

      deletePtr(pT_fsr[iZ][iferm]);
      deletePtr(lambda_fsr[iZ][iferm]);
      deletePtr(phi_fsr[iZ][iferm]);

      deletePtr(pT_ferm[iZ][iferm]);
      deletePtr(lambda_ferm[iZ][iferm]);
      deletePtr(phi_ferm[iZ][iferm]);
      deletePtr(massbar_ferm[iZ][iferm]);
    }
  }

  for (unsigned int iv=0; iv<3; iv++) deletePtr(mManip[iv]);

  deletePtr(m1highcut);
  deletePtr(m2highcut);
  deletePtr(m1lowcut);
  deletePtr(m2lowcut);
  deletePtr(mFFOScut);
  deletePtr(mFFSScut);

  deletePtr(varOne);
  deletePtr(varZero);
}
void HiggsMassConstraint::destroySplinePDF(){ deletePtr(spline1DFactory); deletePtr(spline2DFactory); }
void HiggsMassConstraint::destroyPdfFactory(){
  // Delete the bwProdPDF and simpleBWPDF first since the measurables and parameters come from the pdfFactory
  for (unsigned int ip=0; ip<simpleBWPDF.size(); ip++) deletePtr(simpleBWPDF.at(ip));
  deletePtr(bwProdPDF);

  // Only one of these is true; no need to delete pdfFactory since it is simply a mother-pointer to either of these.
  pdfFactory=0;
  deletePtr(hvvFastFactory); deletePtr(tgint);
  deletePtr(hvvFactory);
  deletePtr(xvvFactory);
}
void HiggsMassConstraint::destroyConstraintPdfs(){
  deletePtr(constraintsPDF);
  deletePtr(auxilliaryConstraintsPDF);
  deletePtr(massCuts);
  for (int iZ=1; iZ>=0; iZ--){
    for (int iferm=1; iferm>=0; iferm--){
      for (int ifsr=1; ifsr>=0; ifsr--) deletePtr(gausConstraintsPDF[iZ][iferm][ifsr]);
    }
  }
}
void HiggsMassConstraint::destroyCompoundPdf(){ deletePtr(PDF); }
void HiggsMassConstraint::destroyCompoundFastPdf(){ deletePtr(fastPDF); }

void HiggsMassConstraint::setFastPDF(bool useFastPDF_){ useFastPDF = useFastPDF_; }
void HiggsMassConstraint::setPtEtaCuts(
  Double_t pTcut_muon_,
  Double_t etacut_muon_,
  Double_t pTcut_electron_,
  Double_t etacut_electron_,
  Double_t pTcut_jet_,
  Double_t etacut_jet_,
  Double_t pTcut_fsr_,
  Double_t etacut_fsr_
  ){
  pTcut_muon=pTcut_muon_;
  pTcut_electron=pTcut_electron_;
  pTcut_jet=pTcut_jet_;
  pTcut_fsr=pTcut_fsr_;

  if (etacut_muon_>0.) lambdacut_muon=piovertwo_val-2.*atan(exp(-etacut_muon_));
  else lambdacut_muon=piovertwo_val;
  if (etacut_electron_>0.) lambdacut_electron=piovertwo_val-2.*atan(exp(-etacut_electron_));
  else lambdacut_electron=piovertwo_val;
  if (etacut_jet_>0.) lambdacut_jet=piovertwo_val-2.*atan(exp(-etacut_jet_));
  else lambdacut_jet=piovertwo_val;
  if (etacut_fsr_>0.) lambdacut_fsr=piovertwo_val-2.*atan(exp(-etacut_fsr_));
  else lambdacut_fsr=piovertwo_val;
  // Actual setting of ranges is done per-event
}
void HiggsMassConstraint::setM1M2Cuts(
  Double_t m1lowcut_,
  Double_t m2lowcut_,
  Double_t m1highcut_,
  Double_t m2highcut_,
  Double_t mFFOScut_,
  Double_t mFFSScut_
  ){
  m1lowcut->setConstant(false);
  m2lowcut->setConstant(false);
  m1highcut->setConstant(false);
  m2highcut->setConstant(false);
  mFFOScut->setConstant(false);
  mFFSScut->setConstant(false);
  if (m1lowcut_>=0.) m1lowcut->setVal(m1lowcut_); else m1lowcut->setVal(0.);
  if (m2lowcut_>=0.) m2lowcut->setVal(m2lowcut_); else m2lowcut->setVal(0.);
  if (m1highcut_>=0.) m1highcut->setVal(m1highcut_); else m1highcut->setVal(sqrts);
  if (m2highcut_>=0.) m2highcut->setVal(m2highcut_); else m2highcut->setVal(sqrts);
  if (mFFOScut_>=0.) mFFOScut->setVal(mFFOScut_); else mFFOScut->setVal(0.);
  if (mFFSScut_>=0.) mFFSScut->setVal(mFFSScut_); else mFFSScut->setVal(0.);
  m1lowcut->setConstant(true);
  m2lowcut->setConstant(true);
  m1highcut->setConstant(true);
  m2highcut->setConstant(true);
  mFFOScut->setConstant(true);
  mFFSScut->setConstant(true);
}
void HiggsMassConstraint::setFastIntegrationGraph(TString strfname, TString strtgname){
  TFile* fin = TFile::Open(strfname, "read");
  tgint = (TGraph*)fin->Get(strtgname);
  fin->Close();
}


HiggsMassConstraint::FitMomentumStrategy HiggsMassConstraint::getFitMomentumStrategy(){ return fitMomStrategy_final; }
void HiggsMassConstraint::setWorkingFitMomentumStrategy(HiggsMassConstraint::FitMomentumStrategy fitMomStrategy_){ fitMomStrategy_final=fitMomStrategy_; }
void HiggsMassConstraint::setFitMomentumStrategy(HiggsMassConstraint::FitMomentumStrategy fitMomStrategy_){ fitMomStrategy=fitMomStrategy_; setWorkingFitMomentumStrategy(fitMomStrategy_); }
void HiggsMassConstraint::testFitMomentumStrategy(Int_t& useFullCov, Int_t& FermFSRType, Int_t& fitpT, Int_t& fitlambda, Int_t& fitphi) const{
  if (
    fitMomStrategy_final==HiggsMassConstraint::FullCov_All_pTLambdaPhi ||
    fitMomStrategy_final==HiggsMassConstraint::FullCov_All_pTLambda ||
    fitMomStrategy_final==HiggsMassConstraint::FullCov_All_pTPhi ||
    fitMomStrategy_final==HiggsMassConstraint::FullCov_All_LambdaPhi ||
    fitMomStrategy_final==HiggsMassConstraint::FullCov_NoFSR_pTLambdaPhi ||
    fitMomStrategy_final==HiggsMassConstraint::FullCov_NoFSR_pTLambda ||
    fitMomStrategy_final==HiggsMassConstraint::FullCov_NoFSR_pTPhi ||
    fitMomStrategy_final==HiggsMassConstraint::FullCov_NoFSR_LambdaPhi ||
    fitMomStrategy_final==HiggsMassConstraint::FullCov_OnlyFSR_pTLambdaPhi ||
    fitMomStrategy_final==HiggsMassConstraint::FullCov_OnlyFSR_pTLambda ||
    fitMomStrategy_final==HiggsMassConstraint::FullCov_OnlyFSR_pTPhi ||
    fitMomStrategy_final==HiggsMassConstraint::FullCov_OnlyFSR_LambdaPhi
    ) useFullCov=1;
  else useFullCov=0;

  if (
    fitMomStrategy_final==HiggsMassConstraint::FullCov_All_pTLambdaPhi ||
    fitMomStrategy_final==HiggsMassConstraint::FullCov_All_pTLambda ||
    fitMomStrategy_final==HiggsMassConstraint::FullCov_All_pTPhi ||
    fitMomStrategy_final==HiggsMassConstraint::FullCov_NoFSR_pTLambdaPhi ||
    fitMomStrategy_final==HiggsMassConstraint::FullCov_NoFSR_pTLambda ||
    fitMomStrategy_final==HiggsMassConstraint::FullCov_NoFSR_pTPhi ||
    fitMomStrategy_final==HiggsMassConstraint::FullCov_OnlyFSR_pTLambdaPhi ||
    fitMomStrategy_final==HiggsMassConstraint::FullCov_OnlyFSR_pTLambda ||
    fitMomStrategy_final==HiggsMassConstraint::FullCov_OnlyFSR_pTPhi ||
    fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_All_pTLambdaPhi ||
    fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_All_pTLambda ||
    fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_All_pTPhi ||
    fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_All_pT ||
    fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_NoFSR_pTLambdaPhi ||
    fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_NoFSR_pTLambda ||
    fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_NoFSR_pTPhi ||
    fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_NoFSR_pT ||
    fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_OnlyFSR_pTLambdaPhi ||
    fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_OnlyFSR_pTLambda ||
    fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_OnlyFSR_pTPhi ||
    fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_OnlyFSR_pT
    ) fitpT=1;
  else fitpT=0;

  if (
    fitMomStrategy_final==HiggsMassConstraint::FullCov_All_pTLambdaPhi ||
    fitMomStrategy_final==HiggsMassConstraint::FullCov_All_pTLambda ||
    fitMomStrategy_final==HiggsMassConstraint::FullCov_All_LambdaPhi ||
    fitMomStrategy_final==HiggsMassConstraint::FullCov_NoFSR_pTLambdaPhi ||
    fitMomStrategy_final==HiggsMassConstraint::FullCov_NoFSR_pTLambda ||
    fitMomStrategy_final==HiggsMassConstraint::FullCov_NoFSR_LambdaPhi ||
    fitMomStrategy_final==HiggsMassConstraint::FullCov_OnlyFSR_pTLambdaPhi ||
    fitMomStrategy_final==HiggsMassConstraint::FullCov_OnlyFSR_pTLambda ||
    fitMomStrategy_final==HiggsMassConstraint::FullCov_OnlyFSR_LambdaPhi ||
    fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_All_pTLambdaPhi ||
    fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_All_pTLambda ||
    fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_All_LambdaPhi ||
    fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_All_Lambda ||
    fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_NoFSR_pTLambdaPhi ||
    fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_NoFSR_pTLambda ||
    fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_NoFSR_LambdaPhi ||
    fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_NoFSR_Lambda ||
    fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_OnlyFSR_pTLambdaPhi ||
    fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_OnlyFSR_pTLambda ||
    fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_OnlyFSR_LambdaPhi ||
    fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_OnlyFSR_Lambda
    ) fitlambda=1;
  else fitlambda=0;

  if (
    fitMomStrategy_final==HiggsMassConstraint::FullCov_All_pTLambdaPhi ||
    fitMomStrategy_final==HiggsMassConstraint::FullCov_All_pTPhi ||
    fitMomStrategy_final==HiggsMassConstraint::FullCov_All_LambdaPhi ||
    fitMomStrategy_final==HiggsMassConstraint::FullCov_NoFSR_pTLambdaPhi ||
    fitMomStrategy_final==HiggsMassConstraint::FullCov_NoFSR_pTPhi ||
    fitMomStrategy_final==HiggsMassConstraint::FullCov_NoFSR_LambdaPhi ||
    fitMomStrategy_final==HiggsMassConstraint::FullCov_OnlyFSR_pTLambdaPhi ||
    fitMomStrategy_final==HiggsMassConstraint::FullCov_OnlyFSR_pTPhi ||
    fitMomStrategy_final==HiggsMassConstraint::FullCov_OnlyFSR_LambdaPhi ||
    fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_All_pTLambdaPhi ||
    fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_All_pTPhi ||
    fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_All_LambdaPhi ||
    fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_All_Phi ||
    fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_NoFSR_pTLambdaPhi ||
    fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_NoFSR_pTPhi ||
    fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_NoFSR_LambdaPhi ||
    fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_NoFSR_Phi ||
    fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_OnlyFSR_pTLambdaPhi ||
    fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_OnlyFSR_pTPhi ||
    fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_OnlyFSR_LambdaPhi ||
    fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_OnlyFSR_Phi
    ) fitphi=1;
  else fitphi=0;

  if (
    fitMomStrategy_final==HiggsMassConstraint::FullCov_All_pTLambdaPhi ||
    fitMomStrategy_final==HiggsMassConstraint::FullCov_All_pTLambda ||
    fitMomStrategy_final==HiggsMassConstraint::FullCov_All_pTPhi ||
    fitMomStrategy_final==HiggsMassConstraint::FullCov_All_LambdaPhi ||
    fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_All_pTLambdaPhi ||
    fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_All_pTLambda ||
    fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_All_pTPhi ||
    fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_All_LambdaPhi ||
    fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_All_pT ||
    fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_All_Lambda ||
    fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_All_Phi
    ) FermFSRType=2;
  else if (
    fitMomStrategy_final==HiggsMassConstraint::FullCov_OnlyFSR_pTLambdaPhi ||
    fitMomStrategy_final==HiggsMassConstraint::FullCov_OnlyFSR_pTLambda ||
    fitMomStrategy_final==HiggsMassConstraint::FullCov_OnlyFSR_pTPhi ||
    fitMomStrategy_final==HiggsMassConstraint::FullCov_OnlyFSR_LambdaPhi ||
    fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_OnlyFSR_pTLambdaPhi ||
    fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_OnlyFSR_pTLambda ||
    fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_OnlyFSR_pTPhi ||
    fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_OnlyFSR_LambdaPhi ||
    fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_OnlyFSR_pT ||
    fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_OnlyFSR_Lambda ||
    fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_OnlyFSR_Phi
    ) FermFSRType=1;
  else FermFSRType=0;

  /*
  fitMomStrategy_final==HiggsMassConstraint::FullCov_All_pTLambdaPhi ||
  fitMomStrategy_final==HiggsMassConstraint::FullCov_All_pTLambda ||
  fitMomStrategy_final==HiggsMassConstraint::FullCov_All_pTPhi ||
  fitMomStrategy_final==HiggsMassConstraint::FullCov_All_LambdaPhi ||
  fitMomStrategy_final==HiggsMassConstraint::FullCov_NoFSR_pTLambdaPhi ||
  fitMomStrategy_final==HiggsMassConstraint::FullCov_NoFSR_pTLambda ||
  fitMomStrategy_final==HiggsMassConstraint::FullCov_NoFSR_pTPhi ||
  fitMomStrategy_final==HiggsMassConstraint::FullCov_NoFSR_LambdaPhi ||
  fitMomStrategy_final==HiggsMassConstraint::FullCov_OnlyFSR_pTLambdaPhi ||
  fitMomStrategy_final==HiggsMassConstraint::FullCov_OnlyFSR_pTLambda ||
  fitMomStrategy_final==HiggsMassConstraint::FullCov_OnlyFSR_pTPhi ||
  fitMomStrategy_final==HiggsMassConstraint::FullCov_OnlyFSR_LambdaPhi ||
  fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_All_pTLambdaPhi ||
  fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_All_pTLambda ||
  fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_All_pTPhi ||
  fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_All_LambdaPhi ||
  fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_All_pT ||
  fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_All_Lambda ||
  fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_All_Phi ||
  fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_NoFSR_pTLambdaPhi ||
  fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_NoFSR_pTLambda ||
  fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_NoFSR_pTPhi ||
  fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_NoFSR_LambdaPhi ||
  fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_NoFSR_pT ||
  fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_NoFSR_Lambda ||
  fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_NoFSR_Phi ||
  fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_OnlyFSR_pTLambdaPhi ||
  fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_OnlyFSR_pTLambda ||
  fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_OnlyFSR_pTPhi ||
  fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_OnlyFSR_LambdaPhi ||
  fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_OnlyFSR_pT ||
  fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_OnlyFSR_Lambda ||
  fitMomStrategy_final==HiggsMassConstraint::CovDiagonals_OnlyFSR_Phi
  */
}
void HiggsMassConstraint::decrementMomentumStrategy(HiggsMassConstraint::FitMomentumStrategy& strategy_){
  if (strategy_==HiggsMassConstraint::FullCov_All_pTLambdaPhi) strategy_ = HiggsMassConstraint::FullCov_All_pTLambda;
  else if (strategy_==HiggsMassConstraint::FullCov_All_pTLambda) strategy_ = HiggsMassConstraint::FullCov_All_pTPhi;
  else if (strategy_==HiggsMassConstraint::FullCov_All_pTPhi) strategy_ = HiggsMassConstraint::FullCov_All_LambdaPhi;
  else if (strategy_==HiggsMassConstraint::FullCov_All_LambdaPhi) strategy_ = HiggsMassConstraint::FullCov_NoFSR_pTLambdaPhi;
  else if (strategy_==HiggsMassConstraint::FullCov_NoFSR_pTLambdaPhi) strategy_ = HiggsMassConstraint::FullCov_NoFSR_pTLambda;
  else if (strategy_==HiggsMassConstraint::FullCov_NoFSR_pTLambda) strategy_ = HiggsMassConstraint::FullCov_NoFSR_pTPhi;
  else if (strategy_==HiggsMassConstraint::FullCov_NoFSR_pTPhi) strategy_ = HiggsMassConstraint::FullCov_NoFSR_LambdaPhi;
  else if (strategy_==HiggsMassConstraint::FullCov_NoFSR_LambdaPhi) strategy_ = HiggsMassConstraint::FullCov_OnlyFSR_pTLambdaPhi;
  else if (strategy_==HiggsMassConstraint::FullCov_OnlyFSR_pTLambdaPhi) strategy_ = HiggsMassConstraint::FullCov_OnlyFSR_pTLambda;
  else if (strategy_==HiggsMassConstraint::FullCov_OnlyFSR_pTLambda) strategy_ = HiggsMassConstraint::FullCov_OnlyFSR_pTPhi;
  else if (strategy_==HiggsMassConstraint::FullCov_OnlyFSR_pTPhi) strategy_ = HiggsMassConstraint::FullCov_OnlyFSR_LambdaPhi;
  else if (strategy_==HiggsMassConstraint::FullCov_OnlyFSR_LambdaPhi) strategy_ = HiggsMassConstraint::CovDiagonals_All_pTLambdaPhi;
  else if (strategy_==HiggsMassConstraint::CovDiagonals_All_pTLambdaPhi) strategy_ = HiggsMassConstraint::CovDiagonals_All_pTLambda;
  else if (strategy_==HiggsMassConstraint::CovDiagonals_All_pTLambda) strategy_ = HiggsMassConstraint::CovDiagonals_All_pTPhi;
  else if (strategy_==HiggsMassConstraint::CovDiagonals_All_pTPhi) strategy_ = HiggsMassConstraint::CovDiagonals_All_LambdaPhi;
  else if (strategy_==HiggsMassConstraint::CovDiagonals_All_LambdaPhi) strategy_ = HiggsMassConstraint::CovDiagonals_All_pT;
  else if (strategy_==HiggsMassConstraint::CovDiagonals_All_pT) strategy_ = HiggsMassConstraint::CovDiagonals_All_Lambda;
  else if (strategy_==HiggsMassConstraint::CovDiagonals_All_Lambda) strategy_ = HiggsMassConstraint::CovDiagonals_All_Phi;
  else if (strategy_==HiggsMassConstraint::CovDiagonals_All_Phi) strategy_ = HiggsMassConstraint::CovDiagonals_NoFSR_pTLambdaPhi;
  else if (strategy_==HiggsMassConstraint::CovDiagonals_NoFSR_pTLambdaPhi) strategy_ = HiggsMassConstraint::CovDiagonals_NoFSR_pTLambda;
  else if (strategy_==HiggsMassConstraint::CovDiagonals_NoFSR_pTLambda) strategy_ = HiggsMassConstraint::CovDiagonals_NoFSR_pTPhi;
  else if (strategy_==HiggsMassConstraint::CovDiagonals_NoFSR_pTPhi) strategy_ = HiggsMassConstraint::CovDiagonals_NoFSR_LambdaPhi;
  else if (strategy_==HiggsMassConstraint::CovDiagonals_NoFSR_LambdaPhi) strategy_ = HiggsMassConstraint::CovDiagonals_NoFSR_pT;
  else if (strategy_==HiggsMassConstraint::CovDiagonals_NoFSR_pT) strategy_ = HiggsMassConstraint::CovDiagonals_NoFSR_Lambda;
  else if (strategy_==HiggsMassConstraint::CovDiagonals_NoFSR_Lambda) strategy_ = HiggsMassConstraint::CovDiagonals_NoFSR_Phi;
  else if (strategy_==HiggsMassConstraint::CovDiagonals_NoFSR_Phi) strategy_ = HiggsMassConstraint::CovDiagonals_OnlyFSR_pTLambdaPhi;
  else if (strategy_==HiggsMassConstraint::CovDiagonals_OnlyFSR_pTLambdaPhi) strategy_ = HiggsMassConstraint::CovDiagonals_OnlyFSR_pTLambda;
  else if (strategy_==HiggsMassConstraint::CovDiagonals_OnlyFSR_pTLambda) strategy_ = HiggsMassConstraint::CovDiagonals_OnlyFSR_pTPhi;
  else if (strategy_==HiggsMassConstraint::CovDiagonals_OnlyFSR_pTPhi) strategy_ = HiggsMassConstraint::CovDiagonals_OnlyFSR_LambdaPhi;
  else if (strategy_==HiggsMassConstraint::CovDiagonals_OnlyFSR_LambdaPhi) strategy_ = HiggsMassConstraint::CovDiagonals_OnlyFSR_pT;
  else if (strategy_==HiggsMassConstraint::CovDiagonals_OnlyFSR_pT) strategy_ = HiggsMassConstraint::CovDiagonals_OnlyFSR_Lambda;
  else if (strategy_==HiggsMassConstraint::CovDiagonals_OnlyFSR_Lambda) strategy_ = HiggsMassConstraint::CovDiagonals_OnlyFSR_Phi;
  else if (strategy_==HiggsMassConstraint::CovDiagonals_OnlyFSR_Phi) strategy_ = HiggsMassConstraint::nFitMomentumStrategies;
  else strategy_ = HiggsMassConstraint::nFitMomentumStrategies;
}
void HiggsMassConstraint::incrementMomentumStrategy(HiggsMassConstraint::FitMomentumStrategy& strategy_){
  if (strategy_==HiggsMassConstraint::FullCov_All_pTLambdaPhi) strategy_ = HiggsMassConstraint::nFitMomentumStrategies;
  else if (strategy_==HiggsMassConstraint::FullCov_All_pTLambda) strategy_ = HiggsMassConstraint::FullCov_All_pTLambdaPhi;
  else if (strategy_==HiggsMassConstraint::FullCov_All_pTPhi) strategy_ = HiggsMassConstraint::FullCov_All_pTLambda;
  else if (strategy_==HiggsMassConstraint::FullCov_All_LambdaPhi) strategy_ = HiggsMassConstraint::FullCov_All_pTPhi;
  else if (strategy_==HiggsMassConstraint::FullCov_NoFSR_pTLambdaPhi) strategy_ = HiggsMassConstraint::FullCov_All_LambdaPhi;
  else if (strategy_==HiggsMassConstraint::FullCov_NoFSR_pTLambda) strategy_ = HiggsMassConstraint::FullCov_NoFSR_pTLambdaPhi;
  else if (strategy_==HiggsMassConstraint::FullCov_NoFSR_pTPhi) strategy_ = HiggsMassConstraint::FullCov_NoFSR_pTLambda;
  else if (strategy_==HiggsMassConstraint::FullCov_NoFSR_LambdaPhi) strategy_ = HiggsMassConstraint::FullCov_NoFSR_pTPhi;
  else if (strategy_==HiggsMassConstraint::FullCov_OnlyFSR_pTLambdaPhi) strategy_ = HiggsMassConstraint::FullCov_NoFSR_LambdaPhi;
  else if (strategy_==HiggsMassConstraint::FullCov_OnlyFSR_pTLambda) strategy_ = HiggsMassConstraint::FullCov_OnlyFSR_pTLambdaPhi;
  else if (strategy_==HiggsMassConstraint::FullCov_OnlyFSR_pTPhi) strategy_ = HiggsMassConstraint::FullCov_OnlyFSR_pTLambda;
  else if (strategy_==HiggsMassConstraint::FullCov_OnlyFSR_LambdaPhi) strategy_ = HiggsMassConstraint::FullCov_OnlyFSR_pTPhi;
  else if (strategy_==HiggsMassConstraint::CovDiagonals_All_pTLambdaPhi) strategy_ = HiggsMassConstraint::FullCov_OnlyFSR_LambdaPhi;
  else if (strategy_==HiggsMassConstraint::CovDiagonals_All_pTLambda) strategy_ = HiggsMassConstraint::CovDiagonals_All_pTLambdaPhi;
  else if (strategy_==HiggsMassConstraint::CovDiagonals_All_pTPhi) strategy_ = HiggsMassConstraint::CovDiagonals_All_pTLambda;
  else if (strategy_==HiggsMassConstraint::CovDiagonals_All_LambdaPhi) strategy_ = HiggsMassConstraint::CovDiagonals_All_pTPhi;
  else if (strategy_==HiggsMassConstraint::CovDiagonals_All_pT) strategy_ = HiggsMassConstraint::CovDiagonals_All_LambdaPhi;
  else if (strategy_==HiggsMassConstraint::CovDiagonals_All_Lambda) strategy_ = HiggsMassConstraint::CovDiagonals_All_pT;
  else if (strategy_==HiggsMassConstraint::CovDiagonals_All_Phi) strategy_ = HiggsMassConstraint::CovDiagonals_All_Lambda;
  else if (strategy_==HiggsMassConstraint::CovDiagonals_NoFSR_pTLambdaPhi) strategy_ = HiggsMassConstraint::CovDiagonals_All_Phi;
  else if (strategy_==HiggsMassConstraint::CovDiagonals_NoFSR_pTLambda) strategy_ = HiggsMassConstraint::CovDiagonals_NoFSR_pTLambdaPhi;
  else if (strategy_==HiggsMassConstraint::CovDiagonals_NoFSR_pTPhi) strategy_ = HiggsMassConstraint::CovDiagonals_NoFSR_pTLambda;
  else if (strategy_==HiggsMassConstraint::CovDiagonals_NoFSR_LambdaPhi) strategy_ = HiggsMassConstraint::CovDiagonals_NoFSR_pTPhi;
  else if (strategy_==HiggsMassConstraint::CovDiagonals_NoFSR_pT) strategy_ = HiggsMassConstraint::CovDiagonals_NoFSR_LambdaPhi;
  else if (strategy_==HiggsMassConstraint::CovDiagonals_NoFSR_Lambda) strategy_ = HiggsMassConstraint::CovDiagonals_NoFSR_pT;
  else if (strategy_==HiggsMassConstraint::CovDiagonals_NoFSR_Phi) strategy_ = HiggsMassConstraint::CovDiagonals_NoFSR_Lambda;
  else if (strategy_==HiggsMassConstraint::CovDiagonals_OnlyFSR_pTLambdaPhi) strategy_ = HiggsMassConstraint::CovDiagonals_NoFSR_Phi;
  else if (strategy_==HiggsMassConstraint::CovDiagonals_OnlyFSR_pTLambda) strategy_ = HiggsMassConstraint::CovDiagonals_OnlyFSR_pTLambdaPhi;
  else if (strategy_==HiggsMassConstraint::CovDiagonals_OnlyFSR_pTPhi) strategy_ = HiggsMassConstraint::CovDiagonals_OnlyFSR_pTLambda;
  else if (strategy_==HiggsMassConstraint::CovDiagonals_OnlyFSR_LambdaPhi) strategy_ = HiggsMassConstraint::CovDiagonals_OnlyFSR_pTPhi;
  else if (strategy_==HiggsMassConstraint::CovDiagonals_OnlyFSR_pT) strategy_ = HiggsMassConstraint::CovDiagonals_OnlyFSR_LambdaPhi;
  else if (strategy_==HiggsMassConstraint::CovDiagonals_OnlyFSR_Lambda) strategy_ = HiggsMassConstraint::CovDiagonals_OnlyFSR_pT;
  else if (strategy_==HiggsMassConstraint::CovDiagonals_OnlyFSR_Phi) strategy_ = HiggsMassConstraint::CovDiagonals_OnlyFSR_Lambda;
  else strategy_ = HiggsMassConstraint::nFitMomentumStrategies;
}

void HiggsMassConstraint::setFitVVStrategy(HiggsMassConstraint::FitVVStrategy fitVVStrategy_){ fitVVStrategy=fitVVStrategy_; }
void HiggsMassConstraint::testFitVVStrategy(Int_t& fitV1, Int_t& fitV2) const{
  if (
    fitVVStrategy==HiggsMassConstraint::Fit_All_V1V2 ||
    fitVVStrategy==HiggsMassConstraint::Fit_All_V1
    ) fitV1=1;
  else fitV1=0;
  if (
    fitVVStrategy==HiggsMassConstraint::Fit_All_V1V2 ||
    fitVVStrategy==HiggsMassConstraint::Fit_All_V2
    ) fitV2=1;
  else fitV2=0;
}

void HiggsMassConstraint::sortDaughters(const std::vector<pair<const reco::Candidate*, const pat::PFParticle*>>& FermionWithFSR, std::vector<int>& order) const{
  int nDaughtersBooked=0;
  int nDaughters=(int)FermionWithFSR.size();
  int tmpDindex[2]={ 0 };
  int ordering[4]={ -1, -1, -1, -1 };
  const pair<const reco::Candidate*, const pat::PFParticle*>* df[2]={ 0 };
  const pair<const reco::Candidate*, const pat::PFParticle*>* ds[2]={ 0 };
  if (FermionWithFSR.size()>0){
    df[0] = &(FermionWithFSR.at(0));
    ordering[nDaughtersBooked] = 0;
    nDaughtersBooked++;
    for (int j=1; j<nDaughters; j++){
      const pair<const reco::Candidate*, const pat::PFParticle*>* dtmp = &(FermionWithFSR.at(j));
      if (
        std::abs(dtmp->first->pdgId())==std::abs(df[0]->first->pdgId()) // First daughter in ZZ/ZG/GG/ff requires identical |Q| and |id|.
        ){
        df[1] = dtmp;
        tmpDindex[1] = j;
        ordering[nDaughtersBooked] = j;
        nDaughtersBooked++;
        break;
      }
    }
  }
  int sindex=0;
  for (int j=1; j<nDaughters; j++){
    if (j==tmpDindex[1]) continue;
    const pair<const reco::Candidate*, const pat::PFParticle*>* dtmp = &(FermionWithFSR.at(j));
    ds[sindex] = dtmp;
    ordering[nDaughtersBooked] = j;
    nDaughtersBooked++;
    sindex++;
    if (sindex==2) break;
  }
  
  if (nDaughtersBooked!=nDaughters){
    if (nDaughters>4) std::cout << "HiggsMassConstraint::sortDaughters: Number of daughters passed " << nDaughters << ">4 is currently not supported." << std::endl;
    std::cout << "HiggsMassConstraint::sortDaughters: Number of daughters passed (" << nDaughters << ") is not the same as number of daughters booked (" << nDaughtersBooked << ")! Aborting, no daughters can be recorded." << std::endl;
    assert(0);
  }

  TLorentzVector vdf[2]={ TLorentzVector(0, 0, 0, 0), TLorentzVector(0, 0, 0, 0) };
  TLorentzVector vds[2]={ TLorentzVector(0, 0, 0, 0), TLorentzVector(0, 0, 0, 0) };
  for (int ip=0; ip<2; ip++){
    if (df[ip]!=0){
      if (df[ip]->first!=0){
        TLorentzVector vtmp;
        vtmp.SetXYZT(df[ip]->first->px(), df[ip]->first->py(), df[ip]->first->pz(), df[ip]->first->energy());
        vdf[ip] = vdf[ip] + vtmp;
      }
      if (df[ip]->second!=0){
        TLorentzVector vtmp;
        vtmp.SetXYZT(df[ip]->second->px(), df[ip]->second->py(), df[ip]->second->pz(), df[ip]->second->energy());
        vdf[ip] = vdf[ip] + vtmp;
      }
    }
    if (ds[ip]!=0){
      if (ds[ip]->first!=0){
        TLorentzVector vtmp;
        vtmp.SetXYZT(ds[ip]->first->px(), ds[ip]->first->py(), ds[ip]->first->pz(), ds[ip]->first->energy());
        vds[ip] = vds[ip] + vtmp;
      }
      if (ds[ip]->second!=0){
        TLorentzVector vtmp;
        vtmp.SetXYZT(ds[ip]->second->px(), ds[ip]->second->py(), ds[ip]->second->pz(), ds[ip]->second->energy());
        vds[ip] = vds[ip] + vtmp;
      }
    }
  }
  if (
    (df[0]!=0 && df[1]!=0)
    &&
    (
    // Order by ubar(0)v(1)
    df[0]->first->pdgId()<df[1]->first->pdgId()
    ||
    ((df[0]->first->pdgId()*df[1]->first->pdgId()>0 || (df[0]->first->pdgId()==0 && df[1]->first->pdgId()==0)) && vdf[0].Phi()<vdf[1].Phi())
    )
    ){
    swap(df[0], df[1]);
    swap(vdf[0], vdf[1]);
    swap(ordering[0], ordering[1]);
  }
  if (
    (ds[0]!=0 && ds[1]!=0)
    &&
    (
    // Order by ubar(0)v(1)
    ds[0]->first->pdgId()<ds[1]->first->pdgId()
    ||
    ((ds[0]->first->pdgId()*ds[1]->first->pdgId()>0 || (ds[0]->first->pdgId()==0 && ds[1]->first->pdgId()==0)) && vds[0].Phi()<vds[1].Phi())
    )
    ){
    swap(ds[0], ds[1]);
    swap(vds[0], vds[1]);
    swap(ordering[2], ordering[3]);
  }

  for (int ip=0; ip<4; ip++){
    if (ordering[ip]>=0){
      order.push_back(ordering[ip]);
#if hmc_debug==1
      cout << "HiggsMassConstraint::sortDaughters: Ordering of particle " << ip << ": " << ordering[ip] << endl;
#endif
    }
  }
}
void HiggsMassConstraint::addDaughters(std::vector<pair<const reco::Candidate*, const pat::PFParticle*>>& FermionWithFSR, bool fitRetry){ // Candidate supports jets as well! FSR is also a reco::Candidate daughter.
  // If the current trial is fresh, reset relevant variables.
  if (!fitRetry){
    setWorkingFitMomentumStrategy(fitMomStrategy);
    inputRaw_Fermion_FSR.clear();
  }

  // Reset initial covariance matrix
  resetInitialCovarianceMatrix();

  // Check the fit strategy
  Int_t useFullCov, FermFSRType, fitpT, fitlambda, fitphi;
  testFitMomentumStrategy(useFullCov, FermFSRType, fitpT, fitlambda, fitphi);
  Int_t fitV1, fitV2;
  testFitVVStrategy(fitV1, fitV2);

  int ndaughters=0;
  // Initialize PDG ids and obs. mom.
  for (int ix=0; ix<2; ix++){ for (int iy=0; iy<2; iy++) pdgid_ferm[ix][iy]=pdgUnknown; }

  vector<int> order;
  sortDaughters(FermionWithFSR, order);
  for (unsigned int idau=0; idau<order.size(); idau++){
    int iord = order.at(idau);
    const reco::Candidate* fermion = FermionWithFSR.at(iord).first;
    int iZ = idau/2;
    int iferm = (idau%2);
    if (fermion==0){ cerr << "HiggsMassConstraint::addDaughters : Daughter " << ndaughters << " pointer is 0!" << endl; break; }
    else{
      ndaughters++;
      if (ndaughters>4){ cerr << "HiggsMassConstraint::addDaughters : Number of daughters (" << ndaughters << ") exceeds 4!" << endl; break; }
      else{
        if (!fitRetry) inputRaw_Fermion_FSR.push_back(FermionWithFSR.at(iord));

        // ndaughters=1..4
        // Set PDG id
        int pdgid = fermion->pdgId();
        pdgid_ferm[iZ][iferm] = pdgid;
        if (!(abs(pdgid_ferm[iZ][iferm])==pdgEle || abs(pdgid_ferm[iZ][iferm])==pdgMu || abs(pdgid_ferm[iZ][iferm])==pdgTau)) pdgid_ferm[iZ][iferm]=pdgJet; // Needs to be more thorough if jets are passsed
#if hmc_debug==1
        cout << "HiggsMassConstraint::addDaughters : Daughter with PDG id " << pdgid << " is assigned to iZ=" << iZ << " and iferm=" << iferm << " with effective id " << pdgid_ferm[iZ][iferm] << endl;
        cout << "HiggsMassConstraint::addDaughters : Daughter (px, py, pz, E) = (" << fermion->px() << '\t' << fermion->py() << '\t' << fermion->pz() << '\t' << fermion->energy() << endl;
        cout << "HiggsMassConstraint::addDaughters : Daughter (pT, theta, phi, mass) = (" << fermion->pt() << '\t' << fermion->theta() << " (lambda=" << (piovertwo_val-fermion->theta()) << ")" << '\t' << fermion->phi() << '\t' << fermion->mass() << endl;
#endif

        // Set bar-momenta ranges
        massbar_ferm[iZ][iferm]->setConstant(false);
        massbar_ferm[iZ][iferm]->setRange(-sqrts, sqrts);
        massbar_ferm[iZ][iferm]->setVal(fermion->mass());
        massbar_ferm[iZ][iferm]->setRange(fermion->mass(), fermion->mass());
        massbar_ferm[iZ][iferm]->setConstant(true);

        // Set observed momenta
        pTobs_ferm[iZ][iferm]->setConstant(false);
        lambdaobs_ferm[iZ][iferm]->setConstant(false);
        phiobs_ferm[iZ][iferm]->setConstant(false);
        if (abs(pdgid_ferm[iZ][iferm])==pdgEle){
          pTobs_ferm[iZ][iferm]->setRange(pTcut_electron, sqrts);
          lambdaobs_ferm[iZ][iferm]->setRange(-lambdacut_electron, lambdacut_electron);
        }
        else if (abs(pdgid_ferm[iZ][iferm])==pdgMu){
          pTobs_ferm[iZ][iferm]->setRange(pTcut_muon, sqrts);
          lambdaobs_ferm[iZ][iferm]->setRange(-lambdacut_muon, lambdacut_muon);
        }
        else if (abs(pdgid_ferm[iZ][iferm])==pdgJet){
          pTobs_ferm[iZ][iferm]->setRange(pTcut_jet, sqrts);
          lambdaobs_ferm[iZ][iferm]->setRange(-lambdacut_jet, lambdacut_jet);
        }
        phiobs_ferm[iZ][iferm]->setRange(-pi_val, pi_val);
        pTobs_ferm[iZ][iferm]->setVal(fermion->pt());
        lambdaobs_ferm[iZ][iferm]->setVal(piovertwo_val-fermion->theta());
        phiobs_ferm[iZ][iferm]->setVal(fermion->phi());
        pTobs_ferm[iZ][iferm]->setConstant(true);
        lambdaobs_ferm[iZ][iferm]->setConstant(true);
        phiobs_ferm[iZ][iferm]->setConstant(true);

        // Set refit fermion momenta ranges
        bool fixAll = (iZ==0 && fitV1==0) || (iZ==1 && fitV2==0);
        pT_ferm[iZ][iferm]->setConstant(false);
        lambda_ferm[iZ][iferm]->setConstant(false);
        phi_ferm[iZ][iferm]->setConstant(false);
        if (abs(pdgid_ferm[iZ][iferm])==pdgEle){
          pT_ferm[iZ][iferm]->setRange(pTcut_electron, sqrts);
          lambda_ferm[iZ][iferm]->setRange(-lambdacut_electron, lambdacut_electron);
        }
        else if (abs(pdgid_ferm[iZ][iferm])==pdgMu){
          pT_ferm[iZ][iferm]->setRange(pTcut_muon, sqrts);
          lambda_ferm[iZ][iferm]->setRange(-lambdacut_muon, lambdacut_muon);
        }
        else if (abs(pdgid_ferm[iZ][iferm])==pdgJet){
          pT_ferm[iZ][iferm]->setRange(pTcut_jet, sqrts);
          lambda_ferm[iZ][iferm]->setRange(-lambdacut_jet, lambdacut_jet);
        }
        phi_ferm[iZ][iferm]->setRange(-pi_val, pi_val);
        // Initialize refit fermion momenta to the values of fermion obs-momenta
        pT_ferm[iZ][iferm]->setVal(pTobs_ferm[iZ][iferm]->getVal());
        lambda_ferm[iZ][iferm]->setVal(lambdaobs_ferm[iZ][iferm]->getVal());
        phi_ferm[iZ][iferm]->setVal(phiobs_ferm[iZ][iferm]->getVal());
        if (fixAll){
          pT_ferm[iZ][iferm]->setConstant(true);
          lambda_ferm[iZ][iferm]->setConstant(true);
          phi_ferm[iZ][iferm]->setConstant(true);
        }
        else{
          if (fitpT==0 || FermFSRType==1) pT_ferm[iZ][iferm]->setConstant(true);
          if (fitlambda==0 || FermFSRType==1) lambda_ferm[iZ][iferm]->setConstant(true);
          if (fitphi==0 || FermFSRType==1) phi_ferm[iZ][iferm]->setConstant(true);
        }

        // Get fermion covariance matrices in terms of pT, lambda and phi
        Double_t coefMat_ferm[9] ={ 0 };
        sortGetCovarianceMatrix(coefMat_ferm, fermion);
#if hmc_debug==1
        cout << "HiggsMassConstraint::addDaughters: Daughter Z" << iZ << iferm << " input covariance matrix:" << endl;
        for (int ix=0; ix<3; ix++){ for (int iy=0; iy<3; iy++) cout << coefMat_ferm[3*ix+iy] << '\t'; cout << endl; }
#endif
        setInitialCovarianceMatrix(iZ, iferm, 0, coefMat_ferm);
        if (FermFSRType!=1 && !fixAll){
          if (coefMat_ferm[3*0+0]==0. && coefMat_ferm[3*0+1]==0. && coefMat_ferm[3*0+2]==0.) pT_ferm[iZ][iferm]->setConstant(true);
          if (coefMat_ferm[3*1+1]==0. && coefMat_ferm[3*1+1]==0. && coefMat_ferm[3*1+2]==0.) lambda_ferm[iZ][iferm]->setConstant(true);
          if (coefMat_ferm[3*2+2]==0. && coefMat_ferm[3*1+2]==0. && coefMat_ferm[3*2+2]==0.) phi_ferm[iZ][iferm]->setConstant(true);
          if (!pT_ferm[iZ][iferm]->isConstant()){
            pT_ferm[iZ][iferm]->setRange(
              max(pT_ferm[iZ][iferm]->getMin(), pT_ferm[iZ][iferm]->getVal()-5.*sqrt(coefMat_ferm[3*0+0])),
              min(pT_ferm[iZ][iferm]->getMax(), pT_ferm[iZ][iferm]->getVal()+5.*sqrt(coefMat_ferm[3*0+0]))
              );
            pTobs_ferm[iZ][iferm]->setRange(pT_ferm[iZ][iferm]->getMin(), pT_ferm[iZ][iferm]->getMax());
          }
          if (!lambda_ferm[iZ][iferm]->isConstant()){
            lambda_ferm[iZ][iferm]->setRange(
              max(lambda_ferm[iZ][iferm]->getMin(), lambda_ferm[iZ][iferm]->getVal()-5.*sqrt(coefMat_ferm[3*1+1])),
              min(lambda_ferm[iZ][iferm]->getMax(), lambda_ferm[iZ][iferm]->getVal()+5.*sqrt(coefMat_ferm[3*1+1]))
              );
            lambdaobs_ferm[iZ][iferm]->setRange(lambda_ferm[iZ][iferm]->getMin(), lambda_ferm[iZ][iferm]->getMax());
          }
          if (!phi_ferm[iZ][iferm]->isConstant()){
            phi_ferm[iZ][iferm]->setRange(
              max(phi_ferm[iZ][iferm]->getMin(), phi_ferm[iZ][iferm]->getVal()-5.*sqrt(coefMat_ferm[3*2+2])),
              min(phi_ferm[iZ][iferm]->getMax(), phi_ferm[iZ][iferm]->getVal()+5.*sqrt(coefMat_ferm[3*2+2]))
              );
            phiobs_ferm[iZ][iferm]->setRange(phi_ferm[iZ][iferm]->getMin(), phi_ferm[iZ][iferm]->getMax());
          }

          strategicInvertCovarianceMatrix(useFullCov, fitpT, fitlambda, fitphi, coefMat_ferm);
        }
        else{
          for (int ix=0; ix<3; ix++){ for (int iy=0; iy<3; iy++) coefMat_ferm[3*ix+iy]=0; }
        }
        setInverseCovarianceMatrix(iZ, iferm, 0, coefMat_ferm);

        // Do FSR here
        const pat::PFParticle* gamma = FermionWithFSR.at(iord).second;
        Double_t coefMat_fsr[9] ={ 0 };
        if (gamma!=0){
#if hmc_debug==1
          cout << "HiggsMassConstraint::addDaughters : An FSR is assigned to iZ=" << iZ << " and iferm=" << iferm << endl;
#endif

          // Set observed momenta
          pTobs_fsr[iZ][iferm]->setConstant(false);
          lambdaobs_fsr[iZ][iferm]->setConstant(false);
          phiobs_fsr[iZ][iferm]->setConstant(false);

          pTobs_fsr[iZ][iferm]->setRange(pTcut_fsr, sqrts);
          lambdaobs_fsr[iZ][iferm]->setRange(-lambdacut_fsr, lambdacut_fsr);
          phiobs_fsr[iZ][iferm]->setRange(-pi_val, pi_val);

          pTobs_fsr[iZ][iferm]->setVal(gamma->pt());
          lambdaobs_fsr[iZ][iferm]->setVal(piovertwo_val-gamma->theta());
          phiobs_fsr[iZ][iferm]->setVal(gamma->phi());
          pTobs_fsr[iZ][iferm]->setConstant(true);
          lambdaobs_fsr[iZ][iferm]->setConstant(true);
          phiobs_fsr[iZ][iferm]->setConstant(true);

          // Set fsr ranges within the cuts and initialize
          pT_fsr[iZ][iferm]->setConstant(false); pT_fsr[iZ][iferm]->setRange(pTcut_fsr, sqrts); pT_fsr[iZ][iferm]->setVal(pTobs_fsr[iZ][iferm]->getVal());
          lambda_fsr[iZ][iferm]->setConstant(false); lambda_fsr[iZ][iferm]->setRange(-lambdacut_fsr, lambdacut_fsr); lambda_fsr[iZ][iferm]->setVal(lambdaobs_fsr[iZ][iferm]->getVal());
          phi_fsr[iZ][iferm]->setConstant(false); phi_fsr[iZ][iferm]->setRange(-pi_val, pi_val); phi_fsr[iZ][iferm]->setVal(phiobs_fsr[iZ][iferm]->getVal());
          if (fixAll){
            pT_fsr[iZ][iferm]->setConstant(true);
            lambda_fsr[iZ][iferm]->setConstant(true);
            phi_fsr[iZ][iferm]->setConstant(true);
          }
          else{
            if (fitpT==0 || FermFSRType==0) pT_fsr[iZ][iferm]->setConstant(true);
            if (fitlambda==0 || FermFSRType==0) lambda_fsr[iZ][iferm]->setConstant(true);
            if (fitphi==0 || FermFSRType==0) phi_fsr[iZ][iferm]->setConstant(true);
          }

          // Get fsr covariance matrices in terms of pT, lambda and phi
          sortGetCovarianceMatrix(coefMat_fsr, fermion);
#if hmc_debug==1
          cout << "HiggsMassConstraint::addDaughters : Daughter Z" << iZ << iferm << " FSR input covariance matrix:" << endl;
          for (int ix=0; ix<3; ix++){ for (int iy=0; iy<3; iy++) cout << coefMat_fsr[3*ix+iy] << '\t'; cout << endl; }
#endif
          setInitialCovarianceMatrix(iZ, iferm, 1, coefMat_fsr);
          if (FermFSRType!=0 && !fixAll){
            if (coefMat_fsr[3*0+0]==0. && coefMat_fsr[3*0+1]==0. && coefMat_fsr[3*0+2]==0.) pT_fsr[iZ][iferm]->setConstant(true);
            if (coefMat_fsr[3*1+1]==0. && coefMat_fsr[3*1+1]==0. && coefMat_fsr[3*1+2]==0.) lambda_fsr[iZ][iferm]->setConstant(true);
            if (coefMat_fsr[3*2+2]==0. && coefMat_fsr[3*1+2]==0. && coefMat_fsr[3*2+2]==0.) phi_fsr[iZ][iferm]->setConstant(true);
            if (!pT_fsr[iZ][iferm]->isConstant()){
              pT_fsr[iZ][iferm]->setRange(
                max(pT_fsr[iZ][iferm]->getMin(), pT_fsr[iZ][iferm]->getVal()-5.*sqrt(coefMat_fsr[3*0+0])),
                min(pT_fsr[iZ][iferm]->getMax(), pT_fsr[iZ][iferm]->getVal()+5.*sqrt(coefMat_fsr[3*0+0]))
                );
              pTobs_fsr[iZ][iferm]->setRange(pT_fsr[iZ][iferm]->getMin(), pT_fsr[iZ][iferm]->getMax());
            }
            if (!lambda_fsr[iZ][iferm]->isConstant()){
              lambda_fsr[iZ][iferm]->setRange(
                max(lambda_fsr[iZ][iferm]->getMin(), lambda_fsr[iZ][iferm]->getVal()-5.*sqrt(coefMat_fsr[3*1+1])),
                min(lambda_fsr[iZ][iferm]->getMax(), lambda_fsr[iZ][iferm]->getVal()+5.*sqrt(coefMat_fsr[3*1+1]))
                );
              lambdaobs_fsr[iZ][iferm]->setRange(lambda_fsr[iZ][iferm]->getMin(), lambda_fsr[iZ][iferm]->getMax());
            }
            if (!phi_fsr[iZ][iferm]->isConstant()){
              phi_fsr[iZ][iferm]->setRange(
                max(phi_fsr[iZ][iferm]->getMin(), phi_fsr[iZ][iferm]->getVal()-5.*sqrt(coefMat_fsr[3*2+2])),
                min(phi_fsr[iZ][iferm]->getMax(), phi_fsr[iZ][iferm]->getVal()+5.*sqrt(coefMat_fsr[3*2+2]))
                );
              phiobs_fsr[iZ][iferm]->setRange(phi_fsr[iZ][iferm]->getMin(), phi_fsr[iZ][iferm]->getMax());
            }

            strategicInvertCovarianceMatrix(useFullCov, fitpT, fitlambda, fitphi, coefMat_fsr);
          }
          else{
            for (int ix=0; ix<3; ix++){ for (int iy=0; iy<3; iy++) coefMat_fsr[3*ix+iy]=0; }
          }
          setInverseCovarianceMatrix(iZ, iferm, 1, coefMat_fsr);
        }
        else{
          pTobs_fsr[iZ][iferm]->setConstant(false);
          lambdaobs_fsr[iZ][iferm]->setConstant(false);
          phiobs_fsr[iZ][iferm]->setConstant(false);
          pTobs_fsr[iZ][iferm]->setVal(0.);
          lambdaobs_fsr[iZ][iferm]->setVal(0.);
          phiobs_fsr[iZ][iferm]->setVal(0.);
          pTobs_fsr[iZ][iferm]->setConstant(true);
          lambdaobs_fsr[iZ][iferm]->setConstant(true);
          phiobs_fsr[iZ][iferm]->setConstant(true);

          // setRange below resets range from previous iteration
          pT_fsr[iZ][iferm]->setConstant(false); pT_fsr[iZ][iferm]->setRange(0., sqrts); pT_fsr[iZ][iferm]->setVal(0.); pT_fsr[iZ][iferm]->setRange(0., 0.); pT_fsr[iZ][iferm]->setConstant(true);
          lambda_fsr[iZ][iferm]->setConstant(false); lambda_fsr[iZ][iferm]->setRange(-piovertwo_val, piovertwo_val); lambda_fsr[iZ][iferm]->setVal(0.); lambda_fsr[iZ][iferm]->setRange(0., 0.); lambda_fsr[iZ][iferm]->setConstant(true);
          phi_fsr[iZ][iferm]->setConstant(false); phi_fsr[iZ][iferm]->setRange(-pi_val, pi_val); phi_fsr[iZ][iferm]->setVal(0.); phi_fsr[iZ][iferm]->setRange(0., 0.); phi_fsr[iZ][iferm]->setConstant(true);

          setInverseCovarianceMatrix(iZ, iferm, 1, coefMat_fsr);
        }

      }
    }
  }

  // Set those non-existing particles
  for (int iZ=0; iZ<2; iZ++){
    for (int iferm=0; iferm<2; iferm++){
      if (pdgid_ferm[iZ][iferm]!=pdgUnknown) continue;
#if hmc_debug==1
      cout << "HiggsMassConstraint::addDaughters: Particle at iZ=" << iZ << " and iferm=" << iferm << " is not assigned. Fixing it to 0." << endl;
#endif

      // Set observed momenta
      massbar_ferm[iZ][iferm]->setConstant(false);
      pTobs_ferm[iZ][iferm]->setConstant(false);
      lambdaobs_ferm[iZ][iferm]->setConstant(false);
      phiobs_ferm[iZ][iferm]->setConstant(false);
      massbar_ferm[iZ][iferm]->setVal(0.);
      pTobs_ferm[iZ][iferm]->setVal(0.);
      lambdaobs_ferm[iZ][iferm]->setVal(0.);
      phiobs_ferm[iZ][iferm]->setVal(0.);
      massbar_ferm[iZ][iferm]->setConstant(true);
      pTobs_ferm[iZ][iferm]->setConstant(true);
      lambdaobs_ferm[iZ][iferm]->setConstant(true);
      phiobs_ferm[iZ][iferm]->setConstant(true);

      pTobs_fsr[iZ][iferm]->setConstant(false);
      lambdaobs_fsr[iZ][iferm]->setConstant(false);
      phiobs_fsr[iZ][iferm]->setConstant(false);
      pTobs_fsr[iZ][iferm]->setVal(0.);
      lambdaobs_fsr[iZ][iferm]->setVal(0.);
      phiobs_fsr[iZ][iferm]->setVal(0.);
      pTobs_fsr[iZ][iferm]->setConstant(true);
      lambdaobs_fsr[iZ][iferm]->setConstant(true);
      phiobs_fsr[iZ][iferm]->setConstant(true);

      // setRange below resets range from previous iteration
      pT_ferm[iZ][iferm]->setConstant(false); pT_ferm[iZ][iferm]->setRange(0., sqrts); pT_ferm[iZ][iferm]->setVal(0.); pT_ferm[iZ][iferm]->setRange(0., 0.); pT_ferm[iZ][iferm]->setConstant(true);
      lambda_ferm[iZ][iferm]->setConstant(false); lambda_ferm[iZ][iferm]->setRange(-piovertwo_val, piovertwo_val); lambda_ferm[iZ][iferm]->setVal(0.); lambda_ferm[iZ][iferm]->setRange(0., 0.); lambda_ferm[iZ][iferm]->setConstant(true);
      phi_ferm[iZ][iferm]->setConstant(false); phi_ferm[iZ][iferm]->setRange(-pi_val, pi_val); phi_ferm[iZ][iferm]->setVal(0.); phi_ferm[iZ][iferm]->setRange(0., 0.); phi_ferm[iZ][iferm]->setConstant(true);
      pT_fsr[iZ][iferm]->setConstant(false); pT_fsr[iZ][iferm]->setRange(0., sqrts); pT_fsr[iZ][iferm]->setVal(0.); pT_fsr[iZ][iferm]->setRange(0., 0.); pT_fsr[iZ][iferm]->setConstant(true);
      lambda_fsr[iZ][iferm]->setConstant(false); lambda_fsr[iZ][iferm]->setRange(-piovertwo_val, piovertwo_val); lambda_fsr[iZ][iferm]->setVal(0.); lambda_fsr[iZ][iferm]->setRange(0., 0.); lambda_fsr[iZ][iferm]->setConstant(true);
      phi_fsr[iZ][iferm]->setConstant(false); phi_fsr[iZ][iferm]->setRange(-pi_val, pi_val); phi_fsr[iZ][iferm]->setVal(0.); phi_fsr[iZ][iferm]->setRange(0., 0.); phi_fsr[iZ][iferm]->setConstant(true);

      Double_t coefMat_ferm[9] ={ 0 };
      Double_t coefMat_fsr[9] ={ 0 };
      setInverseCovarianceMatrix(iZ, iferm, 0, coefMat_ferm);
      setInverseCovarianceMatrix(iZ, iferm, 1, coefMat_fsr);
    }
  }

#if hmc_debug==1
  summarizeDaughters();
#endif

  // Set spline variable at the very end
  mManip[2]->setConstant(false);
  for (int iv=0; iv<3; iv++) mManip[iv]->setVal(m[iv]->getVal());
  mManip[2]->setConstant(true);

  if (!fitRetry) setInitialMassErrors();
}
void HiggsMassConstraint::summarizeDaughters()const{
  cout << "=== SUMMARY OF FINAL STATES ===" << endl;

  cout << "|| Fermions ||" << endl;
  for (int iZ=0; iZ<2; iZ++){
    for (int iferm=0; iferm<2; iferm++){
      cout << "Z" << iZ << " daughter " << iferm;
      cout << " ";
      cout << "(px,py,pz,E) = ( ";
      cout << px_ferm[iZ][iferm]->getVal() << " " << py_ferm[iZ][iferm]->getVal() << " " << pz_ferm[iZ][iferm]->getVal() << " " << E_ferm[iZ][iferm]->getVal() << " )";
      cout << " = ";
      cout << "(pT,lambda,phi,m) = ( ";
      cout << pT_ferm[iZ][iferm]->getVal() << " " << lambda_ferm[iZ][iferm]->getVal() << " " << phi_ferm[iZ][iferm]->getVal() << " " << massbar_ferm[iZ][iferm]->getVal() << " )";
      cout << endl;
    }
  }
  cout << "|| FSR ||" << endl;
  for (int iZ=0; iZ<2; iZ++){
    for (int iferm=0; iferm<2; iferm++){
      cout << "Z" << iZ << " daughter " << iferm << " FSR";
      cout << " ";
      cout << "(px,py,pz,E) = ( ";
      cout << px_fsr[iZ][iferm]->getVal() << " " << py_fsr[iZ][iferm]->getVal() << " " << pz_fsr[iZ][iferm]->getVal() << " " << E_fsr[iZ][iferm]->getVal() << " )";
      cout << " = ";
      cout << "(pT,lambda,phi,m) = ( ";
      cout << pT_fsr[iZ][iferm]->getVal() << " " << lambda_fsr[iZ][iferm]->getVal() << " " << phi_fsr[iZ][iferm]->getVal() << " " << 0 << " )";
      cout << endl;
    }
  }
  cout << "|| Fermion + FSR ||" << endl;
  for (int iZ=0; iZ<2; iZ++){
    for (int iferm=0; iferm<2; iferm++){
      cout << "Z" << iZ << " daughter " << iferm << " sum";
      cout << " ";
      cout << "(px,py,pz,E,m) = ( ";
      cout << px_Hdaughter[iZ][iferm]->getVal() << " " << py_Hdaughter[iZ][iferm]->getVal() << " " << pz_Hdaughter[iZ][iferm]->getVal() << " " << E_Hdaughter[iZ][iferm]->getVal() << " " << m_Hdaughter[iZ][iferm]->getVal() << " )";
      cout << endl;
    }
  }
  cout << "|| V1/V2 functions ||" << endl;
  for (int iZ=0; iZ<2; iZ++) cout << "beta(V" << iZ+1 << ") = " << beta_Vdaughter[iZ]->getVal() << endl;
  for (int iZ=0; iZ<2; iZ++) cout << "mAB-OS(V" << iZ+1 << ") = " << mAB[iZ]->getVal() << endl;
  for (int iZ=0; iZ<2; iZ++) cout << "mAB-SS(V" << iZ+1 << ") = " << mAB[iZ+2]->getVal() << endl;
  for (int iZ=0; iZ<2; iZ++) cout << "m" << iZ+1 << " = " << m[iZ]->getVal() << endl;
  cout << "m12 = " << m[2]->getVal() << endl;

  cout << "===============================" << endl;
}

void HiggsMassConstraint::fitTo(std::vector<pair<const reco::Candidate*, const pat::PFParticle*>>& FermionWithFSR){
#if hmc_debug==1
  cout << "\n\nHiggsMassConstraint::fitTo is called." << endl;
#endif
  addDaughters(FermionWithFSR);
  fit();
#if hmc_debug==1
  cout << "HiggsMassConstraint::fitTo is terminating.\n\n" << endl;
#endif
}
void HiggsMassConstraint::getDataVariables(RooArgSet* fitVars, RooArgSet* intVars, RooArgSet* condVars) const{
  if (intVars!=0){
    if (intCodeStart%RooSpin::prime_h1 != 0) intVars->add(*(h1));
    if (intCodeStart%RooSpin::prime_h2 != 0) intVars->add(*(h2));
    if (intCodeStart%RooSpin::prime_hs != 0) intVars->add(*(hs));
    if (intCodeStart%RooSpin::prime_Phi != 0) intVars->add(*(Phi));
    if (intCodeStart%RooSpin::prime_Phi1 != 0) intVars->add(*(Phi1));
    //intVars->add(*(Y));
  }

  for (int iZ=0; iZ<2; iZ++){
    for (int iferm=0; iferm<2; iferm++){
      bool refitvarconst;

      refitvarconst = pT_ferm[iZ][iferm]->isConstant();
      if (!refitvarconst){
        if (fitVars!=0) fitVars->add(*(pT_ferm[iZ][iferm]));
        if (intVars!=0){
          pTobs_ferm[iZ][iferm]->setConstant(false);
          intVars->add(*(pTobs_ferm[iZ][iferm]));
        }
      }
      else if (refitvarconst && condVars!=0) condVars->add(*(pTobs_ferm[iZ][iferm]));

      refitvarconst = lambda_ferm[iZ][iferm]->isConstant();
      if (!refitvarconst){
        if (fitVars!=0) fitVars->add(*(lambda_ferm[iZ][iferm]));
        if (intVars!=0){
          lambdaobs_ferm[iZ][iferm]->setConstant(false);
          intVars->add(*(lambdaobs_ferm[iZ][iferm]));
        }
      }
      else if (refitvarconst && condVars!=0) condVars->add(*(lambdaobs_ferm[iZ][iferm]));

      refitvarconst = phi_ferm[iZ][iferm]->isConstant();
      if (!refitvarconst){
        if (fitVars!=0) fitVars->add(*(phi_ferm[iZ][iferm]));
        if (intVars!=0){
          phiobs_ferm[iZ][iferm]->setConstant(false);
          intVars->add(*(phiobs_ferm[iZ][iferm]));
        }
      }
      else if (refitvarconst && condVars!=0) condVars->add(*(phiobs_ferm[iZ][iferm]));

      refitvarconst = pT_fsr[iZ][iferm]->isConstant();
      if (!refitvarconst){
        if (fitVars!=0) fitVars->add(*(pT_fsr[iZ][iferm]));
        if (intVars!=0){
          pTobs_fsr[iZ][iferm]->setConstant(false);
          intVars->add(*(pTobs_fsr[iZ][iferm]));
        }
      }
      else if (refitvarconst && condVars!=0) condVars->add(*(pTobs_fsr[iZ][iferm]));

      refitvarconst = lambda_fsr[iZ][iferm]->isConstant();
      if (!refitvarconst){
        if (fitVars!=0) fitVars->add(*(lambda_fsr[iZ][iferm]));
        if (intVars!=0){
          lambdaobs_fsr[iZ][iferm]->setConstant(false);
          intVars->add(*(lambdaobs_fsr[iZ][iferm]));
        }
      }
      else if (refitvarconst && condVars!=0) condVars->add(*(lambdaobs_fsr[iZ][iferm]));

      refitvarconst = phi_fsr[iZ][iferm]->isConstant();
      if (!refitvarconst){
        if (fitVars!=0) fitVars->add(*(phi_fsr[iZ][iferm]));
        if (intVars!=0){
          phiobs_fsr[iZ][iferm]->setConstant(false);
          intVars->add(*(phiobs_fsr[iZ][iferm]));
        }
      }
      else if (refitvarconst && condVars!=0) condVars->add(*(phiobs_fsr[iZ][iferm]));
    }
  }
}
RooDataSet* HiggsMassConstraint::getDataset(RooArgSet* fitVars, RooArgSet* condVars) const{
  RooArgSet fit_args;
  RooArgSet data_args;
  RooArgSet cond_args;
  getDataVariables(&fit_args, &data_args, &cond_args);
  if (fitVars!=0){
    fitVars->removeAll();
    fitVars->add(fit_args);
  }
  if (condVars!=0){
    condVars->removeAll();
    condVars->add(cond_args);
  }
  RooDataSet* data = 0;
  if (data_args.getSize()<=15 && data_args.getSize()>0){
    data = new RooDataSet("fittedHdaughters", "", data_args);
    data->add(data_args);
#if hmc_debug==1
    cout << "HiggsMassConstraint::getDataset: Number of arguments: " << data_args.getSize() << endl;
    cout << "HiggsMassConstraint::getDataset: Data:" << endl;
    data->Print("v");
    cout << endl;
#endif
  }
  return data;
}
void HiggsMassConstraint::fit(){
#if hmc_debug==1
  cout << "Begin HiggsMassConstraint::fit" << endl;
#endif

  // Get the data to fit
  RooArgSet fitVars, conditionals;
  RooDataSet* data = getDataset(&fitVars, &conditionals);
  conditionals.add(*(mManip[2]));
  // Get the PDF to use
  RooAbsPdf* activePDF;
  if (useFastPDF){
    destroyCompoundFastPdf();
    constructCompoundFastPdf();
    activePDF = fastPDF;
#if hmc_debug==1
    cout << "HiggsMassConstraint::fit: FastPDF option is active. The PDF being used is " << activePDF->GetName() << "." << endl;
#endif
  }
  else{
    activePDF = PDF;
#if hmc_debug==1
    cout << "HiggsMassConstraint::fit: FastPDF option is inactive. The PDF being used is " << activePDF->GetName() << "." << endl;
#endif
  }
#if hmc_debug==1
  activePDF->Print("v");

  RooArgSet* pdfPars;
  TIterator* parIter;
  RooAbsArg* thePar;
  cout << "HiggsMassConstraint::fit: PDF parameters:" << endl;
  pdfPars = activePDF->getParameters((RooArgSet*)0, true);
  parIter = pdfPars->createIterator();
  while ((thePar = (RooAbsArg*)parIter->Next())) cout << thePar->GetName() << endl;
  cout << endl;
  delete parIter;
  delete pdfPars;

  cout << "HiggsMassConstraint::fit: PDF observables:" << endl;
  pdfPars = activePDF->getObservables((RooArgSet*)0, true);
  parIter = pdfPars->createIterator();
  while ((thePar = (RooAbsArg*)parIter->Next())) cout << thePar->GetName() << endl;
  cout << endl;
  delete parIter;
  delete pdfPars;

  cout << "HiggsMassConstraint::fit: PDF dependents:" << endl;
  pdfPars = activePDF->getDependents((RooArgSet*)0);
  parIter = pdfPars->createIterator();
  while ((thePar = (RooAbsArg*)parIter->Next())) cout << thePar->GetName() << endl;
  cout << endl;
  delete parIter;
  delete pdfPars;
#endif

  // Conditional variables (i.e. m12)
  //RooArgSet conditionals;
  //conditionals.add(*(m[2]));

  // Delete the fitResult in case it exists, just to avoid unwanted memory leaks.
  deletePtr(fitResult);

  // Fix the factory parameters
  if (hvvFactory!=0){ hvvFactory->makeParamsConst(true); hvvFactory->makeCouplingsConst(true); }
  if (hvvFastFactory!=0){ hvvFastFactory->makeParamsConst(true); hvvFastFactory->makeCouplingsConst(true); }
  if (xvvFactory!=0){ xvvFactory->makeParamsConst(true); xvvFactory->makeCouplingsConst(true); }

  /******************************************** BEGIN FIT **************************************************/
  const Int_t minimizerSuccess=0;
  Int_t fitStatus=-999;
  // Set the fit commands
#if hmc_debug==1
  cout << "HiggsMassConstraint::fit: Attempting first fit." << endl;
#endif
  RooLinkedList cmdList;
  RooCmdArg condObsArg = RooFit::ConditionalObservables(conditionals); cmdList.Add((TObject*)&condObsArg);
  RooCmdArg constrArg = RooFit::Constrain(fitVars); cmdList.Add((TObject*)&constrArg); // All fit variables should be constrained!
  RooCmdArg saveArg = RooFit::Save(true); cmdList.Add((TObject*)&saveArg);
  RooCmdArg hesseArg = RooFit::Hesse(true); cmdList.Add((TObject*)&hesseArg);
  //RooCmdArg minimizerArg = RooFit::Minimizer("Minuit", "migrad"); cmdList.Add((TObject*)&minimizerArg);
  //RooCmdArg minimizerStrategyArg = RooFit::Strategy(0); cmdList.Add((TObject*)&minimizerStrategyArg);
  // Misc. options
#if hmc_debug==1
  RooCmdArg timerArg = RooFit::Timer(true); cmdList.Add((TObject*)&timerArg);
  RooCmdArg printlevelArg = RooFit::PrintLevel(3); cmdList.Add((TObject*)&printlevelArg);
#else
  RooCmdArg printlevelArg = RooFit::PrintLevel(-1); cmdList.Add((TObject*)&printlevelArg);
#endif
  // Try with the default strategy
  if (data!=0){
    fitResult = activePDF->fitTo(*data, cmdList);
    fitStatus = fitResult->status();
  }
#if hmc_debug==1
  cout << "HiggsMassConstraint::fit: First fit attempted." << endl;
  if (fitResult!=0){
    cout << "Fit status is " << fitStatus << endl;
    cout << "Fit properties:" << endl;
    fitResult->Print("v");
  }
#endif
  // If the default strategy fails, decrement it until there is no strategy.
  while (fitStatus!=minimizerSuccess){
    decrementMomentumStrategy(fitMomStrategy_final);
    if (fitMomStrategy_final==HiggsMassConstraint::nFitMomentumStrategies) break;
    cerr << "HiggsMassConstraint::fit: Fit did not converge or was not valid. Status changed from " << fitMomStrategy << " to " << fitMomStrategy_final << " to retry." << endl;

    deletePtr(data); deletePtr(fitResult); cmdList.Clear();
    vector<pair<const reco::Candidate*, const pat::PFParticle*>> pseudoinput = inputRaw_Fermion_FSR;
    addDaughters(pseudoinput, true);
    data = getDataset(&fitVars, 0);
    condObsArg = RooFit::ConditionalObservables(conditionals); cmdList.Add((TObject*)&condObsArg);
    constrArg = RooFit::Constrain(fitVars); cmdList.Add((TObject*)&constrArg); // All fit variables should be constrained!
    cmdList.Add((TObject*)&saveArg);
    cmdList.Add((TObject*)&hesseArg);
    //cmdList.Add((TObject*)&minimizerArg);
    //cmdList.Add((TObject*)&minimizerStrategyArg);
    cmdList.Add((TObject*)&printlevelArg);
    if (data!=0){
      fitResult = activePDF->fitTo(*data, cmdList);
      fitStatus = fitResult->status();
    }
    else fitStatus=-999;
  }
  // If decrementing the strategy fails, increment it instead until there is no strategy.
  if (fitStatus!=minimizerSuccess) setWorkingFitMomentumStrategy(fitMomStrategy);
  while (fitStatus!=minimizerSuccess){
    incrementMomentumStrategy(fitMomStrategy_final);
    if (fitMomStrategy_final==HiggsMassConstraint::nFitMomentumStrategies) break;
    cerr << "HiggsMassConstraint::fit: Fit did not converge. Status changed from " << fitMomStrategy << " to " << fitMomStrategy_final << " to retry." << endl;

    deletePtr(data); deletePtr(fitResult); cmdList.Clear();
    vector<pair<const reco::Candidate*, const pat::PFParticle*>> pseudoinput = inputRaw_Fermion_FSR;
    addDaughters(pseudoinput, true);
    data = getDataset(&fitVars, 0);
    condObsArg = RooFit::ConditionalObservables(conditionals); cmdList.Add((TObject*)&condObsArg);
    constrArg = RooFit::Constrain(fitVars); cmdList.Add((TObject*)&constrArg); // All fit variables should be constrained!
    cmdList.Add((TObject*)&saveArg);
    cmdList.Add((TObject*)&hesseArg);
    //cmdList.Add((TObject*)&minimizerArg);
    //cmdList.Add((TObject*)&minimizerStrategyArg);
    cmdList.Add((TObject*)&printlevelArg);
    if (data!=0){
      fitResult = activePDF->fitTo(*data, cmdList);
      fitStatus = fitResult->status();
    }
    else fitStatus=-999;
  }
  /********************************************  END FIT  **************************************************/

  Int_t nDimCovMat=0;
  Int_t nFinalPars=0;
  bool fitFinalSuccess = (fitResult!=0 && fitStatus==minimizerSuccess);

  if (fitFinalSuccess){
    TMatrixDSym mat_tmp = fitResult->covarianceMatrix();;
    nDimCovMat = mat_tmp.GetNcols();
#if hmc_debug==1
    cout << "Number of columns in the unprocessed covariance matrix is " << nDimCovMat << ". The covariance matrix is" << endl;
    for (int ix=0; ix<nDimCovMat; ix++){
      for (int iy=0; iy<nDimCovMat; iy++) cout << mat_tmp[ix][iy] << '\t';
      cout << endl;
    }
    cout << endl;
#endif
    fitCovMatrix.ResizeTo(nDimCovMat, nDimCovMat);
    fitCovMatrix = mat_tmp;
    fitFinalSuccess = fitFinalSuccess && (nDimCovMat>0);
  }
  if (!fitFinalSuccess) cout << "Fit did not converge after all trials. Default parameters are to be used." << endl;

  if (fitFinalSuccess){
    const RooArgList pars = fitResult->floatParsFinal();
    nFinalPars = pars.getSize();
    for (int ip=0; ip<nFinalPars; ip++){
      const RooAbsArg* arg = pars.at(ip);
      if (dynamic_cast<const RooRealVar*>(arg)==0){ cerr << "Parameter " << ip << " (" << arg->GetName() << ") is not a RooRealVar!" << endl; assert(0); }
#if hmc_debug==1
      else{
        cout << "Parameter " << arg->GetName() << " = " << dynamic_cast<const RooRealVar*>(arg)->getVal() << " +- " << dynamic_cast<const RooRealVar*>(arg)->getError() << endl;
      }
#endif
    }

    if (nFinalPars!=nDimCovMat){ cerr << "Number of columns in the covariance matrix (" << nDimCovMat << ") does not match with the number of final floating variables (" << nFinalPars << ")!" << endl; assert(0); }
    else{
#if hmc_debug==1
      cout << "Number of columns in the covariance matrix is " << nDimCovMat << ". The covariance matrix is" << endl;
      for (int ix=0; ix<nDimCovMat; ix++){
        for (int iy=0; iy<nDimCovMat; iy++) cout << fitCovMatrix[ix][iy] << '\t';
        cout << endl;
      }
      cout << endl;
#endif
      // Re-order the covariance matrix here
      standardOrderedFinalCovarianceMatrix(pars);
    }
  }

  // Relax the factory parameters and return to normal
  if (hvvFactory!=0){ hvvFactory->makeParamsConst(false); hvvFactory->makeCouplingsConst(false); }
  if (hvvFastFactory!=0){ hvvFastFactory->makeParamsConst(false); hvvFastFactory->makeCouplingsConst(false); }
  if (xvvFactory!=0){ xvvFactory->makeParamsConst(false); xvvFactory->makeCouplingsConst(false); }
  deletePtr(data);
}

void HiggsMassConstraint::sortGetCovarianceMatrix(double (&momCov)[9], const reco::Candidate* particle){
  const reco::GsfElectron* electron = dynamic_cast<const reco::GsfElectron*>(particle);
  const reco::PFCandidate* pfcand = dynamic_cast<const reco::PFCandidate*>(particle);
  const pat::Muon* muon = dynamic_cast<const pat::Muon*>(particle);
  //const pat::Electron* pat_electron = dynamic_cast<const pat::Electron*>(particle);
  const pat::Jet* jet = dynamic_cast<const pat::Jet*>(particle);

  if (jet!=0) return getCovarianceMatrix(momCov, jet);
  else if (muon!=0) return getCovarianceMatrix(momCov, muon);
  else if (pfcand!=0) return getCovarianceMatrix(momCov, pfcand); // This is a general PFCandidate, which is mostly for use as a photon
  //else if (pat_electron!=0) return getCovarianceMatrix(momCov, pat_electron);
  else if (electron!=0) return getCovarianceMatrix(momCov, electron);
  else{
#if hmc_debug==1
    cout << "HiggsMassConstraint::sortGetCovarianceMatrix: Could not determine the type of particle " << particle << ". Setting all covariance matrices to 0." << endl;
#endif
    for(int i=0;i<9;i++) momCov[i]=0.;
  }
}
void HiggsMassConstraint::getCovarianceMatrix(double (&momCov)[9], const reco::GsfElectron* particle){
  for(int i=0;i<9;i++) momCov[i]=0.;

#if hmc_debug==1
  cout << "Begin HiggsMassConstraint::getCovarianceMatrix(const reco::GsfElectron* particle) with argument " << particle << endl;
#endif

  double lambda = piovertwo_val - particle->theta();
  double energyerr;
  if (particle->ecalDriven()) energyerr = particle->p4Error(reco::GsfElectron::P4_COMBINATION);
  else{
    double ecalEnergy = particle->correctedEcalEnergy();
    double err2 = 0.;
    if (particle->isEB()){
      err2 += (5.24e-02*5.24e-02)/ecalEnergy;
      err2 += (2.01e-01*2.01e-01)/(ecalEnergy*ecalEnergy);
      err2 += 1.00e-02*1.00e-02;
    }
    else if (particle->isEE()){
      err2 += (1.46e-01*1.46e-01)/ecalEnergy;
      err2 += (9.21e-01*9.21e-01)/(ecalEnergy*ecalEnergy);
      err2 += 1.94e-03*1.94e-03;
    }
    energyerr = ecalEnergy * sqrt(err2);
  }
#if hmc_debug==1
  cout << "HiggsMassConstraint::getCovarianceMatrix(const reco::GsfElectron* particle): Energy error before GsfTrack: " << energyerr << endl;
#endif

  double pterr = energyerr*cos(lambda);
  const GsfTrack* gsftrack = &(*(particle->gsfTrack()));
  if (gsftrack!=0){
#if hmc_debug==1
    cout << "HiggsMassConstraint::getCovarianceMatrix(const reco::GsfElectron* particle): GsfTrack " << gsftrack << " found!" << endl;
#endif
    double pterr_uncorrected = gsftrack->ptModeError();
    if (pterr_uncorrected==0.) pterr_uncorrected = pterr;
    double correction = 1.;
    if (pterr_uncorrected!=0.) correction = pow(pterr/pterr_uncorrected, 2);

    double trackCov[GsfTrack::dimensionMode*GsfTrack::dimensionMode];
    for (int ix=0; ix<GsfTrack::dimensionMode; ix++){
      for (int iy=ix; iy<GsfTrack::dimensionMode; iy++){
        if (iy>=ix) trackCov[GsfTrack::dimensionMode*ix+iy] = gsftrack->covarianceMode(ix, iy);
        trackCov[GsfTrack::dimensionMode*iy+ix] = trackCov[GsfTrack::dimensionMode*ix+iy];
      }
    }

    double q = particle->charge();
    double qoverp = q/particle->p();
    double d_pT_d_qoverp;
    if (q==0.){
      q=1.;
      qoverp = q/particle->p();
      d_pT_d_qoverp = 1e12; // (1000 TeV)**2
    }
    else d_pT_d_qoverp = -q*cos(lambda)/pow(qoverp, 2); // ==-p*pT/q
    const double d_pT_d_lambda = -q*sin(lambda)/qoverp; // == -pz
    const double d_pT_d_phi = 0;
    const double d_lambda_d_qoverp = 0.;
    const double d_lambda_d_lambda = 1;
    const double d_lambda_d_phi = 0;
    const double d_phi_d_qoverp = 0;
    const double d_phi_d_lambda = 0;
    const double d_phi_d_phi = 1;

    momCov[3*0+0] = pterr_uncorrected; // pT, pT, no need to re-calculate
    momCov[3*0+1] = // pT, lambda
      d_pT_d_qoverp*d_lambda_d_qoverp * trackCov[GsfTrack::dimensionMode*TrackBase::i_qoverp + TrackBase::i_qoverp] +
      d_pT_d_lambda*d_lambda_d_lambda * trackCov[GsfTrack::dimensionMode*TrackBase::i_lambda + TrackBase::i_lambda] +
      d_pT_d_phi*d_lambda_d_phi * trackCov[GsfTrack::dimensionMode*TrackBase::i_phi + TrackBase::i_phi] +
      (d_pT_d_qoverp*d_lambda_d_lambda + d_pT_d_lambda*d_lambda_d_qoverp) * trackCov[GsfTrack::dimensionMode*TrackBase::i_qoverp + TrackBase::i_lambda] +
      (d_pT_d_qoverp*d_lambda_d_phi + d_pT_d_phi*d_lambda_d_qoverp) * trackCov[GsfTrack::dimensionMode*TrackBase::i_qoverp + TrackBase::i_phi] +
      (d_pT_d_lambda*d_lambda_d_phi + d_pT_d_phi*d_lambda_d_lambda) * trackCov[GsfTrack::dimensionMode*TrackBase::i_lambda + TrackBase::i_phi]
      ;
    momCov[3*0+2] = // pT, phi
      d_pT_d_qoverp*d_phi_d_qoverp * trackCov[GsfTrack::dimensionMode*TrackBase::i_qoverp + TrackBase::i_qoverp] +
      d_pT_d_lambda*d_phi_d_lambda * trackCov[GsfTrack::dimensionMode*TrackBase::i_lambda + TrackBase::i_lambda] +
      d_pT_d_phi*d_phi_d_phi * trackCov[GsfTrack::dimensionMode*TrackBase::i_phi + TrackBase::i_phi] +
      (d_pT_d_qoverp*d_phi_d_lambda + d_pT_d_lambda*d_phi_d_qoverp) * trackCov[GsfTrack::dimensionMode*TrackBase::i_qoverp + TrackBase::i_lambda] +
      (d_pT_d_qoverp*d_phi_d_phi + d_pT_d_phi*d_phi_d_qoverp) * trackCov[GsfTrack::dimensionMode*TrackBase::i_qoverp + TrackBase::i_phi] +
      (d_pT_d_lambda*d_phi_d_phi + d_pT_d_phi*d_phi_d_lambda) * trackCov[GsfTrack::dimensionMode*TrackBase::i_lambda + TrackBase::i_phi]
      ;

    momCov[3*1+0] = momCov[3*0+1];// lambda, pT
    momCov[3*1+1] = // lambda, lambda
      d_lambda_d_qoverp*d_lambda_d_qoverp * trackCov[GsfTrack::dimensionMode*TrackBase::i_qoverp + TrackBase::i_qoverp] +
      d_lambda_d_lambda*d_lambda_d_lambda * trackCov[GsfTrack::dimensionMode*TrackBase::i_lambda + TrackBase::i_lambda] +
      d_lambda_d_phi*d_lambda_d_phi * trackCov[GsfTrack::dimensionMode*TrackBase::i_phi + TrackBase::i_phi] +
      2.*d_lambda_d_qoverp*d_lambda_d_lambda * trackCov[GsfTrack::dimensionMode*TrackBase::i_qoverp + TrackBase::i_lambda] +
      2.*d_lambda_d_qoverp*d_lambda_d_phi * trackCov[GsfTrack::dimensionMode*TrackBase::i_qoverp + TrackBase::i_phi] +
      2.*d_lambda_d_lambda*d_lambda_d_phi * trackCov[GsfTrack::dimensionMode*TrackBase::i_lambda + TrackBase::i_phi]
      ;
    momCov[3*1+2] = // lambda, phi
      d_pT_d_qoverp*d_phi_d_qoverp * trackCov[GsfTrack::dimensionMode*TrackBase::i_qoverp + TrackBase::i_qoverp] +
      d_pT_d_lambda*d_phi_d_lambda * trackCov[GsfTrack::dimensionMode*TrackBase::i_lambda + TrackBase::i_lambda] +
      d_pT_d_phi*d_phi_d_phi * trackCov[GsfTrack::dimensionMode*TrackBase::i_phi + TrackBase::i_phi] +
      (d_pT_d_qoverp*d_phi_d_lambda + d_pT_d_lambda*d_phi_d_qoverp) * trackCov[GsfTrack::dimensionMode*TrackBase::i_qoverp + TrackBase::i_lambda] +
      (d_pT_d_qoverp*d_phi_d_phi + d_pT_d_phi*d_phi_d_qoverp) * trackCov[GsfTrack::dimensionMode*TrackBase::i_qoverp + TrackBase::i_phi] +
      (d_pT_d_lambda*d_phi_d_phi + d_pT_d_phi*d_phi_d_lambda) * trackCov[GsfTrack::dimensionMode*TrackBase::i_lambda + TrackBase::i_phi]
      ;

    momCov[3*2+0] = momCov[3*0+2];// phi, pT
    momCov[3*2+1] = momCov[3*1+2];// phi, lambda
    momCov[3*2+2] = // phi, phi
      d_phi_d_qoverp*d_phi_d_qoverp * trackCov[GsfTrack::dimensionMode*TrackBase::i_qoverp + TrackBase::i_qoverp] +
      d_phi_d_lambda*d_phi_d_lambda * trackCov[GsfTrack::dimensionMode*TrackBase::i_lambda + TrackBase::i_lambda] +
      d_phi_d_phi*d_phi_d_phi * trackCov[GsfTrack::dimensionMode*TrackBase::i_phi + TrackBase::i_phi] +
      2.* d_phi_d_qoverp*d_phi_d_lambda * trackCov[GsfTrack::dimensionMode*TrackBase::i_qoverp + TrackBase::i_lambda] +
      2.* d_phi_d_qoverp*d_phi_d_phi * trackCov[GsfTrack::dimensionMode*TrackBase::i_qoverp + TrackBase::i_phi] +
      2.* d_phi_d_lambda*d_phi_d_phi * trackCov[GsfTrack::dimensionMode*TrackBase::i_lambda + TrackBase::i_phi]
      ;

    // Scale covariance matrix for pT error correction
    for (unsigned int ix=0; ix<3; ix++){
      for (unsigned int iy=0; iy<3; iy++){
        if (ix==iy && ix==0) momCov[3*ix+iy] *= correction;
        else if (ix==0 || iy==0) momCov[3*ix+iy] *= sqrt(correction);
      }
    }
  }
  else{
#if hmc_debug==1
    cout << "HiggsMassConstraint::getCovarianceMatrix(const reco::GsfElectron* particle): No GsfTrack present." << endl;
#endif
    momCov[3*0+0] = pow(pterr, 2);
    // Everything else remains 0.
  }
}
void HiggsMassConstraint::getCovarianceMatrix(double (&momCov)[9], const pat::Muon* particle){
  for(int i=0;i<9;i++) momCov[i]=0.;

  const double lambda = piovertwo_val - particle->theta();

  double pterr_uncorrected = particle->muonBestTrack()->ptError();
  double pterr = pterr_uncorrected;
  if (particle->hasUserFloat("correctedPtError")) pterr = particle->userFloat("correctedPtError");
  double correction = 1.;
  if (pterr_uncorrected!=0.) correction = pow(pterr/pterr_uncorrected, 2);

  double trackCov[TrackBase::dimension*TrackBase::dimension];
  for (int ix=0; ix<TrackBase::dimension; ix++){
    for (int iy=ix; iy<TrackBase::dimension; iy++){
      trackCov[TrackBase::dimension*ix+iy] = particle->muonBestTrack()->covariance(ix, iy);
      trackCov[TrackBase::dimension*iy+ix] = trackCov[TrackBase::dimension*ix+iy];
    }
  }

  double q = particle->charge();
  double qoverp = q/particle->p();
  double d_pT_d_qoverp;
  if (q==0.){
    q=1.;
    qoverp = q/particle->p();
    d_pT_d_qoverp = 1e12; // (1000 TeV)**2
  }
  else d_pT_d_qoverp = -q*cos(lambda)/pow(qoverp, 2); // ==-p*pT/q
  const double d_pT_d_lambda = -q*sin(lambda)/qoverp; // == -pz
  const double d_pT_d_phi = 0;
  const double d_lambda_d_qoverp = 0.;
  const double d_lambda_d_lambda = 1;
  const double d_lambda_d_phi = 0;
  const double d_phi_d_qoverp = 0;
  const double d_phi_d_lambda = 0;
  const double d_phi_d_phi = 1;

  momCov[3*0+0] = pterr_uncorrected; // pT, pT, no need to re-calculate
  momCov[3*0+1] = // pT, lambda
    d_pT_d_qoverp*d_lambda_d_qoverp * trackCov[TrackBase::dimension*TrackBase::i_qoverp + TrackBase::i_qoverp] +
    d_pT_d_lambda*d_lambda_d_lambda * trackCov[TrackBase::dimension*TrackBase::i_lambda + TrackBase::i_lambda] +
    d_pT_d_phi*d_lambda_d_phi * trackCov[TrackBase::dimension*TrackBase::i_phi + TrackBase::i_phi] +
    (d_pT_d_qoverp*d_lambda_d_lambda + d_pT_d_lambda*d_lambda_d_qoverp) * trackCov[TrackBase::dimension*TrackBase::i_qoverp + TrackBase::i_lambda] +
    (d_pT_d_qoverp*d_lambda_d_phi + d_pT_d_phi*d_lambda_d_qoverp) * trackCov[TrackBase::dimension*TrackBase::i_qoverp + TrackBase::i_phi] +
    (d_pT_d_lambda*d_lambda_d_phi + d_pT_d_phi*d_lambda_d_lambda) * trackCov[TrackBase::dimension*TrackBase::i_lambda + TrackBase::i_phi]
    ;
  momCov[3*0+2] = // pT, phi
    d_pT_d_qoverp*d_phi_d_qoverp * trackCov[TrackBase::dimension*TrackBase::i_qoverp + TrackBase::i_qoverp] +
    d_pT_d_lambda*d_phi_d_lambda * trackCov[TrackBase::dimension*TrackBase::i_lambda + TrackBase::i_lambda] +
    d_pT_d_phi*d_phi_d_phi * trackCov[TrackBase::dimension*TrackBase::i_phi + TrackBase::i_phi] +
    (d_pT_d_qoverp*d_phi_d_lambda + d_pT_d_lambda*d_phi_d_qoverp) * trackCov[TrackBase::dimension*TrackBase::i_qoverp + TrackBase::i_lambda] +
    (d_pT_d_qoverp*d_phi_d_phi + d_pT_d_phi*d_phi_d_qoverp) * trackCov[TrackBase::dimension*TrackBase::i_qoverp + TrackBase::i_phi] +
    (d_pT_d_lambda*d_phi_d_phi + d_pT_d_phi*d_phi_d_lambda) * trackCov[TrackBase::dimension*TrackBase::i_lambda + TrackBase::i_phi]
    ;

  momCov[3*1+0] = momCov[3*0+1];// lambda, pT
  momCov[3*1+1] = // lambda, lambda
    d_lambda_d_qoverp*d_lambda_d_qoverp * trackCov[TrackBase::dimension*TrackBase::i_qoverp + TrackBase::i_qoverp] +
    d_lambda_d_lambda*d_lambda_d_lambda * trackCov[TrackBase::dimension*TrackBase::i_lambda + TrackBase::i_lambda] +
    d_lambda_d_phi*d_lambda_d_phi * trackCov[TrackBase::dimension*TrackBase::i_phi + TrackBase::i_phi] +
    2.*d_lambda_d_qoverp*d_lambda_d_lambda * trackCov[TrackBase::dimension*TrackBase::i_qoverp + TrackBase::i_lambda] +
    2.*d_lambda_d_qoverp*d_lambda_d_phi * trackCov[TrackBase::dimension*TrackBase::i_qoverp + TrackBase::i_phi] +
    2.*d_lambda_d_lambda*d_lambda_d_phi * trackCov[TrackBase::dimension*TrackBase::i_lambda + TrackBase::i_phi]
    ;
  momCov[3*1+2] = // lambda, phi
    d_pT_d_qoverp*d_phi_d_qoverp * trackCov[TrackBase::dimension*TrackBase::i_qoverp + TrackBase::i_qoverp] +
    d_pT_d_lambda*d_phi_d_lambda * trackCov[TrackBase::dimension*TrackBase::i_lambda + TrackBase::i_lambda] +
    d_pT_d_phi*d_phi_d_phi * trackCov[TrackBase::dimension*TrackBase::i_phi + TrackBase::i_phi] +
    (d_pT_d_qoverp*d_phi_d_lambda + d_pT_d_lambda*d_phi_d_qoverp) * trackCov[TrackBase::dimension*TrackBase::i_qoverp + TrackBase::i_lambda] +
    (d_pT_d_qoverp*d_phi_d_phi + d_pT_d_phi*d_phi_d_qoverp) * trackCov[TrackBase::dimension*TrackBase::i_qoverp + TrackBase::i_phi] +
    (d_pT_d_lambda*d_phi_d_phi + d_pT_d_phi*d_phi_d_lambda) * trackCov[TrackBase::dimension*TrackBase::i_lambda + TrackBase::i_phi]
    ;

  momCov[3*2+0] = momCov[3*0+2];// phi, pT
  momCov[3*2+1] = momCov[3*1+2];// phi, lambda
  momCov[3*2+2] = // phi, phi
    d_phi_d_qoverp*d_phi_d_qoverp * trackCov[TrackBase::dimension*TrackBase::i_qoverp + TrackBase::i_qoverp] +
    d_phi_d_lambda*d_phi_d_lambda * trackCov[TrackBase::dimension*TrackBase::i_lambda + TrackBase::i_lambda] +
    d_phi_d_phi*d_phi_d_phi * trackCov[TrackBase::dimension*TrackBase::i_phi + TrackBase::i_phi] +
    2.* d_phi_d_qoverp*d_phi_d_lambda * trackCov[TrackBase::dimension*TrackBase::i_qoverp + TrackBase::i_lambda] +
    2.* d_phi_d_qoverp*d_phi_d_phi * trackCov[TrackBase::dimension*TrackBase::i_qoverp + TrackBase::i_phi] +
    2.* d_phi_d_lambda*d_phi_d_phi * trackCov[TrackBase::dimension*TrackBase::i_lambda + TrackBase::i_phi]
    ;

  // Scale covariance matrix for pT error correction
  for (unsigned int ix=0; ix<3; ix++){
    for (unsigned int iy=0; iy<3; iy++){
      if (ix==iy && ix==0) momCov[3*ix+iy] *= correction;
      else if (ix==0 || iy==0) momCov[3*ix+iy] *= sqrt(correction);
    }
  }
}
void HiggsMassConstraint::getCovarianceMatrix(double (&momCov)[9], const reco::PFCandidate* particle){
  for(int i=0;i<9;i++) momCov[i]=0.;

  double lambda = piovertwo_val - particle->theta();

  double energyerr = PFEnergyResolution().getEnergyResolutionEm(particle->energy(), particle->eta());
  double pterr = energyerr*cos(lambda);
  const reco::Track* track = &(*(particle->bestTrack()));
  if (track!=0){
#if hmc_debug==1
    cout << "HiggsMassConstraint::getCovarianceMatrix(const reco::PFCandidate* particle): Track " << track << " found!" << endl;
#endif
    double pterr_uncorrected = track->ptError();
    if (pterr_uncorrected==0.) pterr_uncorrected = pterr;
    double correction = 1.;
    if (pterr_uncorrected!=0.) correction = pow(pterr/pterr_uncorrected, 2);

    double trackCov[TrackBase::dimension*TrackBase::dimension];
    for (int ix=0; ix<TrackBase::dimension; ix++){
      for (int iy=ix; iy<TrackBase::dimension; iy++){
        trackCov[TrackBase::dimension*ix+iy] = track->covariance(ix, iy);
        trackCov[TrackBase::dimension*iy+ix] = trackCov[TrackBase::dimension*ix+iy];
      }
    }

    double q = particle->charge();
    double qoverp = q/particle->p();
    double d_pT_d_qoverp;
    if (q==0.){
      q=1.;
      qoverp = q/particle->p();
      d_pT_d_qoverp = 1e12; // (1000 TeV)**2
    }
    else d_pT_d_qoverp = -q*cos(lambda)/pow(qoverp, 2); // ==-p*pT/q
    const double d_pT_d_lambda = -q*sin(lambda)/qoverp; // == -pz
    const double d_pT_d_phi = 0;
    const double d_lambda_d_qoverp = 0.;
    const double d_lambda_d_lambda = 1;
    const double d_lambda_d_phi = 0;
    const double d_phi_d_qoverp = 0;
    const double d_phi_d_lambda = 0;
    const double d_phi_d_phi = 1;

    momCov[3*0+0] = pterr_uncorrected; // pT, pT, no need to re-calculate
    momCov[3*0+1] = // pT, lambda
      d_pT_d_qoverp*d_lambda_d_qoverp * trackCov[TrackBase::dimension*TrackBase::i_qoverp + TrackBase::i_qoverp] +
      d_pT_d_lambda*d_lambda_d_lambda * trackCov[TrackBase::dimension*TrackBase::i_lambda + TrackBase::i_lambda] +
      d_pT_d_phi*d_lambda_d_phi * trackCov[TrackBase::dimension*TrackBase::i_phi + TrackBase::i_phi] +
      (d_pT_d_qoverp*d_lambda_d_lambda + d_pT_d_lambda*d_lambda_d_qoverp) * trackCov[TrackBase::dimension*TrackBase::i_qoverp + TrackBase::i_lambda] +
      (d_pT_d_qoverp*d_lambda_d_phi + d_pT_d_phi*d_lambda_d_qoverp) * trackCov[TrackBase::dimension*TrackBase::i_qoverp + TrackBase::i_phi] +
      (d_pT_d_lambda*d_lambda_d_phi + d_pT_d_phi*d_lambda_d_lambda) * trackCov[TrackBase::dimension*TrackBase::i_lambda + TrackBase::i_phi]
      ;
    momCov[3*0+2] = // pT, phi
      d_pT_d_qoverp*d_phi_d_qoverp * trackCov[TrackBase::dimension*TrackBase::i_qoverp + TrackBase::i_qoverp] +
      d_pT_d_lambda*d_phi_d_lambda * trackCov[TrackBase::dimension*TrackBase::i_lambda + TrackBase::i_lambda] +
      d_pT_d_phi*d_phi_d_phi * trackCov[TrackBase::dimension*TrackBase::i_phi + TrackBase::i_phi] +
      (d_pT_d_qoverp*d_phi_d_lambda + d_pT_d_lambda*d_phi_d_qoverp) * trackCov[TrackBase::dimension*TrackBase::i_qoverp + TrackBase::i_lambda] +
      (d_pT_d_qoverp*d_phi_d_phi + d_pT_d_phi*d_phi_d_qoverp) * trackCov[TrackBase::dimension*TrackBase::i_qoverp + TrackBase::i_phi] +
      (d_pT_d_lambda*d_phi_d_phi + d_pT_d_phi*d_phi_d_lambda) * trackCov[TrackBase::dimension*TrackBase::i_lambda + TrackBase::i_phi]
      ;

    momCov[3*1+0] = momCov[3*0+1];// lambda, pT
    momCov[3*1+1] = // lambda, lambda
      d_lambda_d_qoverp*d_lambda_d_qoverp * trackCov[TrackBase::dimension*TrackBase::i_qoverp + TrackBase::i_qoverp] +
      d_lambda_d_lambda*d_lambda_d_lambda * trackCov[TrackBase::dimension*TrackBase::i_lambda + TrackBase::i_lambda] +
      d_lambda_d_phi*d_lambda_d_phi * trackCov[TrackBase::dimension*TrackBase::i_phi + TrackBase::i_phi] +
      2.*d_lambda_d_qoverp*d_lambda_d_lambda * trackCov[TrackBase::dimension*TrackBase::i_qoverp + TrackBase::i_lambda] +
      2.*d_lambda_d_qoverp*d_lambda_d_phi * trackCov[TrackBase::dimension*TrackBase::i_qoverp + TrackBase::i_phi] +
      2.*d_lambda_d_lambda*d_lambda_d_phi * trackCov[TrackBase::dimension*TrackBase::i_lambda + TrackBase::i_phi]
      ;
    momCov[3*1+2] = // lambda, phi
      d_pT_d_qoverp*d_phi_d_qoverp * trackCov[TrackBase::dimension*TrackBase::i_qoverp + TrackBase::i_qoverp] +
      d_pT_d_lambda*d_phi_d_lambda * trackCov[TrackBase::dimension*TrackBase::i_lambda + TrackBase::i_lambda] +
      d_pT_d_phi*d_phi_d_phi * trackCov[TrackBase::dimension*TrackBase::i_phi + TrackBase::i_phi] +
      (d_pT_d_qoverp*d_phi_d_lambda + d_pT_d_lambda*d_phi_d_qoverp) * trackCov[TrackBase::dimension*TrackBase::i_qoverp + TrackBase::i_lambda] +
      (d_pT_d_qoverp*d_phi_d_phi + d_pT_d_phi*d_phi_d_qoverp) * trackCov[TrackBase::dimension*TrackBase::i_qoverp + TrackBase::i_phi] +
      (d_pT_d_lambda*d_phi_d_phi + d_pT_d_phi*d_phi_d_lambda) * trackCov[TrackBase::dimension*TrackBase::i_lambda + TrackBase::i_phi]
      ;

    momCov[3*2+0] = momCov[3*0+2];// phi, pT
    momCov[3*2+1] = momCov[3*1+2];// phi, lambda
    momCov[3*2+2] = // phi, phi
      d_phi_d_qoverp*d_phi_d_qoverp * trackCov[TrackBase::dimension*TrackBase::i_qoverp + TrackBase::i_qoverp] +
      d_phi_d_lambda*d_phi_d_lambda * trackCov[TrackBase::dimension*TrackBase::i_lambda + TrackBase::i_lambda] +
      d_phi_d_phi*d_phi_d_phi * trackCov[TrackBase::dimension*TrackBase::i_phi + TrackBase::i_phi] +
      2.* d_phi_d_qoverp*d_phi_d_lambda * trackCov[TrackBase::dimension*TrackBase::i_qoverp + TrackBase::i_lambda] +
      2.* d_phi_d_qoverp*d_phi_d_phi * trackCov[TrackBase::dimension*TrackBase::i_qoverp + TrackBase::i_phi] +
      2.* d_phi_d_lambda*d_phi_d_phi * trackCov[TrackBase::dimension*TrackBase::i_lambda + TrackBase::i_phi]
      ;

    // Scale covariance matrix for pT error correction
    for (unsigned int ix=0; ix<3; ix++){
      for (unsigned int iy=0; iy<3; iy++){
        if (ix==iy && ix==0) momCov[3*ix+iy] *= correction;
        else if (ix==0 || iy==0) momCov[3*ix+iy] *= sqrt(correction);
      }
    }
  }
  else{
#if hmc_debug==1
    cout << "HiggsMassConstraint::getCovarianceMatrix(const reco::PFCandidate* particle): No track present." << endl;
#endif
  momCov[3*0+0] = pow(pterr, 2);
  // Everything else is 0. I know this cannot be correct, but let's work with it for now.
  }
}
void HiggsMassConstraint::getCovarianceMatrix(double (&momCov)[9], const pat::Jet* particle){
  for(int i=0;i<9;i++) momCov[i]=0.;

  double lambda = piovertwo_val - particle->theta();
  double energyerr = particle->userFloat(jecString.Data());

  double C_p_p = energyerr*energyerr;
  // Everything else is 0. I know this cannot be correct, but let's work with it for now.
  // Should loop over track references if the object is PFJet using reco::TrackRefVector PFJet::getTrackRefs() and sum inverse covariance matrices
  momCov[3*0+0] = C_p_p*pow(cos(lambda), 2);
}

void HiggsMassConstraint::invertOneDimensional(Int_t includeIndex, double (&momCov)[9]){
  double momCov_tmp[9]={ 0 };
  if (momCov[3*includeIndex+includeIndex]!=0.) momCov_tmp[3*includeIndex+includeIndex] = 1./momCov[3*includeIndex+includeIndex];
  for (int ix=0; ix<3; ix++){ for (int iy=0; iy<3; iy++) momCov[3*ix+iy]=momCov_tmp[3*ix+iy]; }
}
void HiggsMassConstraint::invertTwoDimensional(Int_t omitIndex, double (&momCov)[9]){
  double momCov_tmp[9]={ 0 };
  double determinant = 0;
  if (omitIndex==0){
    determinant = (momCov[3*1+1]*momCov[3*2+2]) - (momCov[3*1+2]*momCov[3*2+1]);
    momCov_tmp[3*1+1] = momCov[3*2+2];
    momCov_tmp[3*2+2] = momCov[3*1+1];
    momCov_tmp[3*1+2] = -momCov[3*1+2];
    momCov_tmp[3*2+1] = -momCov[3*2+1];
  }
  else if (omitIndex==1){
    determinant = (momCov[3*0+0]*momCov[3*2+2]) - (momCov[3*0+2]*momCov[3*2+0]);
    momCov_tmp[3*0+0] = momCov[3*2+2];
    momCov_tmp[3*2+2] = momCov[3*0+0];
    momCov_tmp[3*0+2] = -momCov[3*0+2];
    momCov_tmp[3*2+0] = -momCov[3*2+0];
  }
  else if (omitIndex==2){
    determinant = (momCov[3*0+0]*momCov[3*1+1]) - (momCov[3*0+1]*momCov[3*1+0]);
    momCov_tmp[3*0+0] = momCov[3*1+1];
    momCov_tmp[3*1+1] = momCov[3*0+0];
    momCov_tmp[3*0+1] = -momCov[3*0+1];
    momCov_tmp[3*1+0] = -momCov[3*1+0];
  }
  for (int ix=0; ix<3; ix++){
    for (int iy=0; iy<3; iy++){
      if (determinant!=0) momCov[3*ix+iy]=momCov_tmp[3*ix+iy]/determinant;
      else momCov[3*ix+iy]=0;
    }
  }
}
void HiggsMassConstraint::invertThreeDimensional(double (&momCov)[9]){
  double momCov_tmp[9]={ 0 };
  double determinant = 0.
    + (momCov[3*0+0]*momCov[3*1+1]*momCov[3*2+2])
    + (momCov[3*0+1]*momCov[3*1+2]*momCov[3*2+0])
    + (momCov[3*0+2]*momCov[3*1+0]*momCov[3*2+1])
    - (momCov[3*0+2]*momCov[3*1+1]*momCov[3*2+0])
    - (momCov[3*0+1]*momCov[3*1+0]*momCov[3*2+2])
    - (momCov[3*0+0]*momCov[3*1+2]*momCov[3*2+1]);
  momCov_tmp[3*0+0] = momCov[3*1+1]*momCov[3*2+2] - momCov[3*1+2]*momCov[3*2+1];
  momCov_tmp[3*1+1] = momCov[3*0+0]*momCov[3*2+2] - momCov[3*0+2]*momCov[3*2+0];
  momCov_tmp[3*2+2] = momCov[3*0+0]*momCov[3*1+1] - momCov[3*0+1]*momCov[3*1+0];
  momCov_tmp[3*0+1] = momCov[3*0+2]*momCov[3*2+1] - momCov[3*0+1]*momCov[3*2+2];
  momCov_tmp[3*1+0] = momCov[3*1+2]*momCov[3*2+0] - momCov[3*2+2]*momCov[3*1+0];
  momCov_tmp[3*1+2] = momCov[3*1+0]*momCov[3*0+2] - momCov[3*0+0]*momCov[3*1+2];
  momCov_tmp[3*2+1] = momCov[3*2+0]*momCov[3*0+1] - momCov[3*2+1]*momCov[3*0+0];
  momCov_tmp[3*2+0] = momCov[3*2+1]*momCov[3*1+0] - momCov[3*1+1]*momCov[3*2+0];
  momCov_tmp[3*0+2] = momCov[3*0+1]*momCov[3*1+2] - momCov[3*0+2]*momCov[3*1+1];
  for (int ix=0; ix<3; ix++){
    for (int iy=0; iy<3; iy++){
      if (determinant!=0) momCov[3*ix+iy]=momCov_tmp[3*ix+iy]/determinant;
      else momCov[3*ix+iy]=0;
    }
  }
}
void HiggsMassConstraint::strategicInvertCovarianceMatrix(Int_t useFullCov, Int_t fitpT, Int_t fitlambda, Int_t fitphi, double (&momCov)[9]){
  // Make sure that the diagonal elements are set to 0 before further computation if useFullCov==0
  if (useFullCov==0){ for (int ix=0; ix<3; ix++){ for (int iy=0; iy<3; iy++){ if(ix!=iy) momCov[3*ix+iy]=0; } } }
  // Fit for only one of the observable types
  if((fitpT+fitlambda+fitphi)==1){
    if(fitpT==1) invertOneDimensional(0, momCov);
    else if(fitlambda==1) invertOneDimensional(1, momCov);
    else if(fitphi==1) invertOneDimensional(2, momCov);
  }
  // Constrain only one variable to initial values
  else if((fitpT+fitlambda+fitphi)==2){
    if(fitpT==0) invertTwoDimensional(0, momCov);
    else if(fitlambda==0) invertTwoDimensional(1, momCov);
    else if(fitphi==0) invertTwoDimensional(2, momCov);
  }
  // Release all three variables
  else invertThreeDimensional(momCov);
}
void HiggsMassConstraint::setInverseCovarianceMatrix(Int_t iZ, Int_t iferm, Int_t fsrindex, Double_t momCov[9]){
  if (fsrindex==0){ for (int i=0; i<9; i++){ invcov_ferm[iZ][iferm][i]->setConstant(false); invcov_ferm[iZ][iferm][i]->setVal(momCov[i]); invcov_ferm[iZ][iferm][i]->setConstant(true); } }
  else{ for (int i=0; i<9; i++){ invcov_fsr[iZ][iferm][i]->setConstant(false); invcov_fsr[iZ][iferm][i]->setVal(momCov[i]); invcov_fsr[iZ][iferm][i]->setConstant(true); } }
#if hmc_debug==1
  if (fsrindex==0){
    cout << "Inverse of the covariance matrix for Z" << iZ+1 << " daughter " << iferm+1 << " is:" << endl;
    for (int i=0; i<3; i++){
      for (int j=0; j<3; j++) cout <<  invcov_ferm[iZ][iferm][3*i+j]->getVal() << '\t';
      cout << endl;
    }
  }
  else{
    cout << "Inverse of the covariance matrix for Z" << iZ+1 << " daughter " << iferm+1 << " FSR is:" << endl;
    for (int i=0; i<3; i++){
      for (int j=0; j<3; j++) cout <<  invcov_fsr[iZ][iferm][3*i+j]->getVal() << '\t';
      cout << endl;
    }
  }
#endif
}

void HiggsMassConstraint::setInitialCovarianceMatrix(Int_t iZ, Int_t iferm, Int_t fsrindex, Double_t momCov[9]){
  const int nDims = 24;
  TMatrixDSym matrix(nDims);
  for (int ix=0; ix<nDims; ix++){
    for (int iy=0; iy<nDims; iy++) matrix[ix][iy] = 0;
  }
  for (int ix=0; ix<3; ix++){
    for (int iy=0; iy<3; iy++){
      matrix[6*(2*iZ+iferm)+3*fsrindex+ix][6*(2*iZ+iferm)+3*fsrindex+iy] = momCov[3*ix+iy];
    }
  }
  initCovMatrix+=matrix;
}
void HiggsMassConstraint::resetInitialCovarianceMatrix(){
  const int nDims = 24;
  initCovMatrix.ResizeTo(nDims, nDims);
  for (int ix=0; ix<nDims; ix++){ for (int iy=0; iy<nDims; iy++) initCovMatrix[ix][iy] = 0; }
}

void HiggsMassConstraint::setInitialMassErrors(){ for (int iv=0; iv<3; iv++) initMassError[iv] = getRefittedMassError(iv); }
void HiggsMassConstraint::resetInitialMassErrors(){ for (int iv=0; iv<3; iv++) initMassError[iv] = 0; }

bool HiggsMassConstraint::standardOrderedFinalCovarianceMatrix(const RooArgList& pars){
  const int nDims = 24;
  TMatrixDSym finalMatrix(nDims);
  for (int ix=0; ix<nDims; ix++){
    for (int iy=0; iy<nDims; iy++) finalMatrix[ix][iy] = 0;
  }

  Int_t nFinalPars = pars.getSize();
  Int_t* order = new Int_t[nFinalPars];
  for (int ip=0; ip<nFinalPars; ip++){
    RooRealVar* arg = dynamic_cast<RooRealVar*>(pars.at(ip));
    Int_t index = fitParameterCorrespondance(arg);
    if (index<0){ delete[] order; return false; }
    order[ip]=index;
  }
  for (int ix=0; ix<nFinalPars; ix++){
    for (int iy=0; iy<nFinalPars; iy++){
      finalMatrix[order[ix]][order[iy]] = fitCovMatrix[ix][iy];
    }
  }
  delete[] order;
  fitCovMatrix.ResizeTo(nDims, nDims);
  fitCovMatrix=finalMatrix;

  // Test the final fit strategy to decide which cov. mat. elements to drop
  Int_t useFullCov, FermFSRType, fitpT, fitlambda, fitphi;
  testFitMomentumStrategy(useFullCov, FermFSRType, fitpT, fitlambda, fitphi);

  // Add the unused initial cov. mat. to calculate the error correctly
  TMatrixDSym addMatrix(nDims);
  addMatrix = initCovMatrix;
  for (int ix=0; ix<nDims; ix++){
    bool doSkip=false;
    for (int iy=ix; iy<nDims; iy++){
      if (fitCovMatrix[ix][iy]!=0.){ doSkip=true; break; }
    }
    if (doSkip){ for (int iy=ix; iy<nDims; iy++){ addMatrix[ix][iy]=0; addMatrix[iy][ix]=0; } }
  }
  for (int ix=0; ix<nDims; ix++){
    for (int iy=ix; iy<nDims; iy++){
      if (useFullCov==0){
        if (ix!=iy){ addMatrix[ix][iy]=0; addMatrix[iy][ix]=0; }
      }
      else{
        // Always include pT errors, decide to include lambda or phi correlations
        if (fitlambda==0 && (ix%3==1 || iy%3==1) && ix!=iy){ addMatrix[ix][iy]=0; addMatrix[iy][ix]=0; }
        if (fitphi==0 && (ix%3==2 || iy%3==2) && ix!=iy){ addMatrix[ix][iy]=0; addMatrix[iy][ix]=0; }
      }
    }
  }
  fitCovMatrix += addMatrix;

#if hmc_debug==1
  cout << "HiggsMassConstraint::standardOrderedFinalCovarianceMatrix: Final covariance matrix is:" << endl;
  for (int ix=0; ix<nDims; ix++){
    for (int iy=0; iy<nDims; iy++) cout << fitCovMatrix[ix][iy] << '\t';
    cout << endl;
  }
#endif
  return true;
}
Int_t HiggsMassConstraint::fitParameterCorrespondance(RooRealVar* par){
  Int_t index=-1;
  for (int iZ=0; iZ<2; iZ++){
    for (int iferm=0; iferm<2; iferm++){
      if (TString(par->GetName())==TString(pT_ferm[iZ][iferm]->GetName())) index = 6*(2*iZ+iferm)+0;
      else if (TString(par->GetName())==TString(lambda_ferm[iZ][iferm]->GetName())) index = 6*(2*iZ+iferm)+1;
      else if (TString(par->GetName())==TString(phi_ferm[iZ][iferm]->GetName())) index = 6*(2*iZ+iferm)+2;
      else if (TString(par->GetName())==TString(pT_fsr[iZ][iferm]->GetName())) index = 6*(2*iZ+iferm)+3;
      else if (TString(par->GetName())==TString(lambda_fsr[iZ][iferm]->GetName())) index = 6*(2*iZ+iferm)+4;
      else if (TString(par->GetName())==TString(phi_fsr[iZ][iferm]->GetName())) index = 6*(2*iZ+iferm)+5;
      if (index>=0) break;
    }
  }
  if (index<0) cerr << "HiggsMassConstraint::fitParameterCorrespondance: Parameter " << par->GetName() << " not found!" << endl;
  return index;
}


Double_t HiggsMassConstraint::d_Ek_d_pTk(Int_t kZ, Int_t kferm, Int_t fsrindex) const{
  Double_t pt_over_E=0;
  Double_t lambda=0;
  Double_t pT=0;
  Double_t E=0;
  if (fsrindex==0){
    E = E_ferm[kZ][kferm]->getVal();
    pT = pT_ferm[kZ][kferm]->getVal();
    lambda = lambda_ferm[kZ][kferm]->getVal();
  }
  else{
    E = E_fsr[kZ][kferm]->getVal();
    pT = pT_fsr[kZ][kferm]->getVal();
    lambda = lambda_fsr[kZ][kferm]->getVal();
  }
  if (E>0.) pt_over_E = pT/E;

  return (pt_over_E/pow(cos(lambda), 2));
}
Double_t HiggsMassConstraint::d_pjk_d_pTk(Int_t kZ, Int_t kferm, Int_t fsrindex, Int_t j) const{
  Double_t lambda=0;
  Double_t phi=0;
  if (fsrindex==0){
    phi = phi_ferm[kZ][kferm]->getVal();
    lambda = lambda_ferm[kZ][kferm]->getVal();
  }
  else{
    phi = phi_fsr[kZ][kferm]->getVal();
    lambda = lambda_fsr[kZ][kferm]->getVal();
  }

  if (j==0) return cos(phi);
  else if (j==1) return sin(phi);
  else return tan(lambda);
}
Double_t HiggsMassConstraint::d_Ek_d_lambdak(Int_t kZ, Int_t kferm, Int_t fsrindex) const{
  Double_t pz=0;
  if (fsrindex==0) pz = pz_ferm[kZ][kferm]->getVal();
  else pz = pz_fsr[kZ][kferm]->getVal();

  return (pz*d_Ek_d_pTk(kZ, kferm, fsrindex));
}
Double_t HiggsMassConstraint::d_pjk_d_lambdak(Int_t kZ, Int_t kferm, Int_t fsrindex, Int_t j) const{
  if (j!=2) return 0; // px and py do not depend on lambda
  Double_t pT=0;
  Double_t lambda=0;
  if (fsrindex==0){
    pT = pT_ferm[kZ][kferm]->getVal();
    lambda = lambda_ferm[kZ][kferm]->getVal();
  }
  else{
    pT = pT_fsr[kZ][kferm]->getVal();
    lambda = lambda_fsr[kZ][kferm]->getVal();
  }
  return (pT/pow(cos(lambda), 2));
}
Double_t HiggsMassConstraint::d_Ek_d_phik(Int_t kZ, Int_t kferm, Int_t fsrindex) const{ return 0; }
Double_t HiggsMassConstraint::d_pjk_d_phik(Int_t kZ, Int_t kferm, Int_t fsrindex, Int_t j) const{
  if (j>=2) return 0;
  Double_t pT=0;
  Double_t phi=0;
  if (fsrindex==0){
    pT = pT_ferm[kZ][kferm]->getVal();
    phi = phi_ferm[kZ][kferm]->getVal();
  }
  else{
    pT = pT_fsr[kZ][kferm]->getVal();
    phi = phi_fsr[kZ][kferm]->getVal();
  }
  if (j==0) return (-pT*sin(phi));
  else return (pT*cos(phi));
}

Double_t HiggsMassConstraint::d_m123_d_pTk(Int_t imass, Int_t kZ, Int_t kferm, Int_t fsrindex) const{
  if((imass<2 && kZ!=imass) || imass<0 || imass>2) return 0;

  Double_t mass=m[imass]->getVal();
  Double_t sum_pj[4]={ 0 };
  for(int iferm=0; iferm<2; iferm++){
    sum_pj[0] += px_Hdaughter[kZ][iferm]->getVal();
    sum_pj[1] += py_Hdaughter[kZ][iferm]->getVal();
    sum_pj[2] += pz_Hdaughter[kZ][iferm]->getVal();
    sum_pj[3] += E_Hdaughter[kZ][iferm]->getVal();
  }
  if(imass==2){
    for(int iferm=0; iferm<2; iferm++){
      sum_pj[0] += px_Hdaughter[1-kZ][iferm]->getVal();
      sum_pj[1] += py_Hdaughter[1-kZ][iferm]->getVal();
      sum_pj[2] += pz_Hdaughter[1-kZ][iferm]->getVal();
      sum_pj[3] += E_Hdaughter[1-kZ][iferm]->getVal();
    }
  }

  Double_t value = d_Ek_d_pTk(kZ, kferm, fsrindex);
  if (mass!=0.){
    value *= sum_pj[3]/mass;
    for (int j=0; j<3; j++) value -= sum_pj[j]/mass*d_pjk_d_pTk(kZ, kferm, fsrindex, j);
  }
  return value;
}
Double_t HiggsMassConstraint::d_m123_d_lambdak(Int_t imass, Int_t kZ, Int_t kferm, Int_t fsrindex) const{
  if((imass<2 && kZ!=imass) || imass<0 || imass>2) return 0;

  Double_t mass=m[imass]->getVal();
  Double_t sum_pj[4]={ 0 };
  for(int iferm=0; iferm<2; iferm++){
    sum_pj[0] += px_Hdaughter[kZ][iferm]->getVal();
    sum_pj[1] += py_Hdaughter[kZ][iferm]->getVal();
    sum_pj[2] += pz_Hdaughter[kZ][iferm]->getVal();
    sum_pj[3] += E_Hdaughter[kZ][iferm]->getVal();
  }
  if(imass==2){
    for(int iferm=0; iferm<2; iferm++){
      sum_pj[0] += px_Hdaughter[1-kZ][iferm]->getVal();
      sum_pj[1] += py_Hdaughter[1-kZ][iferm]->getVal();
      sum_pj[2] += pz_Hdaughter[1-kZ][iferm]->getVal();
      sum_pj[3] += E_Hdaughter[1-kZ][iferm]->getVal();
    }
  }

  Double_t value = d_Ek_d_lambdak(kZ, kferm, fsrindex);
  if (mass!=0.){
    value *= sum_pj[3]/mass;
    for (int j=0; j<3; j++) value -= sum_pj[j]/mass*d_pjk_d_lambdak(kZ, kferm, fsrindex, j);
  }
  return value;
}
Double_t HiggsMassConstraint::d_m123_d_phik(Int_t imass, Int_t kZ, Int_t kferm, Int_t fsrindex) const{
  if((imass<2 && kZ!=imass) || imass<0 || imass>2) return 0;

  Double_t mass=m[imass]->getVal();
  Double_t sum_pj[4]={ 0 };
  for(int iferm=0; iferm<2; iferm++){
    sum_pj[0] += px_Hdaughter[kZ][iferm]->getVal();
    sum_pj[1] += py_Hdaughter[kZ][iferm]->getVal();
    sum_pj[2] += pz_Hdaughter[kZ][iferm]->getVal();
    sum_pj[3] += E_Hdaughter[kZ][iferm]->getVal();
  }
  if(imass==2){
    for(int iferm=0; iferm<2; iferm++){
      sum_pj[0] += px_Hdaughter[1-kZ][iferm]->getVal();
      sum_pj[1] += py_Hdaughter[1-kZ][iferm]->getVal();
      sum_pj[2] += pz_Hdaughter[1-kZ][iferm]->getVal();
      sum_pj[3] += E_Hdaughter[1-kZ][iferm]->getVal();
    }
  }

  Double_t value = d_Ek_d_phik(kZ, kferm, fsrindex);
  if (mass!=0.){
    value *= sum_pj[3]/mass;
    for (int j=0; j<3; j++) value -= sum_pj[j]/mass*d_pjk_d_phik(kZ, kferm, fsrindex, j);
  }
  return value;
}

Double_t HiggsMassConstraint::getRefittedMassError(Int_t imass) const{ // imass==0 is m1, imass==1 is m2, imass==2 is m12.
  Double_t value = 0;
  const Int_t NColCovMat=24;
  if (fitCovMatrix.GetNcols()!=NColCovMat) return value; // This means the fit was not run.
  vector<Double_t> jacArray;
  for (int iZ=0; iZ<2; iZ++){
    for (int iferm=0; iferm<2; iferm++){
      for (int ifsr=0; ifsr<2; ifsr++){
        Double_t d_m_d_pTk_ival = d_m123_d_pTk(imass, iZ, iferm, ifsr);
        Double_t d_m_d_lambdak_ival = d_m123_d_lambdak(imass, iZ, iferm, ifsr);
        Double_t d_m_d_phik_ival = d_m123_d_phik(imass, iZ, iferm, ifsr);

        jacArray.push_back(d_m_d_pTk_ival);
        jacArray.push_back(d_m_d_lambdak_ival);
        jacArray.push_back(d_m_d_phik_ival);
      }
    }
  }

  for (int iZ=0; iZ<2; iZ++){
    for (int iferm=0; iferm<2; iferm++){
      for (int ifsr=0; ifsr<2; ifsr++){
        Int_t ipt = 6*(2*iZ+iferm)+3*ifsr+0;
        Int_t ilambda = 6*(2*iZ+iferm)+3*ifsr+1;
        Int_t iphi = 6*(2*iZ+iferm)+3*ifsr+2;

        for (int jZ=0; jZ<2; jZ++){
          for (int jferm=0; jferm<2; jferm++){
            for (int jfsr=0; jfsr<2; jfsr++){
              Int_t jpt = 6*(2*jZ+jferm)+3*jfsr+0;
              Int_t jlambda = 6*(2*jZ+jferm)+3*jfsr+1;
              Int_t jphi = 6*(2*jZ+jferm)+3*jfsr+2;

              value += (fitCovMatrix[ipt][jpt]) * jacArray.at(ipt) * jacArray.at(jpt);
              value += (fitCovMatrix[ipt][jlambda]) * jacArray.at(ipt) * jacArray.at(jlambda);
              value += (fitCovMatrix[ipt][jphi]) * jacArray.at(ipt) * jacArray.at(jphi);
              value += (fitCovMatrix[ilambda][jpt]) * jacArray.at(ilambda) * jacArray.at(jpt);
              value += (fitCovMatrix[ilambda][jlambda]) * jacArray.at(ilambda) * jacArray.at(jlambda);
              value += (fitCovMatrix[ilambda][jphi]) * jacArray.at(ilambda) * jacArray.at(jphi);
              value += (fitCovMatrix[iphi][jpt]) * jacArray.at(iphi) * jacArray.at(jpt);
              value += (fitCovMatrix[iphi][jlambda]) * jacArray.at(iphi) * jacArray.at(jlambda);
              value += (fitCovMatrix[iphi][jphi]) * jacArray.at(iphi) * jacArray.at(jphi);
            }
          }
        }
      }
    }
  }

  if (value<=0.){
    cout << "HiggsMassConstraint::getRefittedMassError: m" << imass+1 << " error " << value << " is non-positive. Covariance matrix and derivatives used were:" << endl;
    summarizeDaughters();
    for (int ix=0; ix<24; ix++){
      for (int iy=0; iy<24; iy++) cout << fitCovMatrix[ix][iy] << '\t';
      cout << endl;
    }
    for (int iZ=0; iZ<2; iZ++){
      for (int iferm=0; iferm<2; iferm++){
        for (int ifsr=0; ifsr<2; ifsr++){
          Int_t ipt = 6*(2*iZ+iferm)+3*ifsr+0;
          Int_t ilambda = 6*(2*iZ+iferm)+3*ifsr+1;
          Int_t iphi = 6*(2*iZ+iferm)+3*ifsr+2;

          cout << "d_m" << imass+1 << "_d_pTk(" << iZ << "," << iferm << "," << ifsr << ") = " << jacArray.at(ipt) << endl;
          cout << "d_m" << imass+1 << "_d_lambdak(" << iZ << "," << iferm << "," << ifsr << ") = " << jacArray.at(ilambda) << endl;
          cout << "d_m" << imass+1 << "_d_phik(" << iZ << "," << iferm << "," << ifsr << ") = " << jacArray.at(iphi) << endl;
        }
      }
    }
    value=0;
    //assert(0);
  }

  value = sqrt(value);
  return value;
}
Double_t HiggsMassConstraint::getRefittedMass(Int_t imass) const{
  if (imass<3) return m[imass]->getVal();
  else return 0;
}
TLorentzVector HiggsMassConstraint::getRefittedMomentum(Int_t iZ, Int_t iferm, Int_t fsrindex) const{
  TLorentzVector result;
  if (fsrindex==0) result.SetXYZT(px_ferm[iZ][iferm]->getVal(), py_ferm[iZ][iferm]->getVal(), pz_ferm[iZ][iferm]->getVal(), E_ferm[iZ][iferm]->getVal());
  else result.SetXYZT(px_fsr[iZ][iferm]->getVal(),py_fsr[iZ][iferm]->getVal(),pz_fsr[iZ][iferm]->getVal(),E_fsr[iZ][iferm]->getVal());
  return result;
}
Double_t HiggsMassConstraint::getRefittedPtError(Int_t iZ, Int_t iferm, Int_t fsrindex) const{
  Int_t ipt = 6*(2*iZ+iferm)+3*fsrindex+0;
  Double_t res = fitCovMatrix[ipt][ipt]; if (std::isnan(res)) res=0;
  return sqrt(res);
}
Double_t HiggsMassConstraint::getRefittedLambdaError(Int_t iZ, Int_t iferm, Int_t fsrindex) const{
  Int_t ilambda = 6*(2*iZ+iferm)+3*fsrindex+1;
  Double_t res = fitCovMatrix[ilambda][ilambda]; if (std::isnan(res)) res=0;
  return sqrt(res);
}
Double_t HiggsMassConstraint::getRefittedPhiError(Int_t iZ, Int_t iferm, Int_t fsrindex) const{
  Int_t iphi = 6*(2*iZ+iferm)+3*fsrindex+2;
  Double_t res = fitCovMatrix[iphi][iphi]; if (std::isnan(res)) res=0;
  return sqrt(res);
}
Double_t HiggsMassConstraint::getObsMassError(Int_t imass) const{ // imass==0 is m1, imass==1 is m2, imass==2 is m12.
  if (imass<3) return initMassError[imass];
  else return 0;
}
Double_t HiggsMassConstraint::getObsPtError(Int_t iZ, Int_t iferm, Int_t fsrindex) const{
  Int_t ipt = 6*(2*iZ+iferm)+3*fsrindex+0;
  Double_t res = initCovMatrix[ipt][ipt]; if (std::isnan(res)) res=0;
  return sqrt(res);
}
Double_t HiggsMassConstraint::getObsLambdaError(Int_t iZ, Int_t iferm, Int_t fsrindex) const{
  Int_t ilambda = 6*(2*iZ+iferm)+3*fsrindex+1;
  Double_t res = initCovMatrix[ilambda][ilambda]; if (std::isnan(res)) res=0;
  return sqrt(res);
}
Double_t HiggsMassConstraint::getObsPhiError(Int_t iZ, Int_t iferm, Int_t fsrindex) const{
  Int_t iphi = 6*(2*iZ+iferm)+3*fsrindex+2;
  Double_t res = initCovMatrix[iphi][iphi]; if (std::isnan(res)) res=0;
  return sqrt(res);
}

