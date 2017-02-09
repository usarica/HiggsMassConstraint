{

  TString LIBMCFM="libmcfm_703.so";
  gSystem->AddIncludePath("-I$ROOFITSYS/include/");
  gSystem->AddIncludePath("-I$CMSSW_BASE/interface/");
  gSystem->AddIncludePath("-I$CMSSW_BASE/src/");
  gSystem->AddIncludePath("-I$CMSSW_BASE/src/ZZMatrixElement/MELA/interface/");
  gSystem->AddIncludePath("-I$CMSSW_BASE/src/HiggsMassConstraint/HiggsMassConstraint/interface/");
  gSystem->Load("libZZMatrixElementMELA.so");
  gSystem->Load("libHiggsMassConstraintHiggsMassConstraint.so");
  gSystem->Load("$CMSSW_BASE/src/ZZMatrixElement/MELA/data/$SCRAM_ARCH/" + LIBMCFM);

}