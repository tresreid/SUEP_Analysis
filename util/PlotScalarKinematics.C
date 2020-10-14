#include "../util/RootStyle.cc"
#include "../util/Plotter2.C"

void PlotScalarKinematics(){
  // Defining all values needed
  const TString path = "root://cmseos.fnal.gov//store/user/kdipetri/SUEP/Production_v0.1/2018/NTUP/";
  // or locally
  // const TString path = "output";
  const TString variable[] = {""};
    //GenParticles.Pt()
  //};
  const TString model[] = {"darkPho","darkPhoHad"};
  const TString stream[] = {""};//{"offline","scouting"};
  const TString sgnl[] = {
    "mMed-125_mDark-2_temp-2","mMed-400_mDark-2_temp-2",
    "mMed-750_mDark-2_temp-2","mMed-1000_mDark-2_temp-2"
  };
  const TString bkg[] = {""};
    //"QCD_HT200to300","QCD_HT300to500","QCD_HT500to700","QCD_HT700to1000",
    //"QCD_HT1000to1500","QCD_HT1500to2000","QCD_HT2000toInf"
  //};
  const Float_t xs[] = {1559000,311900,29070,5962,1207,119.9,25.24};
  const Float_t scale[] = {1,1,1,1,1};//{5.,10.,10.,10.,10.};

  const int n_sgnl = sizeof(sgnl)/sizeof(sgnl[0]);
  const int n_bkg = sizeof(bkg)/sizeof(bkg[0]);

  for (int i = 0; i < 0; ++i) {
    Plotter2(variable[i], model[0], path, n_sgnl, sgnl, n_bkg, bkg, stream[0], scale, xs);
    Plotter2(variable[i], model[1], path, n_sgnl, sgnl, n_bkg, bkg, stream[0], scale, xs);
  }
}
