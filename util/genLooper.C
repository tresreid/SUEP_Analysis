{
  TStopwatch sw;
  sw.Start();
  TFile *f = TFile::Open("/Users/chrispap/PrivateSamples.SUEP_2018_mMed-1000_mDark-2_temp-2_decay-darkPho_13TeV-pythia8_n-100_0_RA2AnalysisTree.root");
  ROOT::RDataFrame d("TreeMaker2/PreSelection",f);
  auto getX = [](vector<TVector3> &v){RVec<double> out; for( auto p: v) out.push_back(p.X()); return out;};
  auto h2=d.Define("X",getX,{"vec_list"}).Histo1D("X")
  h2->Draw();
  //auto df = d.Filter("GenParticles.size() > 0").Define("Pt","GenParticles.pT()");
  //auto h = df.Histo1D("Pt");
  //h->Draw();
  sw.Stop();
  sw.Print();
}
