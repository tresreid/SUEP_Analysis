void Plotter2(const TString variable, const TString model, const TString path, const int n_sgnl, const TString sgnl[],
                    const int n_bkg, const TString bkg[], const TString stream,
                    const Float_t scale[], const Float_t xs[]) {
  set_root_style();
  int rebin_ratio = 1;

  TCanvas* c = new TCanvas(model,model,1);
  TFile* fs[n_sgnl];
  TTree* tr[n_sgnl];
  TH1F* hs[n_sgnl];
  THStack* h_signal = new THStack();
  TLegend* lg_s = new TLegend(0.27,0.75,0.87,0.92);

  for (int i = 0; i < n_sgnl; ++i) {
    fs[i] = TFile::Open(path+"/PrivateSamples.SUEP_2018_"+sgnl[i]+"_decay-"+model+"_13TeV-pythia8_n-100_0_RA2AnalysisTree.root");
    tr[i] = (TTree*)fs[i]->Get("TreeMaker2/PreSelection");
    hs[i] = (TH1F*)tr[i]->Draw("GenParticles.Pt()");
    h_signal->Add(hs[i]);
    hs[i]->SetLineWidth(3);
    lg_s->AddEntry(hs[i],sgnl[i]+"_decay-"+model,"l");
  }

  h_signal->Draw("plc hist nostack");
  h_signal->SetTitle(";"+variable+";");
  lg_s->SetBorderSize(0);
  lg_s->SetTextSize(0.03);
  lg_s->Draw();
}
