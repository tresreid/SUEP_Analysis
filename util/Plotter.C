TH1F* plotQCD(const TString variable, const TString path, const int n_bkg, const TString bkg[], const TString stream, const Float_t xs[]) {
  TFile* f_bkg[n_bkg];
  TH1F* hb[n_bkg];
  TH1F* h_bkg = new TH1F("h_bkg","QCD "+variable,100,0,1);
  for (int i = 0; i < n_bkg; ++i) {
    f_bkg[i] = TFile::Open(path+"/"+bkg[i]+".root");
    hb[i] = (TH1F*)f_bkg[i]->Get(bkg[i]+"_"+stream+"_evtshape_"+variable);
    h_bkg->Add(hb[i],xs[i]);
  }
  return h_bkg;
}

void Plotter(const TString variable, const TString model, const TString path, const int n_sgnl, const TString sgnl[],
                    const int n_bkg, const TString bkg[], const TString stream,
                    const Float_t scale[], const Float_t xs[]) {
  set_root_style();

  TCanvas* c = new TCanvas(model,model,1);
  TFile* fs[n_sgnl];
  TH1F* hs[n_sgnl];
  THStack* h_signal = new THStack();
  TH1F* hQCD = plotQCD(variable, path, n_bkg, bkg, stream, xs);
  TLegend* lg_s = new TLegend(0.27,0.75,0.87,0.92);
  TLegend* lg_bkg = new TLegend(0.27,0.7,0.87,0.75);
  TLatex ltx;

  for (int i = 0; i < n_sgnl; ++i) {
    fs[i] = TFile::Open(path+"/"+sgnl[i]+"_decay-"+model+".root");
    hs[i] = (TH1F*)fs[i]->Get(sgnl[i]+"_decay-"+model+"_"+stream+"_evtshape_"+variable);
    hs[i]->Scale(scale[i]/hs[i]->Integral());
    h_signal->Add(hs[i]);
    hs[i]->SetLineWidth(2);
    lg_s->AddEntry(hs[i],sgnl[i]+"_decay-"+model,"l");
  }

  h_signal->Draw("plc hist nostack");
  h_signal->SetTitle(";"+variable+";");
  lg_s->SetBorderSize(0);
  lg_s->SetTextSize(0.03);
  lg_s->Draw();
  hQCD->Scale(scale[n_sgnl]/hQCD->Integral());
  hQCD->SetLineWidth(2);
  hQCD->Draw("hist same");
  lg_bkg->AddEntry(hQCD,"QCD ","l");
  lg_bkg->SetBorderSize(0);
  lg_bkg->SetTextSize(0.03);
  lg_bkg->Draw();

  ltx.SetTextSize(0.028);
  if(stream == "offline") ltx.DrawLatexNDC(0.4,0.67,"Offline (H_{T}>1200GeV or leading jet p_{t}>500GeV)");
  else if(stream == "scouting") ltx.DrawLatexNDC(0.4,0.67,"Scouting stream (H_{T}>500GeV)");
  c->SaveAs("plots/"+variable+"_decay-"+model+"_"+stream+".png");
}
