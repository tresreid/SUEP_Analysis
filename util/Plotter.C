TH1F* plotQCD(const TString variable, const TString path, const int n_bkg,
              const TString bkg[], const TString stream, const Float_t xs[],
              const int bins, const float xmin, const float xmax) {
  TFile* f_bkg[n_bkg];
  TH1F* hb[n_bkg];
  TH1F* h_bkg = new TH1F("h_bkg","QCD "+variable,bins,xmin,xmax);
  for (int i = 0; i < n_bkg; ++i) {
    f_bkg[i] = TFile::Open(path+"/"+bkg[i]+".root");
    hb[i] = (TH1F*)f_bkg[i]->Get(bkg[i]+"_"+stream+"_"+variable);
    h_bkg->Add(hb[i],xs[i]);
  }
  return h_bkg;
}

void Plotter(const TString variable, const TString model, const TString path, const int n_sgnl, const TString sgnl[],
                    const int n_bkg, const TString bkg[], const TString stream,
                    const Float_t scale[], const Float_t xs[]) {
  set_root_style();
  int rebin_ratio = 1;

  TCanvas* c = new TCanvas(model,model,1);
  TFile* fs[n_sgnl];
  TH1F* hs[n_sgnl];
  THStack* h_signal = new THStack();
  TH1F* hQCD;
  TLegend* lg_s = new TLegend(0.27,0.75,0.87,0.92);
  TLegend* lg_bkg = new TLegend(0.27,0.7,0.87,0.75);
  TLatex ltx;

  for (int i = 0; i < n_sgnl; ++i) {
    fs[i] = TFile::Open(path+"/"+sgnl[i]+"_decay-"+model+".root");
    hs[i] = (TH1F*)fs[i]->Get(sgnl[i]+"_decay-"+model+"_"+stream+"_"+variable);
    hs[i]->Scale(scale[i]/hs[i]->Integral());
    hs[i]->Rebin(rebin_ratio);
    h_signal->Add(hs[i]);
    hs[i]->SetLineWidth(3);
    lg_s->AddEntry(hs[i],sgnl[i]+"_decay-"+model,"l");
  }

  h_signal->Draw("plc hist nostack");
  h_signal->SetTitle(";"+variable+";");
  lg_s->SetBorderSize(0);
  lg_s->SetTextSize(0.03);
  lg_s->Draw();

  hQCD = plotQCD(variable, path, n_bkg, bkg, stream, xs, hs[0]->GetNbinsX(), hs[0]->GetXaxis()->GetXmin(), hs[0]->GetXaxis()->GetXmax());
  hQCD->Scale(scale[n_sgnl]/hQCD->Integral());
  hQCD->Rebin(rebin_ratio);
  hQCD->SetLineWidth(3);
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
