void plotSphericity() {
  TString s1 = "mMed-125_mDark-2_temp-2_decay-generic";
  TString s2 = "mMed-400_mDark-2_temp-2_decay-generic";
  TString s3 = "mMed-750_mDark-2_temp-2_decay-generic";
  TString s4 = "mMed-1000_mDark-2_temp-2_decay-generic";
  TFile* f1 = TFile::Open("output/"+s1+".root");
  TFile* f2 = TFile::Open("output/"+s2+".root");
  TFile* f3 = TFile::Open("output/"+s3+".root");
  TFile* f4 = TFile::Open("output/"+s4+".root");
  TH1F* h1 = (TH1F*)f1->Get(s1+"_offline_evtshape_sphericity");
  TH1F* h2 = (TH1F*)f2->Get(s2+"_offline_evtshape_sphericity");
  TH1F* h3 = (TH1F*)f3->Get(s3+"_offline_evtshape_sphericity");
  TH1F* h4 = (TH1F*)f4->Get(s4+"_offline_evtshape_sphericity");

  TCanvas* c = new TCanvas("c","c",1);
  THStack* hs = new THStack();
  hs->Add(h1);
  hs->Add(h2);
  hs->Add(h3);
  hs->Add(h4);
  hs->Draw("plc hist nostack")
  c->BuildLegend();
}
