#include "../util/RootStyle.cc"

// Defining all values needed
const TString stream = "_offline";
const TString bkg[] = {
  "QCD_HT200to300","QCD_HT300to500","QCD_HT500to700","QCD_HT700to1000",
  "QCD_HT1000to1500","QCD_HT1500to2000","QCD_HT2000toInf"
};
const Float_t xs[] = {1559000,311900,29070,5962,1207,119.9,25.24};

TH1F* plotQCD() {
  unsigned int n_bkg = sizeof(bkg)/sizeof(bkg[0]);
  TFile* f_bkg[n_bkg];
  TH1F* hb[n_bkg];
  TH1F* h_bkg = new TH1F("h_bkg","QCD sphericity",100,0,1);
  for (int i = 0; i < n_bkg; ++i) {
    cout << "output/"+bkg[i]+".root" << "\n";
    f_bkg[i] = TFile::Open("output/"+bkg[i]+".root");
    hb[i] = (TH1F*)f_bkg[i]->Get(bkg[i]+stream+"_evtshape_sphericity");
    h_bkg->Add(hb[i],xs[i]);
  }
  return h_bkg;
}


void plotSphericity() {
  TString s1 = "mMed-125_mDark-2_temp-2_decay-generic";
  TString s2 = "mMed-400_mDark-2_temp-2_decay-generic";
  TString s3 = "mMed-750_mDark-2_temp-2_decay-generic";
  TString s4 = "mMed-1000_mDark-2_temp-2_decay-generic";
  TFile* f1 = TFile::Open("output/"+s1+".root");
  TFile* f2 = TFile::Open("output/"+s2+".root");
  TFile* f3 = TFile::Open("output/"+s3+".root");
  TFile* f4 = TFile::Open("output/"+s4+".root");
  TH1F* h1 = (TH1F*)f1->Get(s1+stream+"_evtshape_sphericity");
  TH1F* h2 = (TH1F*)f2->Get(s2+stream+"_evtshape_sphericity");
  TH1F* h3 = (TH1F*)f3->Get(s3+stream+"_evtshape_sphericity");
  TH1F* h4 = (TH1F*)f4->Get(s4+stream+"_evtshape_sphericity");
  h1->Scale(5./h1->Integral());
  h2->Scale(10./h2->Integral());
  h3->Scale(10./h3->Integral());
  h4->Scale(10./h4->Integral());

  set_root_style();
  TCanvas* c = new TCanvas("c","c",1);
  THStack* hs = new THStack();
  hs->Add(h1);
  hs->Add(h2);
  hs->Add(h3);
  hs->Add(h4);
  h1->SetLineWidth(2);
  h2->SetLineWidth(2);
  h3->SetLineWidth(2);
  h4->SetLineWidth(2);
  hs->Draw("plc hist nostack");
  hs->SetTitle(";sphericity;");

  TLegend* lg = new TLegend(0.3,0.65,0.89,0.89);
  lg->AddEntry(h1,s1,"l");
  lg->AddEntry(h2,s2,"l");
  lg->AddEntry(h3,s3,"l");
  lg->AddEntry(h4,s4,"l");
  lg->SetBorderSize(0);
  lg->SetTextSize(0.032);
  lg->Draw();

  TH1F* hQCD = plotQCD();
  hQCD->Scale(10./hQCD->Integral());
  hQCD->SetLineWidth(2);
  hQCD->Draw("hist same");
  TLegend* lg_bkg = new TLegend(0.3,0.6,0.89,0.65);
  lg_bkg->AddEntry(hQCD,"QCD ","l");
  lg_bkg->SetBorderSize(0);
  lg_bkg->SetTextSize(0.032);
  lg_bkg->Draw();
}
