#include "TROOT.h"
#include "omp.h"
#include "TTree.h"
#include "TFile.h"
#include "TBranch.h"
#include "Math/PtEtaPhiE4D.h"
#include "Math/LorentzVector.h"
#include "Math/Vector4Dfwd.h"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"
#include <cmath>
#include <tuple>
#include <cstdlib>
#include <iostream>
#include <stdio.h>
#include "TLorentzVector.h"
using namespace fastjet;
using namespace std;

#include <vector>

tuple<double,double> jet_width(PseudoJet jet){
  float girth = 0;
  float pt_avg =0;
  int count =0;
  vector<PseudoJet> constituents = jet.constituents();
  for(int i=0; i< constituents.size(); i++){
    float phi = constituents[i].phi();
    float eta = constituents[i].eta();
    float pt = constituents[i].pt();
    float dPhi = abs(jet.phi() - phi);
    if(dPhi > M_PI){dPhi = dPhi - 2*M_PI;}
    float dEta = jet.eta() - eta;
    float dR = sqrt(dEta*dEta + dPhi*dPhi);
    girth += pt * dR;
    pt_avg += pt;
    count += 1;
  }
  return {(girth / jet.pt()),((pt_avg/count)/jet.pt())};



}



tuple<int,int> jet_constituents(PseudoJet jet, vector<ROOT::Math::PtEtaPhiEVector> scalars, vector<ROOT::Math::PtEtaPhiEVector> isrs){
    int scalar_part = 0;
    int isr_part = 0;
    vector<PseudoJet> constituents = jet.constituents();
    for(int i=0; i< constituents.size(); i++){
    int min_dR = 9999;
    bool is_scalar = true;
    float phi = constituents[i].phi();
    float eta = constituents[i].eta();
    float pt = constituents[i].pt();
    for(int j=0; j<scalars.size();j++){
      if(abs(scalars[j].pt() - pt)/scalars[j].pt() > 0.3) {continue;}
      float dPhi = abs(scalars[j].phi() - phi);
      if(dPhi > M_PI){dPhi = dPhi - 2*M_PI;}
      float dEta = scalars[j].eta() - eta;
      float dR = sqrt(dEta*dEta + dPhi*dPhi);
      if(dR < min_dR){
        is_scalar = true;
        min_dR = dR;
      }
    }
    for(int j=0; j<isrs.size();j++){
      if(abs(isrs[j].pt() - pt)/isrs[j].pt() > 0.3) {continue;}
      float dPhi = abs(isrs[j].phi() - phi);
      if(dPhi > M_PI){dPhi = dPhi - 2*M_PI;}
      float dEta = isrs[j].eta() - eta;
      float dR = sqrt(dEta*dEta + dPhi*dPhi);
      if(dR < min_dR){
        is_scalar = false;
        min_dR = dR;
      }
    }
    if(is_scalar){ scalar_part++;}
    else {isr_part++;}
  }
  return {scalar_part, isr_part};
}

    //run different clustering algorithms
inline JetDefinition getJetDef(int algo, float Ri){
    if(algo == 0){JetDefinition jet_def(cambridge_algorithm,Ri);
      return jet_def;
    }
    if(algo == 1){JetDefinition jet_def(kt_algorithm,Ri);
      return jet_def;
    }
    else{JetDefinition jet_def(antikt_algorithm,Ri);
      return jet_def;
    }
}

int main(int argc,char** argv){
  printf("Starting\n");
  ROOT::EnableThreadSafety();  
  ROOT::EnableImplicitMT();
  //TFile *infile = TFile::Open("root://cmseos.fnal.gov//store/user/kdipetri/SUEP/Production_v0.2/2018/NTUP/PrivateSamples.SUEP_2018_mMed-1000_mDark-2_temp-2_decay-darkPhoHad_13TeV-pythia8_n-100_0_RA2AnalysisTree.root");
  string file = "root://cmseos.fnal.gov//store/user/kdipetri/SUEP/Production_v0.2/2018/NTUP/";
  int data=0;
  vector<int> skip;
  string sample;
  if(argc!=1){data = atoi(argv[1]);}
  switch(data){
  case 0: file += "PrivateSamples.SUEP_2018_mMed-1000_mDark-2_temp-2_decay-darkPhoHad_13TeV-pythia8_n-100_0_RA2AnalysisTree.root";sample="sig_1000"; break;
  case 1: file += "PrivateSamples.SUEP_2018_mMed-750_mDark-2_temp-2_decay-darkPhoHad_13TeV-pythia8_n-100_0_RA2AnalysisTree.root"; sample="sig_750";break;
  case 2: file += "PrivateSamples.SUEP_2018_mMed-400_mDark-2_temp-2_decay-darkPhoHad_13TeV-pythia8_n-100_0_RA2AnalysisTree.root"; sample="sig_400";break;
  case 3: file += "Autumn18.QCD_HT300to500_TuneCP5_13TeV-madgraphMLM-pythia8_RA2AnalysisTree.root";   sample="qcd_300";break;
  case 4: file += "Autumn18.QCD_HT500to700_TuneCP5_13TeV-madgraphMLM-pythia8_RA2AnalysisTree.root";   sample="qcd_500";break;
  case 5: file += "Autumn18.QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8_RA2AnalysisTree.root";  sample="qcd_700";break;
  case 6: file += "Autumn18.QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8_RA2AnalysisTree.root"; sample="qcd_1000"; break;
  case 7: file += "Autumn18.QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8_RA2AnalysisTree.root"; sample="qcd_1500";break;
  case 8: file += "Autumn18.QCD_HT2000toInf_TuneCP5_13TeV-madgraphMLM-pythia8_RA2AnalysisTree.root";  sample="qcd_2000"; break;
  default: file += "PrivateSamples.SUEP_2018_mMed-1000_mDark-2_temp-2_decay-darkPhoHad_13TeV-pythia8_n-100_0_RA2AnalysisTree.root"; sample="sig_1000";break;
  }
  int batch = 0;
  if(argc>=3){batch = atoi(argv[2]);}
  //printf("opening file: %s\n",file.c_str());
  FILE* OutFile = fopen(("data/"+sample+"_v"+to_string(batch)+".txt").c_str(),"w");
  TFile *infile = TFile::Open(file.c_str());

  TTree *t1= (TTree*)infile->Get("TreeMaker2/PreSelection");
  vector<ROOT::Math::PtEtaPhiEVector>* genPart = 0;
  vector<ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>>>* trkPart =0;
  double ht;
  vector<int>* gen_parentId=0;
  vector<int>* gen_charge =0;
  vector<int>* gen_pdgId =0;
  vector<int>* gen_status=0;
  vector<int>* trk_pv =0;
  vector<bool>* trk_matched =0;
  int nvtx =0;
  int numInteractions =0;
  t1->SetBranchAddress("GenParticles",&genPart);
  t1->SetBranchAddress("Tracks",&trkPart);
  t1->SetBranchAddress("HT",&ht);
  t1->SetBranchAddress("GenParticles_ParentId",&gen_parentId);
  t1->SetBranchAddress("GenParticles_Charge",&gen_charge);
  t1->SetBranchAddress("GenParticles_PdgId",&gen_pdgId);
  t1->SetBranchAddress("GenParticles_Status",&gen_status);
  t1->SetBranchAddress("Tracks_fromPV0",&trk_pv);
  t1->SetBranchAddress("Tracks_matchedToPFCandidate",&trk_matched);
  t1->SetBranchAddress("NVtx",&nvtx);
  t1->SetBranchAddress("NumInteractions",&numInteractions);
//#pragma omp parallel for 

  int nentries=0;// = 10000;
  if(batch==0) {
    nentries=10000;
    fprintf(OutFile,"event algo R jetid pt eta phi multiplicity girth mass trkpt suep isr suep_tot isr_tot nvtx numInteractions\n");
  }
  else{nentries=t1->GetEntries();}
  //printf("b %d e %d\n",batch,nentries);
  for(int entry=10000*batch; entry<nentries; entry++){
  //for(int entry=24; entry<25; entry++){
    if(entry%1000==0) {printf("entry %d/%d\n",entry,nentries);}
    //if(find(skip.begin(),skip.end(),entry) != skip.end()){printf("skipping entry: %d\n",entry); continue;}
    t1->GetEntry(entry);
    if(ht <1000){ continue;}
    vector<ROOT::Math::PtEtaPhiEVector> scalars;
    vector<ROOT::Math::PtEtaPhiEVector> isrs;

    // sort gen particles into isr and sueps
    for(int igen=0;igen<genPart->size();igen++){
      if(abs(gen_charge->at(igen)) != 1 || gen_status->at(igen) != 1 || genPart->at(igen).eta() >= 2.5 || genPart->at(igen).pt() <= 1){continue;}
      int pdgid = abs(gen_pdgId->at(igen));
      if((pdgid != 11) && (pdgid != 13) && (pdgid != 22) && (pdgid < 100 )){continue;}
      if(gen_parentId->at(igen) == 999998){ scalars.push_back(genPart->at(igen));} 
      else{ isrs.push_back(genPart->at(igen));} 

    }
    vector<PseudoJet> particles;
    for(int itrk=0;itrk<trkPart->size();itrk++){
      if ((trk_pv->at(itrk) <2) || (trk_matched->at(itrk) ==0)){continue;}

    double trkx = trkPart->at(itrk).x();
    double trky = trkPart->at(itrk).y();
    double trkz = trkPart->at(itrk).z();
    //ROOT::Math::PxPyPzEVector trk;
    TLorentzVector trk;
    //trk.SetPxPyPzE(trkx,trky,trkz,trkx*trkx+trky*trky+trkz*trkz+0.13957*0.13957);
    trk.SetXYZM(trkx,trky,trkz,0.13957);
    if(trk.Pt() <= 1 || trk.Eta() >= 2.5){continue;}
    particles.push_back(PseudoJet(trkx,trky,trkz,trkx*trkx+trky*trky+trkz*trkz+0.13957*0.13957));
    //printf("x: %f, y: %f, z: %f, E: %f, pt: %f, eta: %f, phi: %f\n",trk.X(),trk.Y(),trk.Z(),trk.E(),trk.Pt(),trk.Eta(),trk.Phi());
    //printf("x: %f, y: %f, z: %f, E: %f, pt: %f, eta: %f, phi: %f\n",trk.x(),trk.y(),trk.z(),trk.E(),trk.Pt(),trk.Eta(),trk.Phi());
    }

    ////run different clustering algorithms
    int algos[] = {-1,0,1};
    float Rs[] = {0.8,1.0,1.5,2.0};
//#pragma omp parallel for
    for(int algo : algos){
//#pragma omp parallel for
      for(float Ri : Rs){
    
        JetDefinition jet_def = getJetDef(algo,Ri);
        ClusterSequence cs(particles,jet_def);
        vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets(30));
        //vector<PseudoJet> jets = sorted_by_pt(cs.exclusive_jets(5));
        for(unsigned i = 0; i < jets.size(); i++) {
          PseudoJet jet = jets[i];
          auto [suep_part, isr_part] = jet_constituents(jet,scalars,isrs);
          auto [girth, trackpt] = jet_width(jet);
          //printf("%d %d %.1f %d %f %f %f %d %f %f %f %d %d %d %d %d %d\n",entry,algo, Ri,i,jet.pt(),jet.eta(),jet.phi(),jet.constituents().size(),girth,sqrt(jet.m()),trackpt, suep_part, isr_part,scalars.size(), isrs.size(), nvtx, numInteractions);   
          fprintf(OutFile,"%d %d %.1f %d %f %f %f %d %f %f %f %d %d %d %d %d %d\n",entry,algo, Ri,i,jet.pt(),jet.eta(),jet.phi(),jet.constituents().size(),girth,sqrt(jet.m()),trackpt, suep_part, isr_part,scalars.size(), isrs.size(), nvtx, numInteractions); 
        }
      }
    }

  }
fclose(OutFile);
}
  
