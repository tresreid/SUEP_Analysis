//#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TBranch.h"
#include "Math/PtEtaPhiE4D.h"
#include "Math/LorentzVector.h"
#include "Math/Vector4Dfwd.h"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"
using namespace fastjet;

#include <vector>

void jetAlgoComp(){
  printf("Starting\n");
  TFile *infile = TFile::Open("root://cmseos.fnal.gov//store/user/kdipetri/SUEP/Production_v0.2/2018/NTUP/PrivateSamples.SUEP_2018_mMed-1000_mDark-2_temp-2_decay-darkPhoHad_13TeV-pythia8_n-100_0_RA2AnalysisTree.root");

  TTree *t1= (TTree*)infile->Get("TreeMaker2/PreSelection");
  t1->Print();
  vector<ROOT::Math::PtEtaPhiEVector>* genPart = 0;
  //vector<ROOT::Math::Cartesian3D<double>>* trkPart = 0;
  vector<ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>>>* trkPart =0;
  double ht;
  vector<int>* gen_parentId;
  vector<int>* gen_charge;
  vector<int>* gen_pdgId;
  vector<int>* gen_status;
  vector<int>* trk_pv;
  vector<bool>* trk_matched;
  int nvtx;
  int numInteractions;
  t1->SetBranchAddress("GenParticles",&genPart);
  t1->SetBranchAddress("Tracks",&trkPart);
  t1->SetBranchAddress("HT",&ht);
  t1->SetBranchAddress("GenParticles_ParentId",&gen_parentId);
  t1->SetBranchAddress("GenParticles_Charge",&gen_charge);
  //t1->SetBranchAddress("GenParticles_PdgId",&gen_pdgId);
  //t1->SetBranchAddress("GenParticles_Status",&gen_status);
  //t1->SetBranchAddress("Tracks_fromPV0",&trk_pv);
  //t1->SetBranchAddress("Tracks_matchedToPFCandidate",&trk_matched);
  //t1->SetBranchAddress("Tracks_quality",&trk_quality);
  t1->SetBranchAddress("NVtx",&nvtx);
  t1->SetBranchAddress("NumInteractions",&numInteractions);
  
  for(Int_t entry=0; entry<50;entry++){
  //Int_t entry = 24;
    t1->GetEntry(entry);
    if(ht <1000){ continue;}
    vector<PseudoJet> particles;
    for(int itrk=0;itrk<10/*trkPart->size()*/;itrk++){
      //if (trk_pv->at(itrk) <2 || trk_matched->at(itrk) <=0){continue;}
      //ROOT::Math::LorentzVector trk;
      //printf("test %e\n", trkPart->at(100).x());
      //ROOT::Math::PtEtaPhiEVector trk(trkPart->at(itrk).x(),trkPart->at(itrk).y(),trkPart->at(itrk).z(),trkPart->at(itrk).x()*trkPart->at(itrk).x()+trkPart->at(itrk).y()*trkPart->at(itrk).y()+ trkPart->at(itrk).z()*trkPart->at(itrk).z()+0.13957*0.13957);
      //printf("evt: %d, trk %d: , pt:%f\n",entry,itrk,trk.Pt());
    //printf("event: %d; ht: %f; pt: %e; eta: %e; phi: %e\n",entry,ht,genPart->at(100).Pt(),genPart->at(100).Eta(),genPart->at(100).Phi());
    //printf("event: %d; ht: %f; pt: %e; eta: %e; phi: %e\n",entry,ht,trkPart->at(100).x(),trkPart->at(100).y(),trkPart->at(100).z());
    printf("event: %d; ht: %f; pt: %e; eta: %e; phi: %e\n",entry,ht,trkPart->at(itrk).x(),trkPart->at(itrk).y(),trkPart->at(itrk).z());

    double trkx = trkPart->at(itrk).x();
    double trky = trkPart->at(itrk).y();
    double trkz = trkPart->at(itrk).z();
    ROOT::Math::PxPyPzEVector trk;//(0,0,0,0);
    trk.SetPxPyPzE(trkx,trky,trkz,trkx*trkx+trky*trky+trkz*trkz+0.13957*0.13957);
    particles.push_back(PseudoJet(trkx,trky,trkz,trkx*trkx+trky*trky+trkz*trkz+0.13957*0.13957));
    
    //printf("event: %d, trk %d, pt: %f eta %f phi %f\n",entry,itrk,trk.Pt(),trk.Eta(),trk.Phi());
    }
    //printf("event: %d; ht: %f; pt: %e; eta: %e; phi: %e\n",entry,ht,genPart->at(100).Pt(),genPart->at(100).Eta(),genPart->at(100).Phi());
    //printf("event: %d; ht: %f; pt: %e; eta: %e; phi: %e\n",entry,ht,trkPart->at(100).z(),trkPart->at(100).y(),trkPart->at(100).z());
  }
}
  
