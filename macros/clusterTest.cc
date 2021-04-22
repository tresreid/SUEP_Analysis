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

  vector<float> trkPartx = {9.511365,8.380543,8.128602,7.568648,-2.366235,-2.107425,-3.827330,-0.671908,-0.689713,-0.589932};
  vector<float> trkParty = {-22.425724,-20.914677,-18.222322,-17.835532, -4.444517,-4.525560,-3.168599,-0.765502,-0.744183,-0.824324};
  vector<float> trkPartz = {57.456248, 53.728303,47.222619,45.436706,1.189091,1.081467,0.806888, 0.412438,-0.514610,0.437422};
  vector<PseudoJet> particles;
  for(int itrk; itrk<trkPartx.size(); itrk++){
    double trkx = trkPartx.at(itrk);
    double trky = trkParty.at(itrk);
    double trkz = trkPartz.at(itrk);
    TLorentzVector trk;
    trk.SetXYZM(trkx,trky,trkz,0.13957);
    particles.push_back(PseudoJet(trkx,trky,trkz,trkx*trkx+trky*trky+trkz*trkz+0.13957*0.13957));
    //printf("x: %f, y: %f, z: %f, E: %f, pt: %f, eta: %f, phi: %f\n",trk.X(),trk.Y(),trk.Z(),trk.E(),trk.Pt(),trk.Eta(),trk.Phi());
    //printf("x: %f, y: %f, z: %f, E: %f, pt: %f, eta: %f, phi: %f\n",trk.x(),trk.y(),trk.z(),trk.E(),trk.Pt(),trk.Eta(),trk.Phi());
  }

    ////run different clustering algorithms
    int algos[] = {-1,0,1};
    float Rs[] = {0.8,1.0,1.5,2.0};
    for(int algo : algos){
      for(float Ri : Rs){
    
        JetDefinition jet_def = getJetDef(algo,Ri);
        ClusterSequence cs(particles,jet_def);
        vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets(30));
        for(unsigned i = 0; i < jets.size(); i++) {
          PseudoJet jet = jets[i];
          //printf("%d %.1f %d %f %f %f %d %f\n",algo, Ri,i,jet.pt(),jet.eta(),jet.phi()>M_PI? jet.phi()-2*M_PI:jet.phi(),jet.constituents().size(),sqrt(jet.m()));
          for (unsigned j =0; j< jet.constituents().size(); j++){
            printf("%d %d %f %f %f %f\n",i,j,jet.constituents().at(j).pt(),jet.constituents().at(j).eta(),jet.constituents().at(j).phi(),jet.constituents().at(j).m());
          }
        }
      }
    }

}
  
