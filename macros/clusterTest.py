import numpy as np
from pyjet import cluster, DTYPE_PTEPM
import uproot_methods




print("starting python cluster test")
trkPartx = [9.511365,8.380543,8.128602,7.568648,-2.366235,-2.107425,-3.827330,-0.671908,-0.689713,-0.589932]
trkParty = [-22.425724,-20.914677,-18.222322,-17.835532, -4.444517,-4.525560,-3.168599,-0.765502,-0.744183,-0.824324]
trkPartz = [57.456248, 53.728303,47.222619,45.436706,1.189091,1.081467,0.806888, 0.412438,-0.514610,0.437422]
trke = []
particles = [];
for itrk, trkx in enumerate(trkPartx):
    trky = trkParty[itrk];
    trkz = trkPartz[itrk];
    trke.append(np.sqrt(trkx**2+trky**2+trkz**2+0.13957**2))
trk = uproot_methods.TLorentzVectorArray.from_cartesian(trkPartx,trkParty,trkPartz,trke);
#print(trk)

algos = [-1,0,1]
Rs = [0.8,1.0,1.5,2.0]
for algo in algos:
  for R in Rs:
    event = np.zeros(len(trk),dtype=DTYPE_PTEPM)
    event['pT'] = [p.pt for p in trk]
    event['eta'] = [p.eta for p in trk]
    event['phi'] = [p.phi for p in trk]
    #event['mass'] = [p.mass for p in trk]
    sequence = cluster(event,R=R,p=algo)
    jets = sequence.inclusive_jets(ptmin=30)
    
    for i,jet in enumerate(jets):
      print("%d %.1f %d %f %f %f %d %f"%(algo, R,i,jet.pt,jet.eta,jet.phi,len(jet.constituents()),jet.mass))
      print(jet.constituents())
