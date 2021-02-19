import uproot
import uproot_methods
#import ROOT
#import coffea
import math
#from coffea import hist
import numpy as np
#import seutils
from math import pi
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import colors
from pyjet import cluster, DTYPE_PTEPM
from pyjet.testdata import get_event
from pathlib import Path
import sys
#from joblib import Parallel,delayed

#ROOT.gSystem.Load('libGenVector')

eta_min, eta_max = -4., 4.
extent = [-np.pi, np.pi, eta_min, eta_max]
bins = 200
eta_edges = np.linspace(eta_min, eta_max, bins + 1)
phi_edges = np.linspace(-np.pi, np.pi, bins + 1)

def jet_width(jet):
    girth = 0
    pt_avg =0
    count =0
    for i,constituent in enumerate([x for x in jet.constituents_array() if x['pT'] > 1]):
      if constituent['pT'] < 1:
        continue
      #min_dR = 9999
      #is_scalar = True
      #get scalar constituents
      phi = constituent['phi']
      eta = constituent['eta']
      pt = constituent['pT']
      dPhi = abs(jet.phi-phi) 
      if dPhi > math.pi:
        dPhi = dPhi - 2*math.pi
      dEta = jet.eta-eta
      dR = math.sqrt(dEta*dEta + dPhi*dPhi)
      girth += pt * dR 
      pt_avg += pt
      count += 1
      #if dR < min_dR:
      #  is_scalar = True
      #  min_dR = dR
    return ((girth / jet.pt),((pt_avg/count)/jet.pt))

def jet_constituents(jet,scalars,isrs):
    scalar_part = 0
    isr_part = 0
    for i,constituent in enumerate([x for x in jet.constituents_array() if x['pT'] > 1]):
      if constituent['pT'] < 1:
        continue
      min_dR = 9999
      is_scalar = True
      #get scalar constituents
      phi = constituent['phi']
      eta = constituent['eta']
      for scalar in scalars:
        #if abs(scalar.pt -constituent['pT'])/scalar.pt > 0.3:
        #  continue
        dPhi = abs(scalar.phi-phi) 
        if dPhi > math.pi:
          dPhi = dPhi - 2*math.pi
        dEta = scalar.eta-eta
        dR = math.sqrt(dEta*dEta + dPhi*dPhi)
        if dR < min_dR:
          is_scalar = True
          min_dR = dR
          
      #get isr constituents
      for isr in isrs:
        #if abs(scalar.pt -constituent['pT'])/scalar.pt > 0.3:
        #  continue
        dPhi = abs(isr.phi-phi) 
        if dPhi > math.pi:
          dPhi = dPhi - 2*math.pi
        dEta = isr.eta-eta
        dR = math.sqrt(dEta*dEta + dPhi*dPhi)
        if dR < min_dR:
          is_scalar = False
          min_dR = dR
      if is_scalar:
        scalar_part = scalar_part +1
      else:
        isr_part = isr_part +1
    return (scalar_part, isr_part)

def makeJets(jet_algo,R, particles):
  print(particles)
  cone = R*10
  eta = np.linspace(eta_min, eta_max, bins + 1)[:-1] + (eta_max - eta_min) / (2 * bins)
  phi = np.linspace(-np.pi, np.pi, bins + 1)[:-1] + (np.pi / bins)
  X, Y = np.meshgrid(phi, eta)
  event = np.zeros(len(particles),dtype=DTYPE_PTEPM)
  event['pT'] = [p.pt for p in particles]
  event['eta'] = [p.eta for p in particles]
  event['phi'] = [p.phi for p in particles]
  event['mass'] = [p.mass for p in particles]
  #ghosts = np.zeros(eta.shape[0] * phi.shape[0], dtype=DTYPE_PTEPM)
  #ghosts['pT'] = 1e-8
  #ghosts['eta'] = Y.ravel()
  #ghosts['phi'] = X.ravel()
  
  # add ghosts to the event
  #event = np.concatenate([event, ghosts], axis=0)

  # p = -1 (ak), 0 (CA), 1 (kt)
  sequence = cluster(event,R=R,p=jet_algo)
  jets = sequence.inclusive_jets(ptmin=30)
  #jets = sequence.exclusive_jets(5)
  return jets

# Get the file and import using uproot
mMed = 1000
mDark = 2
temp = 2
#decayMode = 'generic'
decayMode = 'darkPhoHad'
#base = '/Users/chrispap/'
# xrootd is not working properly in Python3 :(
base = 'root://cmseos.fnal.gov//store/user/kdipetri/SUEP/Production_v0.2/2018/NTUP/'
#datasets = [base +
#            'PrivateSamples.SUEP_2018_mMed-%d_mDark-%d_temp-%d_decay-%s'
#            '_13TeV-pythia8_n-100_0_RA2AnalysisTree.root'%(mMed, mDark, temp, decayMode),
#           ]
#datasets=[]
#datasets.append("output/mMed-125_mDark-2_temp-2_decay-generic.root")
#datasets.append("output/mMed-400_mDark-2_temp-2_decay-generic.root")
#datasets.append("output/mMed-750_mDark-2_temp-2_decay-generic.root")
#datasets.append("output/mMed-1000_mDark-2_temp-2_decay-generic.root")
datasets=[
base + 'PrivateSamples.SUEP_2018_mMed-%d_mDark-%d_temp-%d_decay-%s_13TeV-pythia8_n-100_0_RA2AnalysisTree.root'%(1000, mDark, temp, decayMode),
base + 'PrivateSamples.SUEP_2018_mMed-%d_mDark-%d_temp-%d_decay-%s_13TeV-pythia8_n-100_0_RA2AnalysisTree.root'%(400, mDark, temp, decayMode),
base + 'Autumn18.QCD_HT300to500_TuneCP5_13TeV-madgraphMLM-pythia8_RA2AnalysisTree.root',
base + 'Autumn18.QCD_HT500to700_TuneCP5_13TeV-madgraphMLM-pythia8_RA2AnalysisTree.root',
base + 'Autumn18.QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8_RA2AnalysisTree.root',
base + 'Autumn18.QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8_RA2AnalysisTree.root',
base + 'Autumn18.QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8_RA2AnalysisTree.root',
base + 'Autumn18.QCD_HT2000toInf_TuneCP5_13TeV-madgraphMLM-pythia8_RA2AnalysisTree.root'
]
if len(sys.argv) == 2:
  ntup = int(sys.argv[1])
  if ntup >= len(datasets):
    print("input greater than max dataset, setting to 0")
    ntup = 0
else:
  ntup = 0
rootfile = datasets[ntup]
#fin = ROOT.TFile.Open(rootfile)
print("opening rootfile: %s"%rootfile)
fin = uproot.open(rootfile)
print("rootfiles opened: ")

# Attach the branches to numpy arrays
tree = fin[b'TreeMaker2/PreSelection']
def get_branch(branchname):
    return tree[branchname].array()
if("SUEP" in datasets[ntup]):
  GenParticles_pt = get_branch('GenParticles.fCoordinates.fPt')
  GenParticles_eta = get_branch('GenParticles.fCoordinates.fEta')
  GenParticles_phi = get_branch('GenParticles.fCoordinates.fPhi')
  GenParticles_E = get_branch('GenParticles.fCoordinates.fE')
  GenParticles_ParentId = get_branch('GenParticles_ParentId')
  GenParticles_Charge = get_branch('GenParticles_Charge')
  GenParticles_PdgId = get_branch('GenParticles_PdgId')
  GenParticles_Status = get_branch('GenParticles_Status')
#extra jet info
print("getting branches")
Jets = get_branch('Jets')
HT = get_branch('HT')
#Tracks = get_branch("Tracks_charge")#uproot.interpret(tree["Tracks"],tobject=False))
#tree.show()
Tracks_x = get_branch("Tracks.fCoordinates.fX")
Tracks_y = get_branch("Tracks.fCoordinates.fY")
Tracks_z = get_branch("Tracks.fCoordinates.fZ")
Tracks_fromPV0 = get_branch("Tracks_fromPV0")
Tracks_matched = get_branch("Tracks_matchedToPFCandidate")
Tracks_quality = get_branch("Tracks_quality")
#Tracks = tree.arrays["Tracks.Eta()"]

NVtx = get_branch("NVtx")
NumInteractions = get_branch("NumInteractions")

print("branches got")

def pass_HT(ievt):
  ht = HT[ievt]
  print("HT: %f"%(ht))
  if ht < 1000:
    return False
  else:
    return True

# Main plotting function
def plot(ievt,outfile,outfile_all, ax=None, boost=False):#, reco=False, jet_algo=-1, R=0.8):
    # Get the particles of ievt event
    nvtx = NVtx[ievt]
    numInteractions = NumInteractions[ievt]
    reco_jets = Jets[ievt]
    tracks_x = Tracks_x[ievt]
    tracks_y = Tracks_y[ievt]
    tracks_z = Tracks_z[ievt]
    tracks_fromPV0 = Tracks_fromPV0[ievt]
    tracks_matched = Tracks_matched[ievt]
    tracks_quality = Tracks_quality[ievt]
    tracks_E = np.sqrt(tracks_x**2+tracks_y**2+tracks_z**2+0.13957**2) # 0.13957 is the pi mass. 1st approximation. Not all tracks come from pions. hard to tell between pions and other objects at low pt?
    tracks = uproot_methods.TLorentzVectorArray.from_cartesian(tracks_x,
                                                               tracks_y,
                                                               tracks_z,
                                                               tracks_E)

    tracks = tracks[(tracks.pt >1.) & (tracks.eta <2.5) & (tracks_fromPV0 >=2) & (tracks_matched >0) ]#& (tracks_quality <2)]

    #if("signal" in ofile_name):
    if True:
      genParticles_pt = GenParticles_pt[ievt]
      genParticles_eta = GenParticles_eta[ievt]
      genParticles_phi = GenParticles_phi[ievt]
      genParticles_E = GenParticles_E[ievt]
      genParticles = uproot_methods.TLorentzVectorArray.from_ptetaphim(genParticles_pt,
                                                                       genParticles_eta,
                                                                       genParticles_phi,
                                                                       genParticles_E)
      genParticles_ParentId = GenParticles_ParentId[ievt]
      genParticles_PdgId = GenParticles_PdgId[ievt]
      genParticles_Status = GenParticles_Status[ievt]
      genParticles_Charge = GenParticles_Charge[ievt]
       #The last copy of the scalar mediator
      scalarParticle = genParticles[(genParticles_PdgId == 25) & (genParticles_Status == 62)]

      # Define mask arrays to select the desired particles
      finalParticles = (genParticles_Status == 1) & (genParticles.pt > 1) & (abs(genParticles_Charge) ==1) & (abs(genParticles.eta) < 2.5)
      fromScalarParticles = genParticles_ParentId == 999998
      isrParticles = genParticles_ParentId != 999998

      # Apply the selection criteria to get the final particle arrays
      # 10 arrays of final particles in total
      # Dividing to e, mu, gamma, pi, all other hadrons
      # for particles that come from the scalar mediator or not
      fromScalarParticles_e = genParticles[finalParticles &
                                           fromScalarParticles &
                                           (abs(genParticles_PdgId) == 11)]
      fromScalarParticles_mu = genParticles[finalParticles &
                                            fromScalarParticles &
                                            (abs(genParticles_PdgId) == 13)]
      fromScalarParticles_gamma = genParticles[finalParticles &
                                               fromScalarParticles &
                                               (abs(genParticles_PdgId) == 22)]
      fromScalarParticles_pi = genParticles[finalParticles &
                                            fromScalarParticles &
                                            (abs(genParticles_PdgId) == 211)]
      fromScalarParticles_hadron = genParticles[finalParticles &
                                                fromScalarParticles &
                                                (abs(genParticles_PdgId) > 100) & (abs(genParticles_PdgId) != 211)]

      isrParticles_e = genParticles[finalParticles & isrParticles &
                                    (abs(genParticles_PdgId) == 11)]
      isrParticles_mu = genParticles[finalParticles & isrParticles &
                                     (abs(genParticles_PdgId) == 13)]
      isrParticles_gamma = genParticles[finalParticles & isrParticles &
                                        (abs(genParticles_PdgId) == 22)]
      isrParticles_pi = genParticles[finalParticles & isrParticles &
                                     (abs(genParticles_PdgId) == 211)]
      isrParticles_hadron = genParticles[finalParticles & isrParticles &
                                         (abs(genParticles_PdgId) > 100)& (abs(genParticles_PdgId) != 211)]

      particles = np.concatenate([fromScalarParticles_e,fromScalarParticles_mu,fromScalarParticles_gamma,fromScalarParticles_pi,fromScalarParticles_hadron ,isrParticles_e , isrParticles_mu , isrParticles_gamma , isrParticles_pi, isrParticles_hadron], axis=0)  
      isrs = np.concatenate([isrParticles_e , isrParticles_mu , isrParticles_gamma , isrParticles_pi, isrParticles_hadron], axis=0)  
      scalars = np.concatenate([fromScalarParticles_e ,fromScalarParticles_mu ,fromScalarParticles_gamma ,fromScalarParticles_pi ,fromScalarParticles_hadron], axis=0)
      #  truthScalar = truthScalar + scalar 
  

    # Function that sets the scaling for markers
    # Two methods for the moment: use tanh or square.
    # Scaling using just the energy is also an option
    def scale(particles, scalar, method=3):
        """Just to scale to a reasonable dot size"""
        #if len(particles) == 0: return []
        energies = particles.energy
        e_max = scalar.energy
        if len(energies) == 0: return []
        if method == 1:
            e_normed = 500000.*np.square(energies/e_max)
        elif method == 2:
            e_normed = 1000.*np.tanh(energies/e_max)
        else:
            e_normed = 2500.*energies/e_max
        return e_normed
    
    for ja in [-1,0,1]:
      if ja == -1:
        jet_algo_title = "AKT"
      elif ja == 0:
        jet_algo_title = "CA"
      elif ja == 1:
        jet_algo_title = "KT"
      for Ri in [0.8,1,1.5,2]:
        lead_jet = [None]*2
        sublead_jet = [None]*2
        multi_jet = [None]*2
        submulti_jet = [None]*2
        #for rec in [1]:
        #  if rec:
        parts = tracks
        reco_title = "Reco"
        #else:
        #  parts = particles
        #  reco_title = "Gen"

        # Initialize plotting
        #if ax is None:
        fig = plt.figure(figsize=(8,8))
        ax = plt.gca()

        # Plot parameters
        ax.set_xlim(-np.pi, np.pi)
        ax.set_ylim(eta_min, eta_max)
        ax.set_xlabel(r'$\phi$', fontsize=18)
        ax.set_ylabel(r'$\eta$', fontsize=18)
        ax.tick_params(axis='both', which='major', labelsize=12)
        # Add scatters to figure
        ax.scatter(tracks.phi, tracks.eta,
                   #s=scale(fromScalarParticles_e,scalarParticle),
                   c='xkcd:blue', marker='o')
        #for trk in parts:
        #  print("x: %f, y: %f, z: %f, E: %f, pt: %f, eta: %f, phi: %f"%(trk.x,trk.y,trk.z,trk.E,trk.pt,trk.eta, trk.phi))
        new_jets = makeJets(ja,Ri,parts)
        for trk in new_jets:
          print("pt: %f, eta: %f, phi: %f"%(trk.pt,trk.eta, trk.phi))
        
        area = np.zeros((phi_edges.shape[0] - 1, eta_edges.shape[0] - 1),
                        dtype=np.float64)
        leading_jetid = -1
        subleading_jetid = -1
        leading_jetpt = -1
        subleading_jetpt = -1
        multi_jetid = -1
        submulti_jetid = -1
        multi_jetntrk = -1
        submulti_jetntrk = -1
        
        for ijet, jet in enumerate(new_jets):
          if jet.pt > leading_jetpt:
            subleading_jetid = leading_jetid
            subleading_jetpt = leading_jetpt
            leading_jetid = ijet
            leading_jetpt = jet.pt
          elif jet.pt >= subleading_jetpt:
            subleading_jetid = ijet
            subleading_jetpt = jet.pt
          if len(jet.constituents_array()) > multi_jetntrk:
            submulti_jetid = multi_jetid
            submulti_jetntrk = multi_jetntrk
            multi_jetid = ijet
            multi_jetntrk = len(jet.constituents_array())
          elif len(jet.constituents_array()) >= submulti_jetntrk:
            submulti_jetid = ijet
            submulti_jetntrk = len(jet.constituents_array())
        for ijet, jet in enumerate(new_jets):
          constit = jet.constituents_array()
          #print(constit)
          jetarea, _, _ = np.histogram2d(constit['phi'],constit['eta'],bins=(phi_edges,eta_edges))
          if ijet== leading_jetid:
            area += (jetarea > 0) *(4)
          elif ijet== subleading_jetid:
            area += (jetarea > 0) *(3)
          elif jet.pt > 100:
            area += (jetarea > 0) *(2)
          else:
            area += (jetarea > 0) *(1)
              
            
        
        
        cmap = colors.ListedColormap(['blue','cyan','green','red'])
        bounds = [0.5,1.5,2.5,3.5,4.5]
        norm = colors.BoundaryNorm(bounds,cmap.N)
        ax.imshow(np.ma.masked_where(area == 0, area).T, cmap=cmap,
                  extent=extent,aspect='auto',norm=norm,# aspect=(eta_max - eta_min) / (2*np.pi),
                  interpolation='none', origin='lower')
        # Create custom legends
        # Legend 1 is particle type
        line7 = ax.scatter([-100], [-100], label='leading jet (pt:.2%f)'%(leading_jetpt), marker='o',c='xkcd:red')
        line8 = ax.scatter([-100], [-100], label='subleading jet (pt:.2%f)'%(subleading_jetpt), marker='o',c='xkcd:green')
        line9 = ax.scatter([-100], [-100], label='jet (pt>100)', marker='o',c='xkcd:cyan')
        line10 = ax.scatter([-100], [-100], label='jet (30< pt < 100)', marker='o',c='xkcd:blue')
        ax.legend(handles=[line7, line8,line9,line10],
                                  loc='upper right', fontsize=12)

        # Legend 2 is about particle origin
        #blue_patch = mpatches.Patch(color='xkcd:blue', label='from scalar')
        #magenta_patch = mpatches.Patch(color='xkcd:magenta', label='not from scalar')
        #plt.legend(handles=[blue_patch, magenta_patch],loc='upper left')

        # build a rectangle in axes coords
        left, width = .0, 1.
        bottom, height = .0, 1.
        center = left + width/2.
        right = left + width
        top = bottom + height

        # axes coordinates are 0,0 is bottom left and 1,1 is upper right
        p = mpatches.Rectangle((left, bottom), width, height,
            fill=False, transform=ax.transAxes, clip_on=False)

        ax.add_patch(p)
        # Print event number
        ax.text(left, top, 'Event %d'%ievt, horizontalalignment='left',
                verticalalignment='bottom', transform=ax.transAxes, fontsize=12)
        # Print sample details
        ax.text(right, top, 'QCD',
                horizontalalignment='right', verticalalignment='bottom',
                transform=ax.transAxes, fontsize=12)
        Path("Results/QCD/%s/%s_%d/"%(reco_title,jet_algo_title,int(Ri*10))).mkdir(parents=True,exist_ok=True)
        fig.savefig('Results/QCD/%s/%s_%d/QCD_Event%d.png'%(reco_title,jet_algo_title,int(Ri*10),ievt))
        fig.clear()
        plt.close(fig)

        if leading_jetid < len(new_jets) and multi_jetid < len(new_jets) and len(lead_jet) > 0 and leading_jetid >=0: 
          lead_jet[0] = new_jets[leading_jetid]
          multi_jet[0] = new_jets[multi_jetid]
          (girth_l, trackpt_l) = jet_width(lead_jet[0])
          (girth_m, trackpt_m) = jet_width(multi_jet[0])
        
          #if("signal" in ofile_name):
          if True:
            (lead_scalars_l, lead_isr_l) = jet_constituents(lead_jet[0], scalars, isrs)
            (lead_scalars_m, lead_isr_m) = jet_constituents(multi_jet[0], scalars, isrs)
            (suep_tracks,isr_tracks) = (len(scalars),len(isrs))
          else:
            (lead_scalars_l, lead_isr_l) = (0,0) 
            (lead_scalars_m, lead_isr_m) = (0,0) 
            (suep_tracks,isr_tracks) = (0,0)
          outfile.write("%d %s %.1f %f %f %d %d %f %f %f %f %f %f %d %d %d %d %d %d\n"%(ievt,jet_algo_title,Ri,
            lead_jet[0].pt, multi_jet[0].pt,
            len([x for x in lead_jet[0].constituents_array() if x['pT'] > 1]),len([x for x in multi_jet[0].constituents_array() if x['pT'] > 1]),
            girth_l, girth_m,
            lead_jet[0].mass ,multi_jet[0].mass,
            trackpt_l, trackpt_m,
            lead_scalars_l, lead_isr_l,lead_scalars_m, lead_isr_m,suep_tracks,isr_tracks
          ))
        #save all jets, not just leading
        for jet_i in range(len(new_jets)):
          (girth_all, trackpt_all) = jet_width(new_jets[jet_i])
        
          #if("signal" in ofile_name):
          if True:
            (lead_scalars_all, lead_isr_all) = jet_constituents(new_jets[jet_i], scalars, isrs)
            (suep_tracks1,isr_tracks1) = (len(scalars),len(isrs))
          else:
            (lead_scalars_all, lead_isr_all) = (0,0) 
            (suep_tracks,isr_tracks) = (0,0)
          outfile_all.write("%d %s %.1f %d %f %d %f %f %f %d %d %d %d %d %d\n"%(ievt,jet_algo_title,Ri, jet_i,
            new_jets[jet_i].pt,
            len([x for x in new_jets[jet_i].constituents_array() if x['pT'] > 1]),
            girth_all,
            new_jets[jet_i].mass,
            trackpt_all,
            lead_scalars_all, lead_isr_all,suep_tracks,isr_tracks,
            nvtx,numInteractions
          ))
         
# The program runs through this loop
# if you want multiple plots then change multi to True
# Warning: multiplot funtion is currently broken
multi = False
if multi:
    figure, axs = plt.subplots(nrows=2, ncols=3, figsize=(16,8))
    axs_flat = axs.flatten()

boost = False
event = 0
ofile_names = ["signal1000","signal400","qcd300","qcd500","qcd700","qcd1000","qcd1500","qcd2000"]
#outfile = open("jet_clustering/%s_comparisonsv2.txt"%ofile_names[ntup],"w")
#outfile.write("Event Jet_algo R pt_l pt_m ntracks_l ntracks_m girth_l girth_m mass_l mass_m trackpt_l trackpt_m suep_tracks_l isr_tracks_l suep_tracks_m isr_tracks_m total_suep total_isr\n") 
#outfile_all = open("jet_clustering/%s_comparisons_allv2.txt"%ofile_names[ntup],"w")

outfile = open("test.txt","w")
outfile_all = open("testall.txt","w")

outfile_all.write("Event Jet_algo R jet_id pt_l ntracks_l girth_l mass_l trackpt_l suep_tracks_l isr_tracks_l total_suep total_isr NVtx NumInteractions\n") 

#def run(i):
for i in [24]:#,45,47]:
#for i in range(0,len(tree)):
  print("Event %d/%d"%(i,len(tree)))
  ht_pass = pass_HT(i)
  if ht_pass:
    plot(i,outfile,outfile_all,boost=boost)#,reco=rec,jet_algo=ja,R=Ri)

#Parallel(n_jobs=4)(delayed(run)(i) for i in [24,45,47])
outfile.close()
outfile_all.close()
#plt.show();
