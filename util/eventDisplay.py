import uproot
import uproot_methods
import ROOT
import coffea
from coffea import hist
import numpy as np
import seutils
from math import pi
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import colors
from pyjet import cluster, DTYPE_PTEPM
from pyjet.testdata import get_event

ROOT.gSystem.Load('libGenVector')

eta_min, eta_max = -4., 4.
extent = [-np.pi, np.pi, eta_min, eta_max]
bins = 200
eta_edges = np.linspace(eta_min, eta_max, bins + 1)
phi_edges = np.linspace(-np.pi, np.pi, bins + 1)

def add_jets(ax,jet_algo,R,particles):
    new_jets = makeJets(jet_algo,R,particles)
    area = np.zeros((phi_edges.shape[0] - 1, eta_edges.shape[0] - 1),
                    dtype=np.float64)
    leading_jetid = -1
    subleading_jetid = -1
    leading_jetpt = -1
    subleading_jetpt = -1
    
    for ijet, jet in enumerate(new_jets):
      if jet.pt > leading_jetpt:
        subleading_jetid = leading_jetid
        subleading_jetpt = leading_jetpt
        leading_jetid = ijet
        leading_jetpt = jet.pt
      elif jet.pt >= subleading_jetpt:
        subleading_jetid = ijet
        subleading_jetpt = jet.pt
    for ijet, jet in enumerate(new_jets):
      constit = jet.constituents_array()
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
    line1 = ax.scatter([-100], [-100], label='$e$',marker='o', c='xkcd:black')
    line2 = ax.scatter([-100], [-100], label='$\mu$', marker='v', c='xkcd:black')
    line3 = ax.scatter([-100], [-100], label='$\gamma$', marker='s', c='xkcd:black')
    line4 = ax.scatter([-100], [-100], label='$\pi$', marker='P', c='xkcd:black')
    line5 = ax.scatter([-100], [-100], label='other hadron', marker='*', c='xkcd:black')
    line6 = ax.scatter([-100], [-100], label='Scalar mediator', marker='o',
                      facecolors='none', edgecolors='r')
    line7 = ax.scatter([-100], [-100], label='leading jet (pt:%f)'%(leading_jetpt), marker='o',c='xkcd:red')
    line8 = ax.scatter([-100], [-100], label='subleading jet (pt:%f)'%(subleading_jetpt), marker='o',c='xkcd:green')
    line9 = ax.scatter([-100], [-100], label='jet (pt>100)', marker='o',c='xkcd:cyan')
    line10 = ax.scatter([-100], [-100], label='jet (30< pt < 100)', marker='o',c='xkcd:blue')
    #first_legend = plt.legend(handles=[line1, line2, line3, line4, line5, line6,line7, line8,line9,line10],
    ax.legend(handles=[line1, line2, line3, line4, line5, line6,line7, line8,line9,line10],
                              loc='upper right', fontsize=12)
    #ax.add_artist(first_legend)
    return (ax,leading_jetpt,subleading_jetpt)


def makeJets(jet_algo,R, particles):
  cone = R*10
  eta = np.linspace(eta_min, eta_max, bins + 1)[:-1] + (eta_max - eta_min) / (2 * bins)
  phi = np.linspace(-np.pi, np.pi, bins + 1)[:-1] + (np.pi / bins)
  X, Y = np.meshgrid(phi, eta)
  event = np.zeros(len(particles),dtype=DTYPE_PTEPM)
  event['pT'] = [p.pt for p in particles]
  event['eta'] = [p.eta for p in particles]
  event['phi'] = [p.phi for p in particles]
  ghosts = np.zeros(eta.shape[0] * phi.shape[0], dtype=DTYPE_PTEPM)
  ghosts['pT'] = 1e-8
  ghosts['eta'] = Y.ravel()
  ghosts['phi'] = X.ravel()
  
  # add ghosts to the event
  event = np.concatenate([event, ghosts], axis=0)

  # p = -1 (ak), 0 (CA), 1 (kt)
  sequence = cluster(event,R=R,p=jet_algo)
  jets = sequence.inclusive_jets(ptmin=30)
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
datasets = [base +
            'PrivateSamples.SUEP_2018_mMed-%d_mDark-%d_temp-%d_decay-%s'
            '_13TeV-pythia8_n-100_0_RA2AnalysisTree.root'%(mMed, mDark, temp, decayMode),
           ]
#datasets=[]
#datasets.append("output/mMed-125_mDark-2_temp-2_decay-generic.root")
#datasets.append("output/mMed-400_mDark-2_temp-2_decay-generic.root")
#datasets.append("output/mMed-750_mDark-2_temp-2_decay-generic.root")
#datasets.append("output/mMed-1000_mDark-2_temp-2_decay-generic.root")
rootfile = datasets[0]
#fin = ROOT.TFile.Open(rootfile)
fin = uproot.open(rootfile)

# Attach the branches to numpy arrays
tree = fin[b'TreeMaker2/PreSelection']
def get_branch(branchname):
    return tree[branchname].array()

GenParticles_pt = get_branch('GenParticles.fCoordinates.fPt')
GenParticles_eta = get_branch('GenParticles.fCoordinates.fEta')
GenParticles_phi = get_branch('GenParticles.fCoordinates.fPhi')
GenParticles_E = get_branch('GenParticles.fCoordinates.fE')
GenParticles_ParentId = get_branch('GenParticles_ParentId')
GenParticles_PdgId = get_branch('GenParticles_PdgId')
GenParticles_Status = get_branch('GenParticles_Status')
#extra jet info
Jets = get_branch('Jets')
HT = get_branch('HT')
#Tracks = get_branch("Tracks_charge")#uproot.interpret(tree["Tracks"],tobject=False))
#tree.show()
Tracks_x = get_branch("Tracks.fCoordinates.fX")
Tracks_y = get_branch("Tracks.fCoordinates.fY")
Tracks_z = get_branch("Tracks.fCoordinates.fZ")
Tracks_fromPV0 = get_branch("Tracks_fromPV0")
Tracks_matched = get_branch("Tracks_matchedToPFCandidate")
#Tracks = tree.arrays["Tracks.Eta()"]

def pass_HT(ievt):
  ht = HT[ievt]
  print("HT: %f"%(ht))
  if ht < 1000:
    return False
  else:
    return True

# Main plotting function
def plot(ievt, ax=None, boost=False):#, reco=False, jet_algo=-1, R=0.8):
    # Get the particles of ievt event
    genParticles_pt = GenParticles_pt[ievt]
    genParticles_eta = GenParticles_eta[ievt]
    genParticles_phi = GenParticles_phi[ievt]
    genParticles_E = GenParticles_E[ievt]
    genParticles = uproot_methods.TLorentzVectorArray.from_ptetaphie(genParticles_pt,
                                                                     genParticles_eta,
                                                                     genParticles_phi,
                                                                     genParticles_E)
    genParticles_ParentId = GenParticles_ParentId[ievt]
    genParticles_PdgId = GenParticles_PdgId[ievt]
    genParticles_Status = GenParticles_Status[ievt]
    reco_jets = Jets[ievt]
    tracks_x = Tracks_x[ievt]
    tracks_y = Tracks_y[ievt]
    tracks_z = Tracks_z[ievt]
    tracks_fromPV0 = Tracks_fromPV0[ievt]
    tracks_matched = Tracks_matched[ievt]
    tracks_E = np.sqrt(tracks_x**2+tracks_y**2+tracks_z**2+0.13957**2) # 0.13957 is the pi mass. 1st approximation. Not all tracks come from pions. hard to tell between pions and other objects at low pt?
    tracks = uproot_methods.TLorentzVectorArray.from_cartesian(tracks_x,
                                                               tracks_y,
                                                               tracks_z,
                                                               tracks_E)
    # The last copy of the scalar mediator
    scalarParticle = genParticles[(genParticles_PdgId == 25) & (genParticles_Status == 62)]

    # Define mask arrays to select the desired particles
    finalParticles = (genParticles_Status == 1) & (genParticles.pt > 1)
    fromScalarParticles = genParticles_ParentId == 999998
    isrParticles = genParticles_ParentId != 999998

    tracks = tracks[(tracks.pt >1.) & (tracks.eta <2.5) & (tracks_fromPV0 >=2) & (tracks_matched >0)]

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
                                              (abs(genParticles_PdgId) > 100)]

    isrParticles_e = genParticles[finalParticles & isrParticles &
                                  (abs(genParticles_PdgId) == 11)]
    isrParticles_mu = genParticles[finalParticles & isrParticles &
                                   (abs(genParticles_PdgId) == 13)]
    isrParticles_gamma = genParticles[finalParticles & isrParticles &
                                      (abs(genParticles_PdgId) == 22)]
    isrParticles_pi = genParticles[finalParticles & isrParticles &
                                   (abs(genParticles_PdgId) == 211)]
    isrParticles_hadron = genParticles[finalParticles & isrParticles &
                                       (abs(genParticles_PdgId) > 100)]

    # Boost everything to scalar's rest frame
    if boost == True:
        fromScalarParticles_e = fromScalarParticles_e.boost(-scalarParticle.p3/scalarParticle.energy)
        fromScalarParticles_mu = fromScalarParticles_mu.boost(-scalarParticle.p3/scalarParticle.energy)
        fromScalarParticles_gamma = fromScalarParticles_gamma.boost(-scalarParticle.p3/scalarParticle.energy)
        fromScalarParticles_pi = fromScalarParticles_pi.boost(-scalarParticle.p3/scalarParticle.energy)
        fromScalarParticles_hadron = fromScalarParticles_hadron.boost(-scalarParticle.p3/scalarParticle.energy)
        isrParticles_e = isrParticles_e.boost(-scalarParticle.p3/scalarParticle.energy)
        isrParticles_mu = isrParticles_mu.boost(-scalarParticle.p3/scalarParticle.energy)
        isrParticles_gamma = isrParticles_gamma.boost(-scalarParticle.p3/scalarParticle.energy)
        isrParticles_pi = isrParticles_pi.boost(-scalarParticle.p3/scalarParticle.energy)
        isrParticles_hadron = isrParticles_hadron.boost(-scalarParticle.p3/scalarParticle.energy)
        scalarParticle = scalarParticle.boost(-scalarParticle.p3/scalarParticle.energy)

        reco_jets = reco_jets.boost(-scalarParticle.p3/scalarParticle.energy)

    particles = np.concatenate([fromScalarParticles_e,fromScalarParticles_mu,fromScalarParticles_gamma,fromScalarParticles_pi,fromScalarParticles_hadron ,isrParticles_e , isrParticles_mu , isrParticles_gamma , isrParticles_pi, isrParticles_hadron], axis=0)  

    # Initialize plotting
    if ax is None:
        fig = plt.figure(figsize=(8,8))
        ax = plt.gca()

    # Plot parameters
    ax.set_xlim(-np.pi, np.pi)
    ax.set_ylim(eta_min, eta_max)
    ax.set_xlabel(r'$\phi$', fontsize=18)
    ax.set_ylabel(r'$\eta$', fontsize=18)
    ax.tick_params(axis='both', which='major', labelsize=12)

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

    # Add scatters to figure
    ax.scatter(fromScalarParticles_e.phi, fromScalarParticles_e.eta,
               s=scale(fromScalarParticles_e,scalarParticle),
               c='xkcd:blue', marker='o')
    ax.scatter(fromScalarParticles_mu.phi, fromScalarParticles_mu.eta,
               s=scale(fromScalarParticles_mu,scalarParticle),
               c='xkcd:blue', marker='v')
    ax.scatter(fromScalarParticles_gamma.phi, fromScalarParticles_gamma.eta,
               s=scale(fromScalarParticles_gamma,scalarParticle),
               c='xkcd:blue', marker='s')
    ax.scatter(fromScalarParticles_pi.phi, fromScalarParticles_pi.eta,
               s=scale(fromScalarParticles_pi,scalarParticle),
               c='xkcd:blue', marker='P')
    ax.scatter(fromScalarParticles_hadron.phi, fromScalarParticles_hadron.eta,
               s=scale(fromScalarParticles_hadron,scalarParticle),
               c='xkcd:blue', marker='*')
    ax.scatter(isrParticles_e.phi, isrParticles_e.eta,
               s=scale(isrParticles_e,scalarParticle),
               c='xkcd:magenta', marker='o')
    ax.scatter(isrParticles_mu.phi, isrParticles_mu.eta,
               s=scale(isrParticles_mu,scalarParticle),
               c='xkcd:magenta', marker='v')
    ax.scatter(isrParticles_gamma.phi, isrParticles_gamma.eta,
               s=scale(isrParticles_gamma,scalarParticle),
               c='xkcd:magenta', marker='s')
    ax.scatter(isrParticles_pi.phi, isrParticles_pi.eta,
               s=scale(isrParticles_pi,scalarParticle),
               c='xkcd:magenta', marker='P')
    ax.scatter(isrParticles_hadron.phi, isrParticles_hadron.eta,
               s=scale(isrParticles_hadron,scalarParticle),
               c='xkcd:magenta', marker='*')

    # Add the scalar mediator to the plot
    ax.scatter(scalarParticle.phi, scalarParticle.eta,
               s=1.*scale(scalarParticle,scalarParticle), facecolors='none',
               edgecolors='r')

    
    for rec in [0,1]:
      if rec:
        parts = tracks
        reco_title = "Reco"
      else:
        parts = particles
        reco_title = "Gen"
      for ja in [-1,0,1]:
        if ja == -1:
          jet_algo_title = "AKT"
        elif ja == 0:
          jet_algo_title = "CA"
        elif ja == 1:
          jet_algo_title = "KT"
        for Ri in [0.8,1,1.5,2]:
          ax1,leading_jetpt,subleading_jetpt = add_jets(ax,ja,Ri,parts)


#          # Create custom legends
#          # Legend 1 is particle type
#          line1 = ax1.scatter([-100], [-100], label='$e$',marker='o', c='xkcd:black')
#          line2 = ax1.scatter([-100], [-100], label='$\mu$', marker='v', c='xkcd:black')
#          line3 = ax1.scatter([-100], [-100], label='$\gamma$', marker='s', c='xkcd:black')
#          line4 = ax1.scatter([-100], [-100], label='$\pi$', marker='P', c='xkcd:black')
#          line5 = ax1.scatter([-100], [-100], label='other hadron', marker='*', c='xkcd:black')
#          line6 = ax1.scatter([-100], [-100], label='Scalar mediator', marker='o',
#                             facecolors='none', edgecolors='r')
#          line7 = ax1.scatter([-100], [-100], label='leading jet (pt:%f)'%(leading_jetpt), marker='o',c='xkcd:red')
#          line8 = ax1.scatter([-100], [-100], label='subleading jet (pt:%f)'%(subleading_jetpt), marker='o',c='xkcd:green')
#          line9 = ax1.scatter([-100], [-100], label='jet (pt>100)', marker='o',c='xkcd:cyan')
#          line10 = ax1.scatter([-100], [-100], label='jet (30< pt < 100)', marker='o',c='xkcd:blue')
#          first_legend = plt.legend(handles=[line1, line2, line3, line4, line5, line6,line7, line8,line9,line10],
#                                    loc='upper right', fontsize=12)
#          ax1.add_artist(first_legend)

          # Legend 2 is about particle origin
          blue_patch = mpatches.Patch(color='xkcd:blue', label='from scalar')
          magenta_patch = mpatches.Patch(color='xkcd:magenta', label='not from scalar')
          plt.legend(handles=[blue_patch, magenta_patch],loc='upper left')

          # build a rectangle in axes coords
          left, width = .0, 1.
          bottom, height = .0, 1.
          center = left + width/2.
          right = left + width
          top = bottom + height

          # axes coordinates are 0,0 is bottom left and 1,1 is upper right
          p = mpatches.Rectangle((left, bottom), width, height,
              fill=False, transform=ax1.transAxes, clip_on=False)

          ax1.add_patch(p)
          # Print event number
          ax1.text(left, top, 'Event %d'%ievt, horizontalalignment='left',
                  verticalalignment='bottom', transform=ax1.transAxes, fontsize=12)
          # Print sample details
          ax1.text(right, top, 'mMed=%d$\,$GeV,mDark=%d$\,$GeV,T=%d$\,$K,%s'%(mMed,mDark,temp,decayMode),
                  horizontalalignment='right', verticalalignment='bottom',
                  transform=ax1.transAxes, fontsize=12)
          # Print details about cuts
          ax1.text(left+0.02, bottom+0.01, 'Final particles have $P_{T}>1\,$GeV',
                  horizontalalignment='left', verticalalignment='bottom',
                  transform=ax1.transAxes, fontsize=12)
          # Print details of scalar mediator
          ax1.text(left+0.02, bottom+0.05, 'Scalar mediator $P_{T}=%d\,$GeV'%(scalarParticle.pt),
                  horizontalalignment='left', verticalalignment='bottom',
                  transform=ax1.transAxes, fontsize=12)
          fig.savefig('Results/%s/%s_%d/mMed%d_mDark%d_temp%d_decay-%s_Event%d.png'%(reco_title,jet_algo_title,int(Ri*10),mMed, mDark, temp, decayMode,ievt))
          fig.clear()
          plt.close(fig)


# The program runs through this loop
# if you want multiple plots then change multi to True
# Warning: multiplot funtion is currently broken
multi = False
if multi:
    figure, axs = plt.subplots(nrows=2, ncols=3, figsize=(16,8))
    axs_flat = axs.flatten()

boost = False
event = 0
#for i in range(event,6+event):
for i in [24,40,45,47]:
    print("Event %s"%(i))
    ht_pass = pass_HT(i)
    if not ht_pass:
      continue
    #if GenParticles.counts[i] == 0: 
    #  continue
    if multi == False:
#        for rec in [0,1]:
#          for ja in [-1,0,1]:
#            for Ri in [0.8,1,1.5,2]:
        plot(i,boost=boost)#,reco=rec,jet_algo=ja,R=Ri)
        #plot(i,boost=boost,reco=False,jet_algo=-1,R=0.8)
        #plot(i,boost=boost,reco=False,jet_algo=-1,R=1)
        #plot(i,boost=boost,reco=False,jet_algo=-1,R=1.5)
        #plot(i,boost=boost,reco=False,jet_algo=-1,R=2)
        #plot(i,boost=boost,reco=False,jet_algo=0,R=0.8)
        #plot(i,boost=boost,reco=False,jet_algo=0,R=1)
        #plot(i,boost=boost,reco=False,jet_algo=0,R=1.5)
        #plot(i,boost=boost,reco=False,jet_algo=0,R=2)
        #plot(i,boost=boost,reco=False,jet_algo=1,R=0.8)
        #plot(i,boost=boost,reco=False,jet_algo=1,R=1)
        #plot(i,boost=boost,reco=False,jet_algo=1,R=1.5)
        #plot(i,boost=boost,reco=False,jet_algo=1,R=2)
        #break;
    else:
        plot(i,ax=axs_flat[i])

#plt.show();
