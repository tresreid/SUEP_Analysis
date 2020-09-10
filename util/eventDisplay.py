import uproot
import coffea
from coffea import hist
import numpy as np
import seutils
from math import pi
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Get the file and import using uproot
mMed = 1000
mDark = 2
temp = 2
decayMode = 'darkPho'
#decayMode = 'darkPhoHad'
base = '/Users/chrispap/'
# xrootd is not working properly in Python3 :(
#base = 'root://cmseos.fnal.gov//store/user/kdipetri/SUEP/Production_v0.1/2018/NTUP/'
datasets = [base +
            'PrivateSamples.SUEP_2018_mMed-%d_mDark-%d_temp-%d_decay-%s'
            '_13TeV-pythia8_n-100_0_RA2AnalysisTree.root'%(mMed, mDark, temp, decayMode),
           ]
rootfile = datasets[0]
fin = uproot.open(rootfile)

# Attach the branches to numpy arrays
tree = fin[b'TreeMaker2'][b'PreSelection']
def get_branch(branchname):
    return tree.arrays(branchname)[branchname]

GenParticles = get_branch(b'GenParticles')
GenParticles_ParentId = get_branch(b'GenParticles_ParentId')
GenParticles_PdgId = get_branch(b'GenParticles_PdgId')
GenParticles_Status = get_branch(b'GenParticles_Status')

# Main plotting function
def plot(ievt, ax=None, boost=False):
    # Get the particles of ievt event
    genParticles = GenParticles[ievt]
    genParticles_ParentId = GenParticles_ParentId[ievt]
    genParticles_PdgId = GenParticles_PdgId[ievt]
    genParticles_Status = GenParticles_Status[ievt]
    # The last copy of the scalar mediator
    scalarParticle = genParticles[(genParticles_PdgId == 25) & (genParticles_Status == 62)]

    # Define mask arrays to select the desired particles
    finalParticles = (genParticles_Status == 1) & (genParticles.pt > 1)
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

    # Initialize plotting
    if ax is None:
        fig = plt.figure(figsize=(8,8))
        ax = plt.gca()

    # Plot parameters
    ax.set_xlim(-pi, pi)
    ax.set_ylim(-4, 4)
    ax.set_xlabel(r'$\phi$', fontsize=18)
    ax.set_ylabel(r'$\eta$', fontsize=18)
    ax.tick_params(axis='both', which='major', labelsize=12)

    # Function that sets the scaling for markers
    # Two methods for the moment: use tanh or square.
    # Scaling using just the energy is also an option
    def scale(particles, scalar, method=3):
        """Just to scale to a reasonable dot size"""
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

    # Create custom legends
    # Legend 1 is particle type
    line1 = ax.scatter([-100], [-100], label='$e$',marker='o', c='xkcd:black')
    line2 = ax.scatter([-100], [-100], label='$\mu$', marker='v', c='xkcd:black')
    line3 = ax.scatter([-100], [-100], label='$\gamma$', marker='s', c='xkcd:black')
    line4 = ax.scatter([-100], [-100], label='$\pi$', marker='P', c='xkcd:black')
    line5 = ax.scatter([-100], [-100], label='other hadron', marker='*', c='xkcd:black')
    line6 = ax.scatter([-100], [-100], label='Scalar mediator', marker='o',
                       facecolors='none', edgecolors='r')
    first_legend = plt.legend(handles=[line1, line2, line3, line4, line5, line6],
                              loc='upper right', fontsize=12)
    ax.add_artist(first_legend)

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
        fill=False, transform=ax.transAxes, clip_on=False)

    ax.add_patch(p)
    # Print event number
    ax.text(left, top, 'Event %d'%ievt, horizontalalignment='left',
            verticalalignment='bottom', transform=ax.transAxes, fontsize=12)
    # Print sample details
    ax.text(right, top, 'mMed=%d$\,$GeV,mDark=%d$\,$GeV,T=%d$\,$K,%s'%(mMed,mDark,temp,decayMode),
            horizontalalignment='right', verticalalignment='bottom',
            transform=ax.transAxes, fontsize=12)
    # Print details about cuts
    ax.text(left+0.02, bottom+0.01, 'Final particles have $P_{T}>1\,$GeV',
            horizontalalignment='left', verticalalignment='bottom',
            transform=ax.transAxes, fontsize=12)
    # Print details of scalar mediator
    ax.text(left+0.02, bottom+0.05, 'Scalar mediator $P_{T}=%d\,$GeV'%(scalarParticle.pt),
            horizontalalignment='left', verticalalignment='bottom',
            transform=ax.transAxes, fontsize=12)

    fig.savefig('Results/mMed%d_mDark%d_temp%d_decay-%s_Event%d.pdf'%(mMed, mDark, temp, decayMode,event))


# The program runs through this loop
# if you want multiple plots then change multi to True
# Warning: multiplot funtion is currently broken
multi = False
if multi:
    figure, axs = plt.subplots(nrows=2, ncols=3, figsize=(16,8))
    axs_flat = axs.flatten()

boost = False
event = 0
for i in range(event,6+event):
    if GenParticles.counts[i] == 0: continue
    if multi == False:
        plot(i,boost=boost)
        break;
    else:
        plot(i,ax=axs_flat[i])

plt.show();
