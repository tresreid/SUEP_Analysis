import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path



lumi = 59.74*1000

jets_1000 = []
jets_400 = []
jets_qcd = []
for qcd in [0,1,2]:
  xsecs = [311900,29070,5962,1207,119.9,25.24] # signal xsec are (125,34.8), (400,5.9), (750,0.5), (1000,0.17)
  files = [300,500,700,1000,1500,2000]
  if qcd == 1:
    for xsec,f in  zip(xsecs,files):
      #fname = "data/qcd%s_comparisonsx.txt"%f
      fname = "../macros/data/qcd_%s_v0.txt"%f
      try:
        firstfile = open(fname)
      except IOError:
        print("file %s doesn't exist"%fname)
        continue
      qcd_tit = "QCD"
      next(firstfile)
      events = 0
      for line in firstfile.readlines():
        events += xsec
        cols = line.rstrip().split(' ')
        #jets_qcd.append({"Event":int(cols[0]), "Jet_algo": cols[1], "R":float(cols[2]), "pt":float(cols[3]), "pt_m":float(cols[4]), "nTracks":int(cols[5]), "nTracks_m":int(cols[6]), "girth":float(cols[7]), "girth_m":float(cols[8]),"mass":float(cols[9]),"mass_m":float(cols[10]),"xsec":xsec*lumi/100000,"trackpt":float(cols[11]),"trackpt_m":float(cols[12]),"suep_tracks":int(0), "isr_tracks":int(0), "suep_tracks_m":int(0), "isr_tracks_m":int(0), "total_suep":int(1), "total_isr":int(0)}) 
      jets_qcd.append({"Event":int(cols[0]), "Jet_algo": int(cols[1]), "R":float(cols[2]),"jet_id":int(cols[3]), "pt":float(cols[4]), "eta":float(cols[5]), "phi":float(cols[6]), "nTracks":int(cols[7]), "girth":float(cols[8]), "mass":float(cols[9]),"xsec":xsec*lumi/10000,"trackpt":float(cols[10]),"suep_tracks":int(cols[11]), "isr_tracks":int(cols[12]), "total_suep":int(cols[13]), "total_isr":int(cols[14]),"NVtx":int(cols[15]),"Num_Interactions":int(cols[16])}) 
      print("QCD %d Events: %d"%(f,events))
      print(len(jets_qcd))
  elif qcd == 2:
    #firstfile = open('data/signal400_comparisons_allv2.txt')
    firstfile = open('../macros/data/sig_400_v0.txt')
    qcd_tit = "Sig400"
    next(firstfile)
    events = 0
    xsec = 5.9
    for line in firstfile.readlines():
      events += 1
      cols = line.rstrip().split(' ')
      #jets_400.append({"Event":int(cols[0]), "Jet_algo": cols[1], "R":float(cols[2]),"jet_id":int(cols[3]), "pt":float(cols[4]), "nTracks":int(cols[5]), "girth":float(cols[6]), "mass":float(cols[7]),"xsec":5.9*lumi/10000,"trackpt":float(cols[8]),"suep_tracks":int(cols[9]), "isr_tracks":int(cols[10]), "total_suep":int(cols[11]), "total_isr":int(cols[12]),"NVtx":int(cols[13]),"Num_Interactions":int(cols[14])}) 
      jets_400.append({"Event":int(cols[0]), "Jet_algo": int(cols[1]), "R":float(cols[2]),"jet_id":int(cols[3]), "pt":float(cols[4]), "eta":float(cols[5]), "phi":float(cols[6]), "nTracks":int(cols[7]), "girth":float(cols[8]), "mass":float(cols[9]),"xsec":xsec*lumi/10000,"trackpt":float(cols[10]),"suep_tracks":int(cols[11]), "isr_tracks":int(cols[12]), "total_suep":int(cols[13]), "total_isr":int(cols[14]),"NVtx":int(cols[15]),"Num_Interactions":int(cols[16])}) 
    print("Signal 400 Events: %d"%(events))
  elif qcd == 0:
    #firstfile = open('data/signal1000_comparisons_allv2.txt')
    firstfile = open('../macros/data/sig_1000_v0.txt')
    #firstfile = open('data/signal1000_comparisons.txt')
    qcd_tit = "Sig"
    next(firstfile)
    events = 0
    xsec = 0.17
    for line in firstfile.readlines():
      events += 1
      cols = line.rstrip().split(' ')
      jets_1000.append({"Event":int(cols[0]), "Jet_algo": int(cols[1]), "R":float(cols[2]),"jet_id":int(cols[3]), "pt":float(cols[4]), "eta":float(cols[5]), "phi":float(cols[6]), "nTracks":int(cols[7]), "girth":float(cols[8]), "mass":float(cols[9]),"xsec":xsec*lumi/10000,"trackpt":float(cols[10]),"suep_tracks":int(cols[11]), "isr_tracks":int(cols[12]), "total_suep":int(cols[13]), "total_isr":int(cols[14]),"NVtx":int(cols[15]),"Num_Interactions":int(cols[16])}) 
    print("Signal Events: %d"%(events))
  
df_jets1000 = pd.DataFrame(jets_1000)
df_jets400 = pd.DataFrame(jets_400)
df_jetsqcd = pd.DataFrame(jets_qcd)
df_jets1000['pt_rank'] = df_jets1000.groupby(["Event","Jet_algo","R"])["pt"].rank(method="first",ascending=False)
df_jets1000['multi_rank'] = df_jets1000.groupby(["Event","Jet_algo","R"])["nTracks"].rank(method="first",ascending=False)
df_jets400['pt_rank'] = df_jets400.groupby(["Event","Jet_algo","R"])["pt"].rank(method="first",ascending=False)
df_jets400['multi_rank'] = df_jets400.groupby(["Event","Jet_algo","R"])["nTracks"].rank(method="first",ascending=False)

df_jets1000["suep_frac"] = df_jets1000["suep_tracks"]/df_jets1000["total_suep"]
df_jets1000["suep_purity"] = df_jets1000["suep_tracks"]/(df_jets1000["isr_tracks"]+df_jets1000["suep_tracks"])
df_jets1000["is_suep"] = df_jets1000["suep_purity"].apply(lambda x: 1 if x > 0.8 else 0) 
df_jets400["suep_frac"] = df_jets400["suep_tracks"]/df_jets400["total_suep"]
df_jets400["suep_purity"] = df_jets400["suep_tracks"]/(df_jets400["isr_tracks"]+df_jets400["suep_tracks"])
df_jets400["is_suep"] = df_jets400["suep_purity"].apply(lambda x: 1 if x > 0.8 else 0) 
df_jetsqcd["suep_frac"] = 0
df_jetsqcd["suep_purity"] = 0
df_jetsqcd["is_suep"] = 0
#print(df_jets1000)
def get_sig(algo,R,var,steps):
  jet1_df = df_jets1000[(df_jets1000["Jet_algo"] == algo) & (df_jets1000["R"] == R) & (df_jets1000["is_suep"] == 1)]
  jet2_df = df_jets400[(df_jets400["Jet_algo"] == algo) & (df_jets400["R"] == R) & (df_jets400["is_suep"] == 1)]

  sig1 = []
  sig2 = []
  tot_sig1 = jet1_df["xsec"].sum() 
  tot_sig2 = jet2_df["xsec"].sum() 
  for i in steps: 
    if 'trackpt' in var:
      sig1.append(jet1_df[jet1_df[var] < i]["xsec"].sum()) 
      sig2.append(jet2_df[jet2_df[var] < i]["xsec"].sum()) 
    else:
      sig1.append(jet1_df[jet1_df[var] > i]["xsec"].sum()) 
      sig2.append(jet2_df[jet2_df[var] > i]["xsec"].sum()) 
  return(sig1,sig2,tot_sig1,tot_sig2)
def get_bkg(algo,R,var,steps):
  jet1_df = df_jets1000[(df_jets1000["Jet_algo"] == algo) & (df_jets1000["R"] == R) & (df_jets1000["is_suep"] == 0)]
  jet2_df = df_jets400[(df_jets400["Jet_algo"] == algo) & (df_jets400["R"] == R) & (df_jets400["is_suep"] == 0)]
  qcd_df = df_jetsqcd[(df_jetsqcd["Jet_algo"] == algo) & (df_jetsqcd["R"] == R)]

  isr1 = []
  isr2 = []
  qcd = []
  bkg1 = []
  bkg2 = []
  tot_isr1 = jet1_df["xsec"].sum() 
  tot_isr2 = jet2_df["xsec"].sum() 
  tot_qcd  = qcd_df["xsec"].sum() 
  tot_bkg1 = qcd_df["xsec"].sum() + jet1_df["xsec"].sum() 
  tot_bkg2 = qcd_df["xsec"].sum() + jet2_df["xsec"].sum() 
  for i in steps: 
    if 'trackpt' in var:
      isr1.append(jet1_df[jet1_df[var] < i]["xsec"].sum()) 
      isr2.append(jet2_df[jet2_df[var] < i]["xsec"].sum()) 
      qcd.append(qcd_df[qcd_df[var] < i]["xsec"].sum())
      bkg1.append((qcd_df[qcd_df[var] < i]["xsec"].sum())+(jet1_df[jet1_df[var] < i]["xsec"].sum()))
      bkg2.append((qcd_df[qcd_df[var] < i]["xsec"].sum())+(jet2_df[jet2_df[var] < i]["xsec"].sum()))
    else:
      isr1.append(jet1_df[jet1_df[var] > i]["xsec"].sum()) 
      isr2.append(jet2_df[jet2_df[var] > i]["xsec"].sum()) 
      qcd.append((qcd_df[qcd_df[var] > i]["xsec"].sum()))
      bkg1.append((qcd_df[qcd_df[var] > i]["xsec"].sum())+(jet1_df[jet1_df[var] > i]["xsec"].sum()))
      bkg2.append((qcd_df[qcd_df[var] > i]["xsec"].sum())+(jet2_df[jet2_df[var] > i]["xsec"].sum()))
  return(qcd,isr1,isr2,bkg1,bkg2,tot_qcd,tot_isr1,tot_isr2,tot_bkg1,tot_bkg2)
def make_all_sig(var,steps):
  algo_set = []
  #for algo in ["AKT","KT","CA"]:
  for algo in [-1,0,1]:
    R_set = []
    for R in [0.8,1.0,1.5,2.0]:
      R_set.append(get_sig(algo,R,var,steps))
    algo_set.append(R_set)
  return algo_set
def make_all_bkg(var,steps):
  algo_set = []
  #for algo in ["AKT","KT","CA"]:
  for algo in [-1,0,1]:
    R_set = []
    for R in [0.8,1.0,1.5,2.0]:
      R_set.append(get_bkg(algo,R,var,steps))
    algo_set.append(R_set)
  return algo_set

def make_eff_algo_roc(var,steps,xtitle):
  sig_set = make_all_sig(var,steps)
  bkg_set = make_all_bkg(var,steps)
  alg_name = ["AKT","KT","CA"]
  #algs = [-1,0,1]
  Rs = [0.8,1.0,1.5,2.0]
  sigs = ["Sig1000","Sig400"]
  bkgs = ["qcd","isr","bkg"]
  for bkg in [0,1,2]:
    for sig in [0,1]:
      for i in [0,1,2]:
        if bkg == 0:
          bkg_x=0
        if bkg == 1:
          bkg_x=1+sig
        if bkg == 2:
          bkg_x=3+sig
        fig, (ax1,ax2) = plt.subplots(1,2)
        fig.suptitle("%s %s %s: %s"%(sigs[sig],bkgs[bkg],alg_name[i],var))
        ax1.set_xlabel(xtitle)
        ax1.set_ylabel("S/sqrt(S+B)")


        ax1.errorbar(steps,sig_set[i][0][sig]/(np.sqrt(np.add(sig_set[i][0][sig],bkg_set[i][0][bkg_x]))),(sig_set[i][0][sig]/(np.sqrt(np.add(sig_set[i][0][sig],bkg_set[i][0][bkg_x]))))*(np.sqrt(np.add(np.reciprocal(sig_set[i][0][sig]),np.reciprocal(4*np.add(sig_set[i][0][sig],bkg_set[i][0][bkg_x]))))),errorevery=int(len(steps)/20),ecolor='lightblue', label="R=0.8", color="blue")
        ax1.errorbar(steps,sig_set[i][1][sig]/(np.sqrt(np.add(sig_set[i][1][sig],bkg_set[i][1][bkg_x]))),(sig_set[i][1][sig]/(np.sqrt(np.add(sig_set[i][1][sig],bkg_set[i][1][bkg_x]))))*(np.sqrt(np.add(np.reciprocal(sig_set[i][1][sig]),np.reciprocal(4*np.add(sig_set[i][1][sig],bkg_set[i][1][bkg_x]))))),errorevery=int(len(steps)/20),ecolor='coral', label="R=1.0", color="orange")
        ax1.errorbar(steps,sig_set[i][2][sig]/(np.sqrt(np.add(sig_set[i][2][sig],bkg_set[i][2][bkg_x]))),(sig_set[i][2][sig]/(np.sqrt(np.add(sig_set[i][2][sig],bkg_set[i][2][bkg_x]))))*(np.sqrt(np.add(np.reciprocal(sig_set[i][2][sig]),np.reciprocal(4*np.add(sig_set[i][2][sig],bkg_set[i][2][bkg_x]))))),errorevery=int(len(steps)/20),ecolor='lightgreen', label="R=1.5", color="green")
        ax1.errorbar(steps,sig_set[i][3][sig]/(np.sqrt(np.add(sig_set[i][3][sig],bkg_set[i][3][bkg_x]))),(sig_set[i][3][sig]/(np.sqrt(np.add(sig_set[i][3][sig],bkg_set[i][3][bkg_x]))))*(np.sqrt(np.add(np.reciprocal(sig_set[i][3][sig]),np.reciprocal(4*np.add(sig_set[i][3][sig],bkg_set[i][3][bkg_x]))))),errorevery=int(len(steps)/20),ecolor='indianred', label="R=2.0", color="red")
        ax1.legend(loc="upper left")

        ax2.set_xlabel(xtitle)
        ax2.set_ylabel("Efficiency")
        ax2.errorbar(steps,sig_set[i][0][sig]/sig_set[i][0][sig+2],(sig_set[i][0][sig]/sig_set[i][0][sig+2])*np.sqrt(np.add(np.reciprocal(sig_set[i][0][sig]),np.reciprocal(sig_set[i][0][sig+2]))), label="Sig R=0.8",errorevery=int(len(steps)/10),ecolor='lightblue', color="blue", linestyle="-")
        ax2.errorbar(steps,sig_set[i][1][sig]/sig_set[i][1][sig+2],(sig_set[i][1][sig]/sig_set[i][1][sig+2])*np.sqrt(np.add(np.reciprocal(sig_set[i][1][sig]),np.reciprocal(sig_set[i][1][sig+2]))), label="Sig R=1.0",errorevery=int(len(steps)/10),ecolor='coral', color="orange", linestyle="-")
        ax2.errorbar(steps,sig_set[i][2][sig]/sig_set[i][2][sig+2],(sig_set[i][2][sig]/sig_set[i][2][sig+2])*np.sqrt(np.add(np.reciprocal(sig_set[i][2][sig]),np.reciprocal(sig_set[i][2][sig+2]))), label="Sig R=1.5",errorevery=int(len(steps)/10),ecolor='lightgreen', color="green", linestyle="-")
        ax2.errorbar(steps,sig_set[i][3][sig]/sig_set[i][3][sig+2],(sig_set[i][3][sig]/sig_set[i][3][sig+2])*np.sqrt(np.add(np.reciprocal(sig_set[i][3][sig]),np.reciprocal(sig_set[i][3][sig+2]))), label="Sig R=2.0",errorevery=int(len(steps)/10),ecolor='indianred', color="red", linestyle="-")


        ax2.errorbar(steps,bkg_set[i][0][bkg_x]/bkg_set[i][0][bkg_x+5],(bkg_set[i][0][bkg_x]/bkg_set[i][0][bkg_x+5])*np.sqrt(np.add(np.reciprocal(bkg_set[i][0][bkg_x]),np.reciprocal(bkg_set[i][0][bkg_x+5]))), errorevery=int(len(steps)/10),ecolor='lightblue', label="Bkg R=0.8", color="blue", linestyle=":")
        ax2.errorbar(steps,bkg_set[i][1][bkg_x]/bkg_set[i][1][bkg_x+5],(bkg_set[i][1][bkg_x]/bkg_set[i][1][bkg_x+5])*np.sqrt(np.add(np.reciprocal(bkg_set[i][1][bkg_x]),np.reciprocal(bkg_set[i][1][bkg_x+5]))), errorevery=int(len(steps)/10),ecolor='coral', label="Bkg R=1.0", color="orange", linestyle=":")
        ax2.errorbar(steps,bkg_set[i][2][bkg_x]/bkg_set[i][2][bkg_x+5],(bkg_set[i][2][bkg_x]/bkg_set[i][2][bkg_x+5])*np.sqrt(np.add(np.reciprocal(bkg_set[i][2][bkg_x]),np.reciprocal(bkg_set[i][2][bkg_x+5]))), errorevery=int(len(steps)/10),ecolor='lightgreen', label="Bkg R=1.5", color="green", linestyle=":")
        ax2.errorbar(steps,bkg_set[i][3][bkg_x]/bkg_set[i][3][bkg_x+5],(bkg_set[i][3][bkg_x]/bkg_set[i][3][bkg_x+5])*np.sqrt(np.add(np.reciprocal(bkg_set[i][3][bkg_x]),np.reciprocal(bkg_set[i][3][bkg_x+5]))), errorevery=int(len(steps)/10),ecolor='indianred', label="Bkg R=2.0", color="red", linestyle=":")
        ax2.legend(loc="upper right")
        fig.tight_layout()
        Path("Plots/ROC/bkg_%s/%s"%(bkgs[bkg],var)).mkdir(parents=True,exist_ok=True)
        fig.savefig("Plots/ROC/bkg_%s/%s/%s_%s.png"%(bkgs[bkg],var,sigs[sig],alg_name[i]))
        plt.close()

      for i in [0,1,2,3]:
        fig, (ax1,ax2) = plt.subplots(1,2)
        fig.suptitle("%s R%s: %s"%(sigs[sig],Rs[i],var))
        ax1.set_xlabel(xtitle)
        ax1.set_ylabel("S/sqrt(S+B)")
        ax1.errorbar(steps,sig_set[0][i][sig]/(np.sqrt(np.add(sig_set[0][i][sig],bkg_set[0][i][bkg_x]))),(sig_set[0][i][sig]/(np.sqrt(np.add(sig_set[0][i][sig],bkg_set[0][i][bkg_x]))))*(np.sqrt(np.add(np.reciprocal(sig_set[0][i][sig]),np.reciprocal(4*np.add(sig_set[0][i][sig],bkg_set[0][i][bkg_x]))))),errorevery=int(len(steps)/20),ecolor='lightblue', label="AKT", color="blue")
        ax1.errorbar(steps,sig_set[1][i][sig]/(np.sqrt(np.add(sig_set[1][i][sig],bkg_set[1][i][bkg_x]))),(sig_set[1][i][sig]/(np.sqrt(np.add(sig_set[1][i][sig],bkg_set[1][i][bkg_x]))))*(np.sqrt(np.add(np.reciprocal(sig_set[1][i][sig]),np.reciprocal(4*np.add(sig_set[1][i][sig],bkg_set[1][i][bkg_x]))))),errorevery=int(len(steps)/20),ecolor='coral', label="KT", color="orange")
        ax1.errorbar(steps,sig_set[2][i][sig]/(np.sqrt(np.add(sig_set[2][i][sig],bkg_set[2][i][bkg_x]))),(sig_set[2][i][sig]/(np.sqrt(np.add(sig_set[2][i][sig],bkg_set[2][i][bkg_x]))))*(np.sqrt(np.add(np.reciprocal(sig_set[2][i][sig]),np.reciprocal(4*np.add(sig_set[2][i][sig],bkg_set[2][i][bkg_x]))))),errorevery=int(len(steps)/20),ecolor='lightgreen', label="CA", color="green")
        ax1.legend(loc="upper left")

        ax2.set_xlabel(xtitle)
        ax2.set_ylabel("Efficiency")
        ax2.errorbar(steps,sig_set[0][i][sig]/sig_set[0][i][sig+2],(sig_set[0][i][sig]/sig_set[0][i][sig+2])*np.sqrt(np.add(np.reciprocal(sig_set[0][i][sig]),np.reciprocal(sig_set[0][i][sig+2]))), errorevery=int(len(steps)/10),ecolor='lightblue', label="Sig AKT", color="blue", linestyle="-")
        ax2.errorbar(steps,sig_set[1][i][sig]/sig_set[1][i][sig+2],(sig_set[1][i][sig]/sig_set[1][i][sig+2])*np.sqrt(np.add(np.reciprocal(sig_set[1][i][sig]),np.reciprocal(sig_set[1][i][sig+2]))), errorevery=int(len(steps)/10),ecolor='coral', label="Sig KT", color="orange", linestyle="-")
        ax2.errorbar(steps,sig_set[2][i][sig]/sig_set[2][i][sig+2],(sig_set[2][i][sig]/sig_set[2][i][sig+2])*np.sqrt(np.add(np.reciprocal(sig_set[2][i][sig]),np.reciprocal(sig_set[2][i][sig+2]))), errorevery=int(len(steps)/10),ecolor='lightgreen', label="Sig CA", color="green", linestyle="-")

        ax2.errorbar(steps,bkg_set[0][i][bkg_x]/bkg_set[0][i][bkg_x+5],(bkg_set[0][i][bkg_x]/bkg_set[0][i][bkg_x+5])*np.sqrt(np.add(np.reciprocal(bkg_set[0][i][bkg_x]),np.reciprocal(bkg_set[0][i][bkg_x+5]))), errorevery=int(len(steps)/10),ecolor='lightblue', label="Bkg AKT", color="blue", linestyle=":")
        ax2.errorbar(steps,bkg_set[1][i][bkg_x]/bkg_set[1][i][bkg_x+5],(bkg_set[1][i][bkg_x]/bkg_set[1][i][bkg_x+5])*np.sqrt(np.add(np.reciprocal(bkg_set[1][i][bkg_x]),np.reciprocal(bkg_set[1][i][bkg_x+5]))), errorevery=int(len(steps)/10),ecolor='coral', label="Bkg KT", color="orange", linestyle=":")
        ax2.errorbar(steps,bkg_set[2][i][bkg_x]/bkg_set[2][i][bkg_x+5],(bkg_set[2][i][bkg_x]/bkg_set[2][i][bkg_x+5])*np.sqrt(np.add(np.reciprocal(bkg_set[2][i][bkg_x]),np.reciprocal(bkg_set[2][i][bkg_x+5]))), errorevery=int(len(steps)/10),ecolor='lightgreen', label="Bkg CA", color="green", linestyle=":")
        ax2.legend(loc="upper right")
        fig.tight_layout()
        Path("Plots/ROC/bkg_%s/%s"%(bkgs[bkg],var)).mkdir(parents=True,exist_ok=True)
        fig.savefig("Plots/ROC/bkg_%s/%s/%s_%s.png"%(bkgs[bkg],var,sigs[sig],Rs[i]))
        plt.close()


def make_eff_algo_combo(var,steps,xtitle):
  sig_set = make_all_sig(var,steps)
  bkg_set = make_all_bkg(var,steps)
  alg_name = ["AKT","KT","CA"]
  algs = [-1,0,1]
  Rs = [0.8,1.0,1.5,2.0]
  sigs = ["Sig1000","Sig400"]
  bkgs = ["qcd","isr","bkg"]
  for bkg in [0,1,2]:
    for sig in [0,1]:
      for i in [0,1,2]:
        for r in [0,1,2,3]:
          if sig == 0:
            jet1_df = df_jets1000[(df_jets1000["Jet_algo"] == algs[i]) & (df_jets1000["R"] == Rs[r]) & (df_jets1000["is_suep"] == 1)]
            isr_df = df_jets1000[(df_jets1000["Jet_algo"] == algs[i]) & (df_jets1000["R"] == Rs[r]) & (df_jets1000["is_suep"] == 0)]
          if sig == 1:
            jet1_df = df_jets400[(df_jets400["Jet_algo"] == algs[i]) & (df_jets400["R"] == Rs[r]) & (df_jets400["is_suep"] == 1)]
            isr_df = df_jets400[(df_jets400["Jet_algo"] == algs[i]) & (df_jets400["R"] == Rs[r]) & (df_jets400["is_suep"] == 0)]
          if bkg == 0:
            bkg_df = df_jetsqcd[(df_jetsqcd["Jet_algo"] == algs[i]) & (df_jetsqcd["R"] == Rs[r])] 
            bkg_x=0
          if bkg == 1:
            bkg_df = isr_df
            bkg_x=1+sig
          if bkg == 2:
            bkg_df = df_jetsqcd[(df_jetsqcd["Jet_algo"] == algs[i]) & (df_jetsqcd["R"] == Rs[r])].append(isr_df) 
            bkg_x=3+sig

          fig, (ax1,ax2) = plt.subplots(1,2)
          fig.suptitle("%s %s %s %s: %s"%(sigs[sig],bkgs[bkg],alg_name[i],Rs[r],var))
          ax1.set_xlabel(xtitle)
          ax1.set_ylabel("S/sqrt(S+B)")


          ax1.errorbar(steps,sig_set[i][r][sig]/(np.sqrt(np.add(sig_set[i][r][sig],bkg_set[i][r][bkg_x]))),(sig_set[i][r][sig]/(np.sqrt(np.add(sig_set[i][r][sig],bkg_set[i][r][bkg_x]))))*(np.sqrt(np.add(np.reciprocal(sig_set[i][r][sig]),np.reciprocal(4*np.add(sig_set[i][r][sig],bkg_set[i][r][bkg_x]))))),errorevery=int(len(steps)/20),ecolor='lightblue', label="R=%s"%Rs[r], color="blue")
          ax1x = ax1.twinx()
          plt.yscale('log')
          ax1x.set_ylabel('Events')
          ax1x.hist(jet1_df[var],range=[steps[0],steps[-1]],weights=jet1_df['xsec'],bins=20,alpha=0.2,color='blue',label=sigs[sig])
          ax1x.hist(bkg_df[var],range=[steps[0],steps[-1]],weights=bkg_df['xsec'],bins=20,alpha=0.2,color='red',label=bkgs[bkg])
          ax1x.legend(loc="upper right")

          ax2.set_xlabel(xtitle)
          ax2.set_ylabel("Efficiency")
          ax2.errorbar(steps,sig_set[i][r][sig]/sig_set[i][r][sig+2],(sig_set[i][r][sig]/sig_set[i][r][sig+2])*np.sqrt(np.add(np.reciprocal(sig_set[i][r][sig]),np.reciprocal(sig_set[i][r][sig+2]))), label="Sig R=%s"%Rs[r],errorevery=int(len(steps)/10),ecolor='lightblue', color="blue", linestyle="-")
          ax2.errorbar(steps,bkg_set[i][r][bkg_x]/bkg_set[i][r][bkg_x+5],(bkg_set[i][r][bkg_x]/bkg_set[i][r][bkg_x+5])*np.sqrt(np.add(np.reciprocal(bkg_set[i][r][bkg_x]),np.reciprocal(bkg_set[i][r][bkg_x+5]))), errorevery=int(len(steps)/10),ecolor='lightblue', label="Bkg R=%s"%Rs[r], color="blue", linestyle=":")
          ax2x = ax2.twinx()
          ax2x.set_ylabel('A.U.')
          ax2x.hist(jet1_df[var],range=[steps[0],steps[-1]],weights=np.ones_like(jet1_df[var])/len(jet1_df[var]),bins=20,alpha=0.2,color='blue')
          ax2x.hist(bkg_df[var],range=[steps[0],steps[-1]],weights=np.ones_like(bkg_df[var])/len(bkg_df[var]),bins=20,alpha=0.2,color='red')
          ax2.legend(loc="upper right")
          fig.tight_layout()
          Path("Plots/COMBO/bkg_%s/%s"%(bkgs[bkg],var)).mkdir(parents=True,exist_ok=True)
          fig.savefig("Plots/COMBO/bkg_%s/%s/%s_%s_%s.png"%(bkgs[bkg],var,sigs[sig],alg_name[i],Rs[r]))
          plt.close()

def make_eff_algo_dist(var,steps,xtitle):
  #algo_set = make_all(var,steps)
  alg_name = ["AKT","KT","CA"]
  algs = [-1,0,1]
  Rs = [0.8,1.0,1.5,2.0]
  sigs = ["Sig1000","Sig400","qcd","isr1000","isr400","bkg1000","bkg400"]
  for sig in [0,1,2,3,4,5,6]:
    if ((("suep" in var) or ("NVtx" in var) or ("Num" in var)) and sig >= 2):
      continue
    for i in [0,1,2]:
        if sig == 0:
          if "suep" in var:
            jet1_df_1 = df_jets1000[(df_jets1000["Jet_algo"] == algs[i]) & (df_jets1000["R"] == Rs[0])]
            jet1_df_2 = df_jets1000[(df_jets1000["Jet_algo"] == algs[i]) & (df_jets1000["R"] == Rs[1])]
            jet1_df_3 = df_jets1000[(df_jets1000["Jet_algo"] == algs[i]) & (df_jets1000["R"] == Rs[2])]
            jet1_df_4 = df_jets1000[(df_jets1000["Jet_algo"] == algs[i]) & (df_jets1000["R"] == Rs[3])]
          else:
            jet1_df_1 = df_jets1000[(df_jets1000["Jet_algo"] == algs[i]) & (df_jets1000["R"] == Rs[0]) & (df_jets1000["is_suep"] == 1)]
            jet1_df_2 = df_jets1000[(df_jets1000["Jet_algo"] == algs[i]) & (df_jets1000["R"] == Rs[1]) & (df_jets1000["is_suep"] == 1)]
            jet1_df_3 = df_jets1000[(df_jets1000["Jet_algo"] == algs[i]) & (df_jets1000["R"] == Rs[2]) & (df_jets1000["is_suep"] == 1)]
            jet1_df_4 = df_jets1000[(df_jets1000["Jet_algo"] == algs[i]) & (df_jets1000["R"] == Rs[3]) & (df_jets1000["is_suep"] == 1)]
        if sig == 1:
          if "suep" in var:
            jet1_df_1 = df_jets400[(df_jets400["Jet_algo"] == algs[i]) & (df_jets400["R"] == Rs[0])]
            jet1_df_2 = df_jets400[(df_jets400["Jet_algo"] == algs[i]) & (df_jets400["R"] == Rs[1])]
            jet1_df_3 = df_jets400[(df_jets400["Jet_algo"] == algs[i]) & (df_jets400["R"] == Rs[2])]
            jet1_df_4 = df_jets400[(df_jets400["Jet_algo"] == algs[i]) & (df_jets400["R"] == Rs[3])]
          else:
            jet1_df_1 = df_jets400[(df_jets400["Jet_algo"] == algs[i]) & (df_jets400["R"] == Rs[0]) & (df_jets400["is_suep"] == 1)]
            jet1_df_2 = df_jets400[(df_jets400["Jet_algo"] == algs[i]) & (df_jets400["R"] == Rs[1]) & (df_jets400["is_suep"] == 1)]
            jet1_df_3 = df_jets400[(df_jets400["Jet_algo"] == algs[i]) & (df_jets400["R"] == Rs[2]) & (df_jets400["is_suep"] == 1)]
            jet1_df_4 = df_jets400[(df_jets400["Jet_algo"] == algs[i]) & (df_jets400["R"] == Rs[3]) & (df_jets400["is_suep"] == 1)]
        if sig == 2:
          jet1_df_1 = df_jetsqcd[(df_jetsqcd["Jet_algo"] == algs[i]) & (df_jetsqcd["R"] == Rs[0])]
          jet1_df_2 = df_jetsqcd[(df_jetsqcd["Jet_algo"] == algs[i]) & (df_jetsqcd["R"] == Rs[1])]
          jet1_df_3 = df_jetsqcd[(df_jetsqcd["Jet_algo"] == algs[i]) & (df_jetsqcd["R"] == Rs[2])]
          jet1_df_4 = df_jetsqcd[(df_jetsqcd["Jet_algo"] == algs[i]) & (df_jetsqcd["R"] == Rs[3])]
        if sig == 3:
          jet1_df_1 = df_jets1000[(df_jets1000["Jet_algo"] == algs[i]) & (df_jets1000["R"] == Rs[0]) & (df_jets1000["is_suep"] == 0)]
          jet1_df_2 = df_jets1000[(df_jets1000["Jet_algo"] == algs[i]) & (df_jets1000["R"] == Rs[1]) & (df_jets1000["is_suep"] == 0)]
          jet1_df_3 = df_jets1000[(df_jets1000["Jet_algo"] == algs[i]) & (df_jets1000["R"] == Rs[2]) & (df_jets1000["is_suep"] == 0)]
          jet1_df_4 = df_jets1000[(df_jets1000["Jet_algo"] == algs[i]) & (df_jets1000["R"] == Rs[3]) & (df_jets1000["is_suep"] == 0)]
        if sig == 4:
          jet1_df_1 = df_jets400[(df_jets400["Jet_algo"] == algs[i]) & (df_jets400["R"] == Rs[0]) & (df_jets400["is_suep"] == 0)]
          jet1_df_2 = df_jets400[(df_jets400["Jet_algo"] == algs[i]) & (df_jets400["R"] == Rs[1]) & (df_jets400["is_suep"] == 0)]
          jet1_df_3 = df_jets400[(df_jets400["Jet_algo"] == algs[i]) & (df_jets400["R"] == Rs[2]) & (df_jets400["is_suep"] == 0)]
          jet1_df_4 = df_jets400[(df_jets400["Jet_algo"] == algs[i]) & (df_jets400["R"] == Rs[3]) & (df_jets400["is_suep"] == 0)]
        if sig == 5:
          jet1_df_1 = df_jetsqcd[(df_jetsqcd["Jet_algo"] == algs[i]) & (df_jetsqcd["R"] == Rs[0])].append(df_jets1000[(df_jets1000["Jet_algo"] == algs[i]) & (df_jets1000["R"] == Rs[0]) & (df_jets1000["is_suep"] == 0)])
          jet1_df_2 = df_jetsqcd[(df_jetsqcd["Jet_algo"] == algs[i]) & (df_jetsqcd["R"] == Rs[1])].append(df_jets1000[(df_jets1000["Jet_algo"] == algs[i]) & (df_jets1000["R"] == Rs[1]) & (df_jets1000["is_suep"] == 0)])
          jet1_df_3 = df_jetsqcd[(df_jetsqcd["Jet_algo"] == algs[i]) & (df_jetsqcd["R"] == Rs[2])].append(df_jets1000[(df_jets1000["Jet_algo"] == algs[i]) & (df_jets1000["R"] == Rs[2]) & (df_jets1000["is_suep"] == 0)])
          jet1_df_4 = df_jetsqcd[(df_jetsqcd["Jet_algo"] == algs[i]) & (df_jetsqcd["R"] == Rs[3])].append(df_jets1000[(df_jets1000["Jet_algo"] == algs[i]) & (df_jets1000["R"] == Rs[3]) & (df_jets1000["is_suep"] == 0)])
        if sig == 6:
          jet1_df_1 = df_jetsqcd[(df_jetsqcd["Jet_algo"] == algs[i]) & (df_jetsqcd["R"] == Rs[0])].append(df_jets400[(df_jets400["Jet_algo"] == algs[i]) & (df_jets400["R"] == Rs[0]) & (df_jets400["is_suep"] == 0)])
          jet1_df_2 = df_jetsqcd[(df_jetsqcd["Jet_algo"] == algs[i]) & (df_jetsqcd["R"] == Rs[1])].append(df_jets400[(df_jets400["Jet_algo"] == algs[i]) & (df_jets400["R"] == Rs[1]) & (df_jets400["is_suep"] == 0)])
          jet1_df_3 = df_jetsqcd[(df_jetsqcd["Jet_algo"] == algs[i]) & (df_jetsqcd["R"] == Rs[2])].append(df_jets400[(df_jets400["Jet_algo"] == algs[i]) & (df_jets400["R"] == Rs[2]) & (df_jets400["is_suep"] == 0)])
          jet1_df_4 = df_jetsqcd[(df_jetsqcd["Jet_algo"] == algs[i]) & (df_jetsqcd["R"] == Rs[3])].append(df_jets400[(df_jets400["Jet_algo"] == algs[i]) & (df_jets400["R"] == Rs[3]) & (df_jets400["is_suep"] == 0)])

        fig, (ax1,ax2) = plt.subplots(1,2)
        fig.suptitle("%s %s: %s"%(sigs[sig],alg_name[i],var))
        ax1.set_xlabel(xtitle)
        plt.yscale('log')
        ax1.set_ylabel('Events')
        ax1.hist(jet1_df_1[var],range=[steps[0],steps[-1]],weights=jet1_df_1['xsec'],bins=20,histtype=u'step',color='blue',label="R=0.8")
        ax1.hist(jet1_df_2[var],range=[steps[0],steps[-1]],weights=jet1_df_2['xsec'],bins=20,histtype=u'step',color='orange',label="R=1.0")
        ax1.hist(jet1_df_3[var],range=[steps[0],steps[-1]],weights=jet1_df_3['xsec'],bins=20,histtype=u'step',color='green',label="R=1.5")
        ax1.hist(jet1_df_4[var],range=[steps[0],steps[-1]],weights=jet1_df_4['xsec'],bins=20,histtype=u'step',color='red',label="R=2.0")

        ax2.set_xlabel(xtitle)
        ax2.set_ylabel('A.U.')
        ax2.hist(jet1_df_1[var],range=[steps[0],steps[-1]],weights=np.ones_like(jet1_df_1[var])/len(jet1_df_1[var]),bins=20,histtype=u'step',color='blue',label="R=0.8")
        ax2.hist(jet1_df_2[var],range=[steps[0],steps[-1]],weights=np.ones_like(jet1_df_2[var])/len(jet1_df_2[var]),bins=20,histtype=u'step',color='orange',label="R=1.0")
        ax2.hist(jet1_df_3[var],range=[steps[0],steps[-1]],weights=np.ones_like(jet1_df_3[var])/len(jet1_df_3[var]),bins=20,histtype=u'step',color='green',label="R=1.5")
        ax2.hist(jet1_df_4[var],range=[steps[0],steps[-1]],weights=np.ones_like(jet1_df_4[var])/len(jet1_df_4[var]),bins=20,histtype=u'step',color='red',label="R=2.0")
        ax2.legend(loc="upper right")
        fig.tight_layout()
        Path("Plots/DIST/%s"%(var)).mkdir(parents=True,exist_ok=True)
        fig.savefig("Plots/DIST/%s/%s_%s.png"%(var,sigs[sig],alg_name[i]))
        plt.close()
    for i in [0,1,2,3]:
        if sig == 0:
          if "suep" in var:
            jet1_df_1 = df_jets1000[(df_jets1000["Jet_algo"] == algs[0]) & (df_jets1000["R"] == Rs[i])]
            jet1_df_2 = df_jets1000[(df_jets1000["Jet_algo"] == algs[1]) & (df_jets1000["R"] == Rs[i])]
            jet1_df_3 = df_jets1000[(df_jets1000["Jet_algo"] == algs[2]) & (df_jets1000["R"] == Rs[i])]
          else:
            jet1_df_1 = df_jets1000[(df_jets1000["Jet_algo"] == algs[0]) & (df_jets1000["R"] == Rs[i]) & (df_jets1000["is_suep"] == 1)]
            jet1_df_2 = df_jets1000[(df_jets1000["Jet_algo"] == algs[1]) & (df_jets1000["R"] == Rs[i]) & (df_jets1000["is_suep"] == 1)]
            jet1_df_3 = df_jets1000[(df_jets1000["Jet_algo"] == algs[2]) & (df_jets1000["R"] == Rs[i]) & (df_jets1000["is_suep"] == 1)]
        if sig == 1:
          if "suep" in var:
            jet1_df_1 = df_jets400[(df_jets400["Jet_algo"] == algs[0]) & (df_jets400["R"] == Rs[i])]
            jet1_df_2 = df_jets400[(df_jets400["Jet_algo"] == algs[1]) & (df_jets400["R"] == Rs[i])]
            jet1_df_3 = df_jets400[(df_jets400["Jet_algo"] == algs[2]) & (df_jets400["R"] == Rs[i])]
          else:
            jet1_df_1 = df_jets400[(df_jets400["Jet_algo"] == algs[0]) & (df_jets400["R"] == Rs[i]) & (df_jets400["is_suep"] == 1)]
            jet1_df_2 = df_jets400[(df_jets400["Jet_algo"] == algs[1]) & (df_jets400["R"] == Rs[i]) & (df_jets400["is_suep"] == 1)]
            jet1_df_3 = df_jets400[(df_jets400["Jet_algo"] == algs[2]) & (df_jets400["R"] == Rs[i]) & (df_jets400["is_suep"] == 1)]
        if sig == 2:
          jet1_df_1 = df_jetsqcd[(df_jetsqcd["Jet_algo"] == algs[0]) & (df_jetsqcd["R"] == Rs[i])]
          jet1_df_2 = df_jetsqcd[(df_jetsqcd["Jet_algo"] == algs[1]) & (df_jetsqcd["R"] == Rs[i])]
          jet1_df_3 = df_jetsqcd[(df_jetsqcd["Jet_algo"] == algs[2]) & (df_jetsqcd["R"] == Rs[i])]
        if sig == 3:
          jet1_df_1 = df_jets1000[(df_jets1000["Jet_algo"] == algs[0]) & (df_jets1000["R"] == Rs[i]) & (df_jets1000["is_suep"] == 0)]
          jet1_df_2 = df_jets1000[(df_jets1000["Jet_algo"] == algs[1]) & (df_jets1000["R"] == Rs[i]) & (df_jets1000["is_suep"] == 0)]
          jet1_df_3 = df_jets1000[(df_jets1000["Jet_algo"] == algs[2]) & (df_jets1000["R"] == Rs[i]) & (df_jets1000["is_suep"] == 0)]
        if sig == 4:
          jet1_df_1 = df_jets400[(df_jets400["Jet_algo"] == algs[0]) & (df_jets400["R"] == Rs[i]) & (df_jets400["is_suep"] == 0)]
          jet1_df_2 = df_jets400[(df_jets400["Jet_algo"] == algs[1]) & (df_jets400["R"] == Rs[i]) & (df_jets400["is_suep"] == 0)]
          jet1_df_3 = df_jets400[(df_jets400["Jet_algo"] == algs[2]) & (df_jets400["R"] == Rs[i]) & (df_jets400["is_suep"] == 0)]
        if sig == 5:
          jet1_df_1 = df_jetsqcd[(df_jetsqcd["Jet_algo"] == algs[0]) & (df_jetsqcd["R"] == Rs[i])].append(df_jets1000[(df_jets1000["Jet_algo"] == algs[0]) & (df_jets1000["R"] == Rs[i]) & (df_jets1000["is_suep"] == 0)])
          jet1_df_2 = df_jetsqcd[(df_jetsqcd["Jet_algo"] == algs[1]) & (df_jetsqcd["R"] == Rs[i])].append(df_jets1000[(df_jets1000["Jet_algo"] == algs[1]) & (df_jets1000["R"] == Rs[i]) & (df_jets1000["is_suep"] == 0)])
          jet1_df_3 = df_jetsqcd[(df_jetsqcd["Jet_algo"] == algs[2]) & (df_jetsqcd["R"] == Rs[i])].append(df_jets1000[(df_jets1000["Jet_algo"] == algs[2]) & (df_jets1000["R"] == Rs[i]) & (df_jets1000["is_suep"] == 0)])
        if sig == 6:
          jet1_df_1 = df_jetsqcd[(df_jetsqcd["Jet_algo"] == algs[0]) & (df_jetsqcd["R"] == Rs[i])].append(df_jets400[(df_jets400["Jet_algo"] == algs[0]) & (df_jets400["R"] == Rs[i]) & (df_jets400["is_suep"] == 0)])
          jet1_df_2 = df_jetsqcd[(df_jetsqcd["Jet_algo"] == algs[1]) & (df_jetsqcd["R"] == Rs[i])].append(df_jets400[(df_jets400["Jet_algo"] == algs[1]) & (df_jets400["R"] == Rs[i]) & (df_jets400["is_suep"] == 0)])
          jet1_df_3 = df_jetsqcd[(df_jetsqcd["Jet_algo"] == algs[2]) & (df_jetsqcd["R"] == Rs[i])].append(df_jets400[(df_jets400["Jet_algo"] == algs[2]) & (df_jets400["R"] == Rs[i]) & (df_jets400["is_suep"] == 0)])

        fig, (ax1,ax2) = plt.subplots(1,2)
        fig.suptitle("%s %s: %s"%(sigs[sig],Rs[i],var))
        ax1.set_xlabel(xtitle)
        plt.yscale('log')
        ax1.set_ylabel('Events')
        ax1.hist(jet1_df_1[var],range=[steps[0],steps[-1]],weights=jet1_df_1['xsec'],bins=20,histtype=u'step',color='blue',label="AKT")
        ax1.hist(jet1_df_2[var],range=[steps[0],steps[-1]],weights=jet1_df_2['xsec'],bins=20,histtype=u'step',color='orange',label="KT")
        ax1.hist(jet1_df_3[var],range=[steps[0],steps[-1]],weights=jet1_df_3['xsec'],bins=20,histtype=u'step',color='green',label="CA")

        ax2.set_xlabel(xtitle)
        ax2.set_ylabel('A.U.')
        ax2.hist(jet1_df_1[var],range=[steps[0],steps[-1]],weights=np.ones_like(jet1_df_1[var])/len(jet1_df_1[var]),bins=20,histtype=u'step',color='blue',label="AKT")
        ax2.hist(jet1_df_2[var],range=[steps[0],steps[-1]],weights=np.ones_like(jet1_df_2[var])/len(jet1_df_2[var]),bins=20,histtype=u'step',color='orange',label="KT")
        ax2.hist(jet1_df_3[var],range=[steps[0],steps[-1]],weights=np.ones_like(jet1_df_3[var])/len(jet1_df_3[var]),bins=20,histtype=u'step',color='green',label="CA")
        ax2.legend(loc="upper right")
        fig.tight_layout()
        Path("Plots/DIST/%s"%(var)).mkdir(parents=True,exist_ok=True)
        fig.savefig("Plots/DIST/%s/%s_%s.png"%(var,sigs[sig],Rs[i]))
        plt.close()

def make_eff_algo_dist_2d(var1,steps1,xtitle1,var2,steps2,xtitle2):
  #algo_set = make_all(var2,steps2)
  alg_name = ["AKT","KT","CA"]
  algs = [-1,0,1]
  Rs = [0.8,1.0,1.5,2.0]
  sigs = ["Sig1000","Sig400"]
  for sig in [0,1]:
    for i in [0,1,2]:
      for r in [0,1,2,3]:
        if sig == 0:
          jet1_df_1 = df_jets1000[(df_jets1000["Jet_algo"] == algs[i]) & (df_jets1000["R"] == Rs[r])]
        if sig == 1:
          jet1_df_1 = df_jets400[(df_jets400["Jet_algo"] == algs[i]) & (df_jets400["R"] == Rs[r])]

        def func(x,y):
          value = 0 
          bin_width = (steps1[-1] - steps1[0])/len(steps1)
          if "trackpt" in var2:
            value = value + jet1_df_1[(jet1_df_1[var2] < y) & (jet1_df_1[var1] >= x ) & (jet1_df_1[var1] < x+bin_width)]["xsec"].sum()
          else:
            value = value + jet1_df_1[(jet1_df_1[var2] > y) & (jet1_df_1[var1] >= x ) & (jet1_df_1[var1] < x+bin_width)]["xsec"].sum()
          if value == 0:
            return np.nan
          return value
    

        Xm,Ym = np.meshgrid(steps1,steps2)
        Z = np.zeros((len(steps2),len(steps1)))
        for xi,x in enumerate(steps1):
          for yi,y in enumerate(steps2):
            Z[yi,xi] = func(x,y)
        fig, (ax1,ax2) = plt.subplots(1,2)
        fig.set_size_inches(12.8,7.2)
        fig.suptitle("%s %s %s: %s vs %s"%(sigs[sig],alg_name[i],Rs[r],var1,var2))
        c = ax1.pcolormesh(Xm,Ym,Z,cmap="RdBu",shading="flat")
        fig.colorbar(c,ax=ax1) 
        ax1.set_xlabel(xtitle1)
        ax1.set_ylabel(xtitle2)
       
         
        d = ax2.hist2d(jet1_df_1[var1],jet1_df_1[var2],range=[[steps1[0],steps1[-1]],[steps2[0],steps2[-1]]],bins=[100,100],weights=jet1_df_1['xsec'],cmin=0.0000001)
        fig.colorbar(d[3],ax=ax2)
        ax2.set_xlabel(xtitle1)
        ax2.set_ylabel(xtitle2)
        Path("Plots/DIST2d/%s_%s"%(var1,var2)).mkdir(parents=True,exist_ok=True)
        fig.savefig("Plots/DIST2d/%s_%s/%s_%s_%s.png"%(var1,var2,sigs[sig],alg_name[i],Rs[r]))
        plt.close()

#print("Starting combo")
#make_eff_algo_combo("nTracks",range(0,500,1),"nTracks")
#make_eff_algo_combo("pt",range(0,1000,10),"Pt [GeV]")
#make_eff_algo_combo("trackpt",[0.005 * x for x in range(0,20,1)],"<track Pt/ jet pt>")
#make_eff_algo_combo("mass",range(0,1000,1),"Mass [GeV")
#make_eff_algo_combo("girth",[0.01 *x for x in range(0,200,1)],"$1/Pt \Sigma_{i} (Pt_{i} * \Delta R_{i,jet})$")
#print("Starting roc")
#make_eff_algo_roc("nTracks",range(0,500,1),"nTracks")
#make_eff_algo_roc("pt",range(0,1000,10),"Pt [GeV]")
#make_eff_algo_roc("trackpt",[0.005 * x for x in range(0,20,1)],"<track Pt/ jet pt>")
#make_eff_algo_roc("mass",range(0,1000,1),"Mass [GeV")
#make_eff_algo_roc("girth",[0.01 *x for x in range(0,200,1)],"$1/Pt \Sigma_{i} (Pt_{i} * \Delta R_{i,jet})$")
print("Starting dist")
make_eff_algo_dist("nTracks",range(0,500,1),"nTracks")
make_eff_algo_dist("pt",range(0,1000,10),"Pt [GeV]")
make_eff_algo_dist("trackpt",[0.005 * x for x in range(0,20,1)],"<track Pt/ jet pt>")
make_eff_algo_dist("mass",range(0,1000,1),"Mass [GeV")
make_eff_algo_dist("girth",[0.01 *x for x in range(0,200,1)],"$1/Pt \Sigma_{i} (Pt_{i} * \Delta R_{i,jet})$")
make_eff_algo_dist("suep_purity",[0.01*x for x in range(0,100,2)],"purity")
make_eff_algo_dist("suep_frac",[0.01*x for x in range(0,100,2)],"Suep Fraction")
make_eff_algo_dist("NVtx",range(0,100,1),"NVtx")
make_eff_algo_dist("Num_Interactions",range(0,100,1),"Num Interactions")
#print("starting purity 2d plots")
#make_eff_algo_dist_2d("suep_purity",[0.01*x for x in range(0,100,2)],"purity","nTracks",range(0,500,25),"nTracks")
#make_eff_algo_dist_2d("suep_purity",[0.01*x for x in range(0,100,2)],"purity","pt",range(0,1000,50),"Pt [GeV]")
#make_eff_algo_dist_2d("suep_purity",[0.01*x for x in range(0,100,2)],"purity","trackpt",[0.005 * x for x in range(0,20,1)],"<track Pt/ jet pt>")
#make_eff_algo_dist_2d("suep_purity",[0.01*x for x in range(0,100,2)],"purity","mass",range(0,1000,50),"Mass [GeV]")
#make_eff_algo_dist_2d("suep_purity",[0.01*x for x in range(0,100,2)],"purity","girth",[0.01 *x for x in range(0,200,20)],"$1/Pt \Sigma_{i} (Pt_{i} * \Delta R_{i,jet})$")
#make_eff_algo_dist_2d("suep_purity",[0.01*x for x in range(0,100,2)],"purity","NVtx",range(0,100,1),"NVtx")
#make_eff_algo_dist_2d("suep_purity",[0.01*x for x in range(0,100,2)],"purity","Num_Interactions",range(0,100,1),"Num Interactions")
#print("starting frac 2d plots")
#make_eff_algo_dist_2d("suep_frac",[0.01*x for x in range(0,100,2)],"Suep Fraction","nTracks",range(0,500,25),"nTracks")
#make_eff_algo_dist_2d("suep_frac",[0.01*x for x in range(0,100,2)],"Suep Fraction","pt",range(0,1000,50),"Pt [GeV]")
#make_eff_algo_dist_2d("suep_frac",[0.01*x for x in range(0,100,2)],"Suep Fraction","trackpt",[0.005 * x for x in range(0,20,1)],"<track Pt/ jet pt>")
#make_eff_algo_dist_2d("suep_frac",[0.01*x for x in range(0,100,2)],"Suep Fraction","mass",range(0,1000,50),"Mass [GeV]")
#make_eff_algo_dist_2d("suep_frac",[0.01*x for x in range(0,100,2)],"Suep Fraction","girth",[0.01 *x for x in range(0,200,20)],"$1/Pt \Sigma_{i} (Pt_{i} * \Delta R_{i,jet})$")
#make_eff_algo_dist_2d("suep_frac",[0.01*x for x in range(0,100,2)],"Suep Fraction","NVtx",range(0,100,1),"NVtx")
#make_eff_algo_dist_2d("suep_frac",[0.01*x for x in range(0,100,2)],"Suep Fraction","Num_Interactions",range(0,100,1),"Num Interactions")
#print("starting nvtx 2d plots")
#make_eff_algo_dist_2d("NVtx",range(0,100,1),"NVtx","nTracks",range(0,500,25),"nTracks")
#make_eff_algo_dist_2d("NVtx",range(0,100,1),"NVtx","pt",range(0,1000,50),"Pt [GeV]")
#make_eff_algo_dist_2d("NVtx",range(0,100,1),"NVtx","trackpt",[0.005 * x for x in range(0,20,1)],"<track Pt/ jet pt>")
#make_eff_algo_dist_2d("NVtx",range(0,100,1),"NVtx","mass",range(0,1000,50),"Mass [GeV]")
#make_eff_algo_dist_2d("NVtx",range(0,100,1),"NVtx","girth",[0.01 *x for x in range(0,200,20)],"$1/Pt \Sigma_{i} (Pt_{i} * \Delta R_{i,jet})$")
#print("starting num interactions 2d plots")
#make_eff_algo_dist_2d("Num_Interactions",range(0,100,1),"Num Interactions","nTracks",range(0,500,25),"nTracks")
#make_eff_algo_dist_2d("Num_Interactions",range(0,100,1),"Num Interactions","pt",range(0,1000,50),"Pt [GeV]")
#make_eff_algo_dist_2d("Num_Interactions",range(0,100,1),"Num Interactions","trackpt",[0.005 * x for x in range(0,20,1)],"<track Pt/ jet pt>")
#make_eff_algo_dist_2d("Num_Interactions",range(0,100,1),"Num Interactions","mass",range(0,1000,50),"Mass [GeV]")
#make_eff_algo_dist_2d("Num_Interactions",range(0,100,1),"Num Interactions","girth",[0.01 *x for x in range(0,200,20)],"$1/Pt \Sigma_{i} (Pt_{i} * \Delta R_{i,jet})$")
#
