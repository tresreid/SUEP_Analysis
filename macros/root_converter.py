import root_pandas
import numpy as np
import pandas as pd
import gc

lumi = 59.74*1000

jets_1000 = []
jets_400 = []
jets_750 = []
jets_300 = []
jets_200 = []
jets_qcd = []
for qcd in [0,1,2,3,4,5]:
  xsecs = [311900,29070,5962,1207,119.9,25.24] # signal xsec are (125,34.8), (400,5.9), (750,0.5), (1000,0.17)
  files = [300,500,700,1000,1500,2000]
  if qcd == 0:
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
        jets_qcd.append({"Event":int(cols[0]), "Jet_algo": int(cols[1]), "R":float(cols[2]),"jet_id":int(cols[3]),
          "pt":float(cols[4]), "eta":float(cols[5]), "phi":float(cols[6]), "nTracks":int(cols[7]), "girth":float(cols[8]), "mass":float(cols[9]),"xsec":xsec*lumi/10000,
           "trackpt":float(cols[10]),"medpt":float(cols[11]),"suep_tracks":int(cols[12]), "isr_tracks":int(cols[13]), "total_suep":int(cols[14]), "total_isr":int(cols[15]),
            "NVtx":int(cols[16]),"Num_Interactions":int(cols[17]),"scalardR":float(cols[18]),"scalarpt":float(cols[19]),"scalareta":float(cols[20]),"scalarphi":float(cols[21]),
            "suep_ptwgt":float(cols[22]),"scalarmass":float(cols[23]),"beta":float(cols[24]),"scalarbeta":float(cols[25]),
            "pt_dispersion":float(cols[26]),"lesHouches":float(cols[27]),"thrust":float(cols[28]),
            "t1":float(cols[29]),"t2":float(cols[30]),"t3":float(cols[31]),"e2":float(cols[32]),"e3":float(cols[33])
            })
      print("QCD %d Events: %d"%(f,events))
      print(len(jets_qcd))
  elif qcd == 1:
    firstfile = open('../macros/data/sig_1000_v0.txt')
    qcd_tit = "Sig"
    next(firstfile)
    events = 0
    xsec = 0.17
    for line in firstfile.readlines():
      events += 1
      cols = line.rstrip().split(' ')
      jets_1000.append({"Event":int(cols[0]), "Jet_algo": int(cols[1]), "R":float(cols[2]),"jet_id":int(cols[3]),
      "pt":float(cols[4]), "eta":float(cols[5]), "phi":float(cols[6]), "nTracks":int(cols[7]), "girth":float(cols[8]), "mass":float(cols[9]),"xsec":xsec*lumi/10000,
           "trackpt":float(cols[10]),"medpt":float(cols[11]),"suep_tracks":int(cols[12]), "isr_tracks":int(cols[13]), "total_suep":int(cols[14]), "total_isr":int(cols[15]),
            "NVtx":int(cols[16]),"Num_Interactions":int(cols[17]),"scalardR":float(cols[18]),"scalarpt":float(cols[19]),"scalareta":float(cols[20]),"scalarphi":float(cols[21]),
            "suep_ptwgt":float(cols[22]),"scalarmass":float(cols[23]),"beta":float(cols[24]),"scalarbeta":float(cols[25]),
            "pt_dispersion":float(cols[26]),"lesHouches":float(cols[27]),"thrust":float(cols[28]),
            "t1":float(cols[29]),"t2":float(cols[30]),"t3":float(cols[31]),"e2":float(cols[32]),"e3":float(cols[33])
            })
    print("Signal Events: %d"%(events))
  elif qcd == 2:
    firstfile = open('../macros/data/sig_400_v0.txt')
    qcd_tit = "Sig400"
    next(firstfile)
    events = 0
    xsec = 5.9
    for line in firstfile.readlines():
      events += 1
      cols = line.rstrip().split(' ')
      jets_400.append({"Event":int(cols[0]), "Jet_algo": int(cols[1]), "R":float(cols[2]),"jet_id":int(cols[3]),
      "pt":float(cols[4]), "eta":float(cols[5]), "phi":float(cols[6]), "nTracks":int(cols[7]), "girth":float(cols[8]), "mass":float(cols[9]),"xsec":xsec*lumi/10000,
           "trackpt":float(cols[10]),"medpt":float(cols[11]),"suep_tracks":int(cols[12]), "isr_tracks":int(cols[13]), "total_suep":int(cols[14]), "total_isr":int(cols[15]),
            "NVtx":int(cols[16]),"Num_Interactions":int(cols[17]),"scalardR":float(cols[18]),"scalarpt":float(cols[19]),"scalareta":float(cols[20]),"scalarphi":float(cols[21]),
            "suep_ptwgt":float(cols[22]),"scalarmass":float(cols[23]),"beta":float(cols[24]),"scalarbeta":float(cols[25]),
            "pt_dispersion":float(cols[26]),"lesHouches":float(cols[27]),"thrust":float(cols[28]),
            "t1":float(cols[29]),"t2":float(cols[30]),"t3":float(cols[31]),"e2":float(cols[32]),"e3":float(cols[33])
            })
    print("Signal 400 Events: %d"%(events))
  elif qcd == 3:
    firstfile = open('../macros/data/sig_750_v0.txt')
    qcd_tit = "Sig750"
    next(firstfile)
    events = 0
    xsec = 0.5
    for line in firstfile.readlines():
      events += 1
      cols = line.rstrip().split(' ')
      jets_750.append({"Event":int(cols[0]), "Jet_algo": int(cols[1]), "R":float(cols[2]),"jet_id":int(cols[3]),
      "pt":float(cols[4]), "eta":float(cols[5]), "phi":float(cols[6]), "nTracks":int(cols[7]), "girth":float(cols[8]), "mass":float(cols[9]),"xsec":xsec*lumi/10000,
           "trackpt":float(cols[10]),"medpt":float(cols[11]),"suep_tracks":int(cols[12]), "isr_tracks":int(cols[13]), "total_suep":int(cols[14]), "total_isr":int(cols[15]),
            "NVtx":int(cols[16]),"Num_Interactions":int(cols[17]),"scalardR":float(cols[18]),"scalarpt":float(cols[19]),"scalareta":float(cols[20]),"scalarphi":float(cols[21]),
            "suep_ptwgt":float(cols[22]),"scalarmass":float(cols[23]),"beta":float(cols[24]),"scalarbeta":float(cols[25]),
            "pt_dispersion":float(cols[26]),"lesHouches":float(cols[27]),"thrust":float(cols[28]),
            "t1":float(cols[29]),"t2":float(cols[30]),"t3":float(cols[31]),"e2":float(cols[32]),"e3":float(cols[33])
            })
    print("Signal 750 Events: %d"%(events))
  elif qcd == 4:
    firstfile = open('../macros/data/sig_200_v0.txt')
    qcd_tit = "Sig200"
    next(firstfile)
    events = 0
    xsec = 13.6
    for line in firstfile.readlines():
      events += 1
      cols = line.rstrip().split(' ')
      jets_200.append({"Event":int(cols[0]), "Jet_algo": int(cols[1]), "R":float(cols[2]),"jet_id":int(cols[3]),
      "pt":float(cols[4]), "eta":float(cols[5]), "phi":float(cols[6]), "nTracks":int(cols[7]), "girth":float(cols[8]), "mass":float(cols[9]),"xsec":xsec*lumi/10000,
           "trackpt":float(cols[10]),"medpt":float(cols[11]),"suep_tracks":int(cols[12]), "isr_tracks":int(cols[13]), "total_suep":int(cols[14]), "total_isr":int(cols[15]),
            "NVtx":int(cols[16]),"Num_Interactions":int(cols[17]),"scalardR":float(cols[18]),"scalarpt":float(cols[19]),"scalareta":float(cols[20]),"scalarphi":float(cols[21]),
            "suep_ptwgt":float(cols[22]),"scalarmass":float(cols[23]),"beta":float(cols[24]),"scalarbeta":float(cols[25]),
            "pt_dispersion":float(cols[26]),"lesHouches":float(cols[27]),"thrust":float(cols[28]),
            "t1":float(cols[29]),"t2":float(cols[30]),"t3":float(cols[31]),"e2":float(cols[32]),"e3":float(cols[33])
            })
    print("Signal 200 Events: %d"%(events))
  elif qcd == 5:
    firstfile = open('../macros/data/sig_300_v0.txt')
    qcd_tit = "Sig300"
    next(firstfile)
    events = 0
    xsec = 8.9
    for line in firstfile.readlines():
      events += 1
      cols = line.rstrip().split(' ')
      jets_300.append({"Event":int(cols[0]), "Jet_algo": int(cols[1]), "R":float(cols[2]),"jet_id":int(cols[3]),
      "pt":float(cols[4]), "eta":float(cols[5]), "phi":float(cols[6]), "nTracks":int(cols[7]), "girth":float(cols[8]), "mass":float(cols[9]),"xsec":xsec*lumi/10000,
           "trackpt":float(cols[10]),"medpt":float(cols[11]),"suep_tracks":int(cols[12]), "isr_tracks":int(cols[13]), "total_suep":int(cols[14]), "total_isr":int(cols[15]),
            "NVtx":int(cols[16]),"Num_Interactions":int(cols[17]),"scalardR":float(cols[18]),"scalarpt":float(cols[19]),"scalareta":float(cols[20]),"scalarphi":float(cols[21]),
            "suep_ptwgt":float(cols[22]),"scalarmass":float(cols[23]),"beta":float(cols[24]),"scalarbeta":float(cols[25]),
            "pt_dispersion":float(cols[26]),"lesHouches":float(cols[27]),"thrust":float(cols[28]),
            "t1":float(cols[29]),"t2":float(cols[30]),"t3":float(cols[31]),"e2":float(cols[32]),"e3":float(cols[33])
            })
    print("Signal 300 Events: %d"%(events))

print("starting pandas")
df_jets1000 = pd.DataFrame(jets_1000)
df_jets400 = pd.DataFrame(jets_400)
df_jetsqcd = pd.DataFrame(jets_qcd)
df_jets750 = pd.DataFrame(jets_750)
df_jets300 = pd.DataFrame(jets_300)
df_jets200 = pd.DataFrame(jets_200)

gc.collect()
df_jets1000["suep_frac"] = df_jets1000["suep_tracks"]/df_jets1000["total_suep"]
df_jets1000["suep_purity"] = df_jets1000["suep_tracks"]/(df_jets1000["isr_tracks"]+df_jets1000["suep_tracks"])
df_jets1000["is_suep"] = df_jets1000["suep_purity"].apply(lambda x: 1 if x > 0.8 else 0)
df_jets400["suep_frac"] = df_jets400["suep_tracks"]/df_jets400["total_suep"]
df_jets400["suep_purity"] = df_jets400["suep_tracks"]/(df_jets400["isr_tracks"]+df_jets400["suep_tracks"])
df_jets400["is_suep"] = df_jets400["suep_purity"].apply(lambda x: 1 if x > 0.8 else 0)
df_jetsqcd["suep_frac"] = 0
df_jetsqcd["suep_purity"] = 0
df_jetsqcd["is_suep"] = 0
df_jets750["suep_frac"] = df_jets750["suep_tracks"]/df_jets750["total_suep"]
df_jets750["suep_purity"] = df_jets750["suep_tracks"]/(df_jets750["isr_tracks"]+df_jets750["suep_tracks"])
df_jets750["is_suep"] = df_jets750["suep_purity"].apply(lambda x: 1 if x > 0.8 else 0)
df_jets300["suep_frac"] = df_jets300["suep_tracks"]/df_jets300["total_suep"]
df_jets300["suep_purity"] = df_jets300["suep_tracks"]/(df_jets300["isr_tracks"]+df_jets300["suep_tracks"])
df_jets300["is_suep"] = df_jets300["suep_purity"].apply(lambda x: 1 if x > 0.8 else 0)
df_jets200["suep_frac"] = df_jets200["suep_tracks"]/df_jets200["total_suep"]
df_jets200["suep_purity"] = df_jets200["suep_tracks"]/(df_jets200["isr_tracks"]+df_jets200["suep_tracks"])
df_jets200["is_suep"] = df_jets200["suep_purity"].apply(lambda x: 1 if x > 0.8 else 0)

df_jets1000["phi"] = df_jets1000["phi"].apply(lambda x: x-2*np.pi if x> np.pi else x)
df_jets400["phi"] = df_jets400["phi"].apply(lambda x: x-2*np.pi if x> np.pi else x)
df_jetsqcd["phi"] = df_jetsqcd["phi"].apply(lambda x: x-2*np.pi if x> np.pi else x)
df_jets750["phi"] = df_jets750["phi"].apply(lambda x: x-2*np.pi if x> np.pi else x)
df_jets300["phi"] = df_jets300["phi"].apply(lambda x: x-2*np.pi if x> np.pi else x)
df_jets200["phi"] = df_jets200["phi"].apply(lambda x: x-2*np.pi if x> np.pi else x)

df_jets1000["pt_res"] = (df_jets1000["pt"] - df_jets1000["scalarpt"])/df_jets1000["scalarpt"]
df_jets1000["mass_res"] = (df_jets1000["mass"] - df_jets1000["scalarmass"])/df_jets1000["scalarmass"]
df_jets1000["beta_res"] = (df_jets1000["beta"] - df_jets1000["scalarbeta"])/df_jets1000["scalarbeta"]
df_jets400["pt_res"] = (df_jets400["pt"] - df_jets400["scalarpt"])/df_jets400["scalarpt"]
df_jets400["mass_res"] = (df_jets400["mass"] - df_jets400["scalarmass"])/df_jets400["scalarmass"]
df_jets400["beta_res"] = (df_jets400["beta"] - df_jets400["scalarbeta"])/df_jets400["scalarbeta"]
df_jetsqcd["pt_res"] = (df_jetsqcd["pt"] - df_jetsqcd["scalarpt"])/df_jetsqcd["scalarpt"]
df_jetsqcd["mass_res"] = (df_jetsqcd["mass"] - df_jetsqcd["scalarmass"])/df_jetsqcd["scalarmass"]
df_jetsqcd["beta_res"] = (df_jetsqcd["beta"] - df_jetsqcd["scalarbeta"])/df_jetsqcd["scalarbeta"]
df_jets750["pt_res"] = (df_jets750["pt"] - df_jets750["scalarpt"])/df_jets750["scalarpt"]
df_jets750["mass_res"] = (df_jets750["mass"] - df_jets750["scalarmass"])/df_jets750["scalarmass"]
df_jets750["beta_res"] = (df_jets750["beta"] - df_jets750["scalarbeta"])/df_jets750["scalarbeta"]
df_jets300["pt_res"] = (df_jets300["pt"] - df_jets300["scalarpt"])/df_jets300["scalarpt"]
df_jets300["mass_res"] = (df_jets300["mass"] - df_jets300["scalarmass"])/df_jets300["scalarmass"]
df_jets300["beta_res"] = (df_jets300["beta"] - df_jets300["scalarbeta"])/df_jets300["scalarbeta"]
df_jets200["pt_res"] = (df_jets200["pt"] - df_jets200["scalarpt"])/df_jets200["scalarpt"]
df_jets200["mass_res"] = (df_jets200["mass"] - df_jets200["scalarmass"])/df_jets200["scalarmass"]
df_jets200["beta_res"] = (df_jets200["beta"] - df_jets200["scalarbeta"])/df_jets200["scalarbeta"]

df_jets1000["t21"] = df_jets1000["t2"]/ df_jets1000["t1"]
df_jets1000["t32"] = df_jets1000["t3"]/ df_jets1000["t2"]
df_jets200["t21"] = df_jets200["t2"]/ df_jets200["t1"]
df_jets200["t32"] = df_jets200["t3"]/ df_jets200["t2"]
df_jets300["t21"] = df_jets300["t2"]/ df_jets300["t1"]
df_jets300["t32"] = df_jets300["t3"]/ df_jets300["t2"]
df_jets400["t21"] = df_jets400["t2"]/ df_jets400["t1"]
df_jets400["t32"] = df_jets400["t3"]/ df_jets400["t2"]
df_jets750["t21"] = df_jets750["t2"]/ df_jets750["t1"]
df_jets750["t32"] = df_jets750["t3"]/ df_jets750["t2"]
df_jetsqcd["t21"] = df_jetsqcd["t2"]/ df_jetsqcd["t1"]
df_jetsqcd["t32"] = df_jetsqcd["t3"]/ df_jetsqcd["t2"]
gc.collect()
print("converting to root")
df_jets1000.to_root("signal_1000.root",key="myTree")
df_jets750.to_root("signal_750.root",key="myTree")
df_jets400.to_root("signal_400.root",key="myTree")
df_jets300.to_root("signal_300.root",key="myTree")
df_jets200.to_root("signal_200.root",key="myTree")
df_jetsqcd.to_root("qcd_1000.root",key="myTree")
