import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

firstfile = open('jet_comparisons_10000v3events.txt')
next(firstfile)
jets = []
for line in firstfile.readlines():
  cols = line.rstrip().split(' ')
  jets.append({"Event":int(cols[0]), "Jet_algo": cols[1], "R":float(cols[2]), "ntracks_lr": int(cols[3]), "ntracks_lg": int(cols[4]), "ntracks_sr": int(cols[5]), "ntracks_sg": int(cols[6]), "ntracks_mr": int(cols[7]), "ntracks_mg": int(cols[8]), "ntracks_scalar": int(cols[9]), "pt_lr": float(cols[10]), "pt_lg": float(cols[11]), "pt_sr": float(cols[12]), "pt_sg": float(cols[13]), "pt_mr": float(cols[14]), "pt_mg": float(cols[15]), "pt_scalar": float(cols[16]), "dR_l": float(cols[17]), "dR_m": float(cols[18]), "scalar_part_lg": int(cols[19]), "isr_part_lg": int(cols[20]), "scalar_part_lr": int(cols[21]), "isr_part_lr": int(cols[22]), "scalar_part_sg": int(cols[23]), "isr_part_sg": int(cols[24]), "scalar_part_sr": int(cols[25]), "isr_part_sr": int(cols[26]), "scalar_part_mg": int(cols[27]), "isr_part_mg": int(cols[28]), "scalar_part_mr": int(cols[29]), "isr_part_mr": int(cols[30])})

#secondfile = open('jet_comparisons_10000v2events.txt')
#next(secondfile)
#for line in secondfile.readlines():
#  cols = line.rstrip().split(' ')
#  jets.append({"Event":int(cols[0]), "Jet_algo": cols[1], "R":float(cols[2]), "ntracks_lr": int(cols[3]), "ntracks_lg": int(cols[4]), "ntracks_sr": int(cols[5]), "ntracks_sg": int(cols[6]), "ntracks_mr": int(cols[7]), "ntracks_mg": int(cols[8]), "ntracks_scalar": int(cols[9]), "pt_lr": float(cols[10]), "pt_lg": float(cols[11]), "pt_sr": float(cols[12]), "pt_sg": float(cols[13]), "pt_mr": float(cols[14]), "pt_mg": float(cols[15]), "pt_scalar": float(cols[16]), "dR_l": float(cols[17]), "dR_m": float(cols[18]), "scalar_part_lg": int(cols[19]), "isr_part_lg": int(cols[20]), "scalar_part_lr": int(cols[21]), "isr_part_lr": int(cols[22]), "scalar_part_sg": int(cols[23]), "isr_part_sg": int(cols[24]), "scalar_part_sr": int(cols[25]), "isr_part_sr": int(cols[26]), "scalar_part_mg": int(cols[27]), "isr_part_mg": int(cols[28]), "scalar_part_mr": int(cols[29]), "isr_part_mr": int(cols[30])})

df_jets = pd.DataFrame(jets)
df_jets["scalar_frac_lr"] = df_jets["scalar_part_lr"]/(df_jets["scalar_part_lr"]+df_jets["isr_part_lr"])
df_jets["scalar_frac_sr"] = df_jets["scalar_part_sr"]/(df_jets["scalar_part_sr"]+df_jets["isr_part_sr"])
df_jets["scalar_frac_mr"] = df_jets["scalar_part_mr"]/(df_jets["scalar_part_mr"]+df_jets["isr_part_mr"])

df_jets["scalar_containment_lr"] = df_jets["scalar_part_lr"]/df_jets["ntracks_scalar"]
df_jets["scalar_containment_sr"] = df_jets["scalar_part_sr"]/df_jets["ntracks_scalar"]
df_jets["scalar_containment_mr"] = df_jets["scalar_part_mr"]/df_jets["ntracks_scalar"]

df_jets["pt_ratio_lr"] = df_jets["pt_lr"]/df_jets["pt_scalar"]
df_jets["pt_ratio_sr"] = df_jets["pt_sr"]/df_jets["pt_scalar"]
df_jets["pt_ratio_mr"] = df_jets["pt_mr"]/df_jets["pt_scalar"]

#print(df_jets)
akt_8 = df_jets[(df_jets["Jet_algo"] == "AKT") & (df_jets["R"] == 0.8)]
akt_10 = df_jets[(df_jets["Jet_algo"] == "AKT") & (df_jets["R"] == 1.0)]
akt_15 = df_jets[(df_jets["Jet_algo"] == "AKT") & (df_jets["R"] == 1.5)]
akt_20 = df_jets[(df_jets["Jet_algo"] == "AKT") & (df_jets["R"] == 2.0)]
ca_8 = df_jets[(df_jets["Jet_algo"] == "CA") & (df_jets["R"] == 0.8)]
ca_10 = df_jets[(df_jets["Jet_algo"] == "CA") & (df_jets["R"] == 1.0)]
ca_15 = df_jets[(df_jets["Jet_algo"] == "CA") & (df_jets["R"] == 1.5)]
ca_20 = df_jets[(df_jets["Jet_algo"] == "CA") & (df_jets["R"] == 2.0)]
kt_8 = df_jets[(df_jets["Jet_algo"] == "KT") & (df_jets["R"] == 0.8)]
kt_10 = df_jets[(df_jets["Jet_algo"] == "KT") & (df_jets["R"] == 1.0)]
kt_15 = df_jets[(df_jets["Jet_algo"] == "KT") & (df_jets["R"] == 1.5)]
kt_20 = df_jets[(df_jets["Jet_algo"] == "KT") & (df_jets["R"] == 2.0)]
#print(akt_8)

range_norm = [0.0,1.0]
range1 = [0.0,2.0]

def make_plot_akt(var,title,xtitle):
  fig, ax = plt.subplots()
  ax.set_title(title)
  plt.xlabel(xtitle)
  plt.ylabel("A.U.")

  if("pt" in title):
    rangex = range1
  else:
    rangex = range_norm

  ax.hist( akt_8[var],range=rangex,weights=np.ones_like(akt_8[var])/len(akt_8[var]), bins=20,histtype=u'step', label="AKT8")
  ax.hist(akt_10[var],range=rangex,weights=np.ones_like(akt_10[var])/len(akt_10[var]), bins=20,histtype=u'step',label="AKT10")
  ax.hist(akt_15[var],range=rangex,weights=np.ones_like(akt_15[var])/len(akt_15[var]), bins=20,histtype=u'step',label="AKT15")
  ax.hist(akt_20[var],range=rangex,weights=np.ones_like(akt_20[var])/len(akt_20[var]), bins=20,histtype=u'step',label="AKT20")
  ax.legend(loc ='upper left')
  fig.savefig('jet_output/%s_AKT.png'%var)
  plt.close()
  #plt.show()
def make_plot_kt(var,title,xtitle):
  fig, ax = plt.subplots()
  ax.set_title(title)
  plt.xlabel(xtitle)
  plt.ylabel("A.U.")

  if("pt" in title):
    rangex = range1
  else:
    rangex = range_norm

  ax.hist( kt_8[var],range=rangex,weights=np.ones_like(kt_8[var])/len(kt_8[var]), bins=20,histtype=u'step', label="KT8")
  ax.hist(kt_10[var],range=rangex,weights=np.ones_like(kt_10[var])/len(kt_10[var]), bins=20,histtype=u'step',label="KT10")
  ax.hist(kt_15[var],range=rangex,weights=np.ones_like(kt_15[var])/len(kt_15[var]), bins=20,histtype=u'step',label="KT15")
  ax.hist(kt_20[var],range=rangex,weights=np.ones_like(kt_20[var])/len(kt_20[var]), bins=20,histtype=u'step',label="KT20")
  ax.legend(loc ='upper left')
  fig.savefig('jet_output/%s_KT.png'%var)
  plt.close()
  #plt.show()
def make_plot_ca(var,title,xtitle):
  fig, ax = plt.subplots()
  ax.set_title(title)
  plt.xlabel(xtitle)
  plt.ylabel("A.U.")

  if("pt" in title):
    rangex = range1
  else:
    rangex = range_norm

  ax.hist( ca_8[var],range=rangex,weights=np.ones_like(ca_8[var])/len(ca_8[var]), bins=20,histtype=u'step', label="CA8")
  ax.hist(ca_10[var],range=rangex,weights=np.ones_like(ca_10[var])/len(ca_10[var]), bins=20,histtype=u'step',label="CA10")
  ax.hist(ca_15[var],range=rangex,weights=np.ones_like(ca_15[var])/len(ca_15[var]), bins=20,histtype=u'step',label="CA15")
  ax.hist(ca_20[var],range=rangex,weights=np.ones_like(ca_20[var])/len(ca_20[var]), bins=20,histtype=u'step',label="CA20")
  ax.legend(loc ='upper left')
  fig.savefig('jet_output/%s_CA.png'%var)
  plt.close()
  #plt.show()
def make_plot_lowR(var,title,xtitle):
  fig, ax = plt.subplots()
  ax.set_title(title)
  plt.xlabel(xtitle)
  plt.ylabel("A.U.")

  if("pt" in title):
    rangex = range1
  else:
    rangex = range_norm

  ax.hist( ca_8[var], range=rangex,weights=np.ones_like(ca_8[var])/len(ca_8[var]), bins=20,histtype=u'step', label="CA8")
  ax.hist(akt_8[var], range=rangex,weights=np.ones_like(akt_8[var])/len(akt_8[var]), bins=20,histtype=u'step',label="AKT8")
  ax.hist(kt_8[var],  range=rangex,weights=np.ones_like(kt_8[var])/len(kt_8[var]), bins=20,histtype=u'step',label="KT8")
  ax.hist( ca_10[var],range=rangex,weights=np.ones_like(ca_10[var])/len(ca_10[var]), bins=20,histtype=u'step', label="CA10")
  ax.hist(akt_10[var],range=rangex,weights=np.ones_like(akt_10[var])/len(akt_10[var]), bins=20,histtype=u'step',label="AKT10")
  ax.hist(kt_10[var], range=rangex,weights=np.ones_like(kt_10[var])/len(kt_10[var]), bins=20,histtype=u'step',label="KT10")
  ax.legend(loc ='upper left')
  fig.savefig('jet_output/%s_lowR.png'%var)
  plt.close()
  #plt.show()
def make_plot_highR(var,title,xtitle):
  fig, ax = plt.subplots()
  ax.set_title(title)
  plt.xlabel(xtitle)
  plt.ylabel("A.U.")

  if("pt" in title):
    rangex = range1
  else:
    rangex = range_norm

  ax.hist( ca_15[var],range=rangex,weights=np.ones_like(ca_15[var])/len(ca_15[var]), bins=20,histtype=u'step', label="CA15")
  ax.hist(akt_15[var],range=rangex,weights=np.ones_like(akt_15[var])/len(akt_15[var]), bins=20,histtype=u'step',label="AKT15")
  ax.hist(kt_15[var], range=rangex,weights=np.ones_like(kt_15[var])/len(kt_15[var]), bins=20,histtype=u'step',label="KT15")
  ax.hist( ca_20[var],range=rangex,weights=np.ones_like(ca_20[var])/len(ca_20[var]), bins=20,histtype=u'step', label="CA20")
  ax.hist(akt_20[var],range=rangex,weights=np.ones_like(akt_20[var])/len(akt_20[var]), bins=20,histtype=u'step',label="AKT20")
  ax.hist(kt_20[var], range=rangex,weights=np.ones_like(kt_20[var])/len(kt_20[var]), bins=20,histtype=u'step',label="KT20")
  ax.legend(loc ='upper left')
  fig.savefig('jet_output/%s_highR.png'%var)
  plt.close()
  #plt.show()

#make_plot_akt("scalar_frac_lr","constituent fraction from scalar daughters (lead reco)","scalar daughter constituents fraction")
#make_plot_kt("scalar_frac_lr","constituent fraction from scalar daughters (lead reco)","scalar daughter constituents fraction")
#make_plot_ca("scalar_frac_lr","constituent fraction from scalar daughters (lead reco)","scalar daughter constituents fraction")
#make_plot_lowR("scalar_frac_lr","constituent fraction from scalar daughters (lead reco)","scalar daughter constituents fraction")
#make_plot_highR("scalar_frac_lr","constituent fraction from scalar daughters (lead reco)","scalar daughter constituents fraction")
#make_plot_akt("scalar_containment_lr","fraction of total scalar daughters (lead reco)","scalar daughters in event fraction")
#make_plot_kt("scalar_containment_lr","fraction of total scalar daughters (lead reco)","scalar daughters in event fraction")
#make_plot_ca("scalar_containment_lr","fraction of total scalar daughters (lead reco)","scalar daughters in event fraction")
#make_plot_lowR("scalar_containment_lr","fraction of total scalar daughters (lead reco)","scalar daughters in event fraction")
#make_plot_highR("scalar_containment_lr","fraction of total scalar daughters (lead reco)","scalar daughters in event fraction")
#make_plot_akt("pt_ratio_lr","pt/ scalar pt (lead reco)","pt/scalar pt")
#make_plot_kt("pt_ratio_lr","pt/ scalar pt (lead reco)","pt/scalar pt")
#make_plot_ca("pt_ratio_lr","pt/ scalar pt (lead reco)","pt/scalar pt")
#make_plot_highR("pt_ratio_lr","pt/ scalar pt (lead reco)","pt/scalar pt")
#make_plot_lowR("pt_ratio_lr","pt/ scalar pt (lead reco)","pt/scalar pt")
#
#
#make_plot_akt("scalar_frac_mr","constituent fraction from scalar daughters (multi reco)","scalar daughter constituents fraction")
#make_plot_kt("scalar_frac_mr","constituent fraction from scalar daughters (multi reco)","scalar daughter constituents fraction")
#make_plot_ca("scalar_frac_mr","constituent fraction from scalar daughters (multi reco)","scalar daughter constituents fraction")
#make_plot_lowR("scalar_frac_mr","constituent fraction from scalar daughters (multi reco)","scalar daughter constituents fraction")
#make_plot_highR("scalar_frac_mr","constituent fraction from scalar daughters (multi reco)","scalar daughter constituents fraction")
#make_plot_akt("scalar_containment_mr","fraction of total scalar daughters (multi reco)","scalar daughters in event fraction")
#make_plot_kt("scalar_containment_mr","fraction of total scalar daughters (multi reco)","scalar daughters in event fraction")
#make_plot_ca("scalar_containment_mr","fraction of total scalar daughters (multi reco)","scalar daughters in event fraction")
#make_plot_lowR("scalar_containment_mr","fraction of total scalar daughters (multi reco)","scalar daughters in event fraction")
#make_plot_highR("scalar_containment_mr","fraction of total scalar daughters (multi reco)","scalar daughters in event fraction")
#make_plot_akt("pt_ratio_mr","pt/ scalar pt (multi reco)","pt/scalar pt")
#make_plot_kt("pt_ratio_mr","pt/ scalar pt (multi reco)","pt/scalar pt")
#make_plot_ca("pt_ratio_mr","pt/ scalar pt (multi reco)","pt/scalar pt")
#make_plot_highR("pt_ratio_mr","pt/ scalar pt (multi reco)","pt/scalar pt")
#make_plot_lowR("pt_ratio_mr","pt/ scalar pt (multi reco)","pt/scalar pt")
#
#
#
#make_plot_akt("dR_l","dR (lead reco)","dR")
#make_plot_kt("dR_l","dR (lead reco)","dR")
#make_plot_ca("dR_l","dR (lead reco)","dR")
#make_plot_highR("dR_l","dR (lead reco)","dR")
#make_plot_lowR("dR_l","dR (lead reco)","dR")
#make_plot_akt("dR_m","dR (multi reco)","dR")
#make_plot_kt("dR_m","dR (multi reco)","dR")
#make_plot_ca("dR_m","dR (multi reco)","dR")
#make_plot_highR("dR_m","dR (multi reco)","dR")
#make_plot_lowR("dR_m","dR (multi reco)","dR")

def get_score(df,tit):
  score = []
  score.append(len(df[df["scalar_frac_lr"] >= 0.8])/ len(df["scalar_frac_lr"]))
  score.append(len(df[df["scalar_frac_mr"] >= 0.8])/ len(df["scalar_frac_mr"]))
  score.append(len(df[df["scalar_containment_lr"] >= 0.8])/ len(df["scalar_containment_lr"]))
  score.append(len(df[df["scalar_containment_mr"] >= 0.8])/ len(df["scalar_containment_mr"]))
  score.append(len(df[(df["pt_ratio_lr"] >= 0.75) & (df["pt_ratio_lr"] <= 1.25)])/ len(df["pt_ratio_lr"]))
  score.append(len(df[(df["pt_ratio_mr"] >= 0.75) & (df["pt_ratio_mr"] <= 1.25)])/ len(df["pt_ratio_mr"]))
  score.append(len(df[(df["dR_l"] <= 0.1)])/ len(df["dR_l"]))
  score.append(len(df[(df["dR_m"] <= 0.1)])/ len(df["dR_m"]))
  print("%s %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f"%(tit,score[0],score[1],score[2],score[3],score[4],score[5],score[6],score[7]))  

print("scalar_frac_lr>= 0.8 scalar_frac_mr>=0.8 scalar_containment_lr>=0.8 scalar_containment_mr>=0.8 pt_ratio_lr(0.75,1.25) pr_ratio_mr(0.75,1.25) dR_l<=0.1 dR_m<=0.1")
get_score(akt_8,"akt8 ")
get_score(akt_10,"akt10")
get_score(akt_15,"akt15")
get_score(akt_20,"akt20")
get_score(ca_8,"ca8  ")
get_score(ca_10,"ca10 ")
get_score(ca_15,"ca15 ")
get_score(ca_20,"ca20 ")
get_score(kt_8,"kt8  ")
get_score(kt_10,"kt10 ")
get_score(kt_15,"kt15 ")
get_score(kt_20,"kt20 ")
