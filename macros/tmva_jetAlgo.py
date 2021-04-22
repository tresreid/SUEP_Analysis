from ROOT import TMVA, TFile, TTree, TCut, TString

TMVA.Tools.Instance()

siginputFile = TFile.Open("signal_1000.root")
bkginputFile = TFile.Open("qcd_1000.root")
outputFile = TFile.Open("TMVA_output1000.root","RECREATE")

factory = TMVA.Factory("TMVAClassification",outputFile, "!V:!Silent:Color:!DrawProgressBar:AnalysisType=Classification")


#declare variables in data loader
loader = TMVA.DataLoader("dataset_dnn")
loader.AddVariable("girth","girth","units","F")
loader.AddVariable("mass","mass","units","F")
loader.AddVariable("pt","pt","units","F")
loader.AddVariable("medpt","medpt","units","F")
loader.AddVariable("trackpt","trackpt","units","F")
loader.AddVariable("nTracks","nTracks","units","F")

#loader.AddVariable("")
#loader.AddVariable("var5 := var1-var2")
#loader.AddVariable("var6 := var1+var2")

#setup datasets
tsignal = siginputFile.Get("myTree")
tbackground = bkginputFile.Get("myTree")
loader.AddSignalTree(tsignal)
loader.AddBackgroundTree(tbackground)
loader.SetSignalWeightExpression("wgt")
loader.SetBackgroundWeightExpression("wgt")
cut_s = TCut("Jet_algo == 1 && R == 1.5")
cut_b = TCut("Jet_algo == 1 && R == 1.5")
#loader.PrepareTrainingAndTestTree(TCut(""),"nTrain_Signal=1000:nTrain_Background=1000:SplitMode=Random:NormMode=NumEvents:!V")
loader.PrepareTrainingAndTestTree(cut_s,cut_b,"nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V")

#configure netork layout
layoutString = TString("Layout=TANH|128,TANH|128,LINEAR")

training0 = TString("LearningRate=1e-1,Momentum=0.9,Repetitions=1,"
                        "ConvergenceSteps=2,BatchSize=256,TestRepetitions=10,"
                        "WeightDecay=1e-4,Regularization=L2,"
                        "DropConfig=0.0+0.5+0.5+0.5, Multithreading=True")
training1 = TString("LearningRate=1e-2,Momentum=0.9,Repetitions=1,"
                        "ConvergenceSteps=2,BatchSize=256,TestRepetitions=10,"
                        "WeightDecay=1e-4,Regularization=L2,"
                        "DropConfig=0.0+0.0+0.0+0.0, Multithreading=True")
trainingStrategyString = TString("TrainingStrategy=")
trainingStrategyString += training0 + TString("|") + training1

dnnOptions = TString("!H:!V:ErrorStrategy=CROSSENTROPY:VarTransform=N:"
        "WeightInitialization=XAVIERUNIFORM")
dnnOptions.Append(":")
dnnOptions.Append(layoutString)
dnnOptions.Append(":")
dnnOptions.Append(trainingStrategyString)

#booking Methods
#stdOptions =  dnnOptions + ":Architecture=CPU"
#factory.BookMethod(loader, TMVA.Types.kDNN, "DNN", stdOptions)
factory.BookMethod(loader,TMVA.Types.kBDT, "BDT",
                   "!V:NTrees=200:MinNodeSize=2.5%:MaxDepth=2:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" )
#factory.BookMethod(loader, TMVA.Types.kMLP, "MLP",
#                   "!H:!V:NeuronType=tanh:VarTransform=N:NCycles=100:HiddenLayers=N+5:TestRate=5:!UseRegulator" )
factory.BookMethod(loader, TMVA.Types.kCuts, "CutsSA",
                    "!H:!V:FitMethod=SA:EffSel")
factory.BookMethod(loader, TMVA.Types.kCuts, "CutsGA",
                    "!H:!V:FitMethod=GA:EffSel")
#factory.BookMethod(loader, TMVA.Types.kCuts, "CutsPCA",
#                    "!H:!V:FitMethod=MC:EffSel:VarProp=FSmart:VarTransform=PCA")
#factory.BookMethod(loader, TMVA.Types.kCuts, "CutsD",
#                    "!H:!V:FitMethod=MC:EffSel:VarProp=FSmart:VarTransform=Decorrelate")


#train
factory.TrainAllMethods()

#test and eval
factory.TestAllMethods()
factory.EvaluateAllMethods()

#c = factory.GetROCCurve(loader)
#c.Draw()
