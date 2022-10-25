#from ROOT import TFile, TLegend, TCanvas, TPad, THStack, TF1, TPaveText, TGaxis, SetOwnership, TObject, gStyle,TH1F
from ROOT import *
import os

import sys
from sys import argv
from optparse import OptionParser

from numpy import log10
from array import array
import math

padRatio = 0.25
padOverlap = 0.15

padGap = 0.01

plotDirectory = "data_pre_plots"

_fileDir ="/nfs/dust/cms/user/hugobg/ZPrime_102X/analysis_output/2016_CHS/muon/"

gROOT.SetBatch(True)

YesLog = True
NoLog=False


# Histogram Information:
# [X-axis title, 
#  Y-axis title,
#  Rebinning factor,
#  [x-min,x-max], -1 means keep as is
#  Extra text about region
#  log plot]
regionText ="loose selection"



import CMS_lumi

from Style import *
thestyle = Style()

HasCMSStyle = False
style = None
if os.path.isfile('tdrstyle.C'):
    ROOT.gROOT.ProcessLine('.L tdrstyle.C')
    ROOT.setTDRStyle()
    print "Found tdrstyle.C file, using this style."
    HasCMSStyle = True
    if os.path.isfile('CMSTopStyle.cc'):
        gROOT.ProcessLine('.L CMSTopStyle.cc+')
        style = CMSTopStyle()
        style.setupICHEPv1()
        print "Found CMSTopStyle.cc file, use TOP style if requested in xml file."
if not HasCMSStyle:
    print "Using default style defined in cuy package."
    thestyle.SetStyle()

ROOT.gROOT.ForceStyle()

#stackList = { "TTbar":[kRed], "DYJets":[kGreen], "QCD":[kYellow],"WJets":[kBlue], "ST":[kOrange], "Diboson":[kTeal]}
#stackList = { "TTToSemiLeptonic_2016_2":[kRed]} 
#stackList = {"TTToHadronic_2016":[kRed-3], "TTToOthers":[kRed+2], "TTToSemiLeptonic_2016":[kRed],"DYJetsToLL_M-50_HT_2016":[kBlue], "QCD_HT_2016":[kTeal], "WJetsToLNu_2016":[kGreen], "ST_2016":[kYellow], "WW_WZ_ZZ_2016v3":[kOrange]}
stackList = { "TTToSemiLeptonic_2016":[kRed], "TTTOthers":[kRed+2], "Others":[kBlue], "WJetsToLNu_2016":[kGreen]}

print stackList
#print stackList[2]
#channel_=argv[5]
year_=argv[3]
_channelText = "#mu+jets %i"%(year)

CMS_lumi.channelText = _channelText
CMS_lumi.writeChannelText = True
CMS_lumi.writeExtraText = True
if year==2016:
	lumi_num=6
elif year==2017:
	lumi_num=5
elif year==2018:
	lumi_num=4

H = 600;
W = 800;


# references for T, B, L, R                                                                                                             
T = 0.08*H
B = 0.12*H
L = 0.12*W
R = 0.1*W


# SetOwnership(canvas, False)
# SetOwnership(canvasRatio, False)
# SetOwnership(pad1, False)
# SetOwnership(pad2, False)



legendHeightPer = 0.04
legList = stackList.keys() 
#legList.reverse()

legendStart = 0.79
legendEnd = 0.97-(R/W)

#legend = TLegend(2*legendStart - legendEnd, 1-T/H-0.01 - legendHeightPer*(len(legList)+1), legendEnd, 0.99-(T/H)-0.01)
legend = TLegend(2*legendStart - legendEnd , 0.99 - (T/H)/(1.-padRatio+padOverlap) - legendHeightPer/(1.-padRatio+padOverlap)*round((len(legList)+1)/2.), legendEnd, 0.99-(T/H)/(1.-padRatio+padOverlap))
legend.SetNColumns(2)


legendR = TLegend(2*legendStart - legendEnd , 0.99 - (T/H)/(1.-padRatio+padOverlap) - 0.60*legendHeightPer/(1.-padRatio+padOverlap)*round((len(legList)+1)/2.)-0.1-0.03, legendEnd, 0.99-(T/H)/(1.-padRatio+padOverlap)-0.03)
#print "x1:",2*legendStart - legendEnd
#print "y1:",0.99 - (T/H)/(1.-padRatio+padOverlap) - 0.60*legendHeightPer/(1.-padRatio+padOverlap)*round((len(legList)+1)/2.)-0.1
#print "x2:", legendEnd
#print "y2:",0.99-(T/H)/(1.-padRatio+padOverlap)
legendR.SetNColumns(1)
#sys.exit()

legendR.SetBorderSize(0)
legendR.SetFillColor(0)



legend.SetBorderSize(0)
legend.SetFillColor(0)

_file = {}




canvas = TCanvas('c1','c1',W,H)
canvas.SetFillColor(0)
canvas.SetBorderMode(0)
canvas.SetFrameFillStyle(0)
canvas.SetFrameBorderMode(0)
canvas.SetLeftMargin( L/W )
canvas.SetRightMargin( R/W )
canvas.SetTopMargin( T/H )
canvas.SetBottomMargin( B/H )
canvas.SetTickx(0)

canvasRatio = TCanvas('c1Ratio','c1Ratio',W,H)
canvasRatio.SetFillColor(0)
canvasRatio.SetBorderMode(0)
canvasRatio.SetFrameFillStyle(0)
canvasRatio.SetFrameBorderMode(0)
canvasRatio.SetLeftMargin( L/W )
canvasRatio.SetRightMargin( R/W )
canvasRatio.SetTopMargin( T/H )
canvasRatio.SetBottomMargin( B/H )
canvasRatio.SetTickx(0)
canvasRatio.SetTicky(0)
canvasRatio.Draw()
canvasRatio.cd()


pad1 = TPad("zxc_p1","zxc_p1",0,padRatio-padOverlap,1,1)
pad2 = TPad("qwe_p2","qwe_p2",0,0,1,padRatio+padOverlap)
pad1.SetLeftMargin( L/W )
pad1.SetRightMargin( R/W )
pad1.SetTopMargin( T/H/(1-padRatio+padOverlap) )
pad1.SetBottomMargin( (padOverlap+padGap)/(1-padRatio+padOverlap) )
pad2.SetLeftMargin( L/W )
pad2.SetRightMargin( R/W )
pad2.SetTopMargin( (padOverlap)/(padRatio+padOverlap) )
pad2.SetBottomMargin( B/H/(padRatio+padOverlap) )
pad1.SetFillColor(0)
pad1.SetBorderMode(0)
pad1.SetFrameFillStyle(0)
pad1.SetFrameBorderMode(0)
pad1.SetTickx(0)
pad1.SetTicky(0)

pad2.SetFillColor(0)
pad2.SetFillStyle(4000)
pad2.SetBorderMode(0)
pad2.SetFrameFillStyle(0)
pad2.SetFrameBorderMode(0)
pad2.SetTickx(0)
pad2.SetTicky(0)


canvasRatio.cd()
pad1.Draw()
pad2.Draw()


canvas.cd()


canvas.ResetDrawn()

stack = THStack("hs","stack")
stack_up = THStack("hs","stack_up")
#SetOwnership(stack,True)
#histName=argv[1]
sum_=0
tree_MC={}
tree_MC_up={}
hist={}
hist1_={}
hist1_up={} 
histo={}

histograms = {
              "rec_chi2" : ["#chi^{2}", "Events", 30, [0,30]],
              "MET_pt"   : ["missing E_{T} [GeV]", "Events", 10, [50,1500]],
              "lep1_pt":[ "p_{T}^{#lep} [GeV]", "Events",50, [50, 500]],
              "Mttbar": ["Mttbar", "Events", 20, [0, 2000]], 
              "Ak8_j1_pt": ["Ak8_j1_pt", "Events", 25, [400,1200]],
              "Ak8_j1_tau32": ["Ak8_j1_tau32", "Events", 20, [0,1]],
              "Ak8_j1_tau21": ["Ak8_j1_tau21", "Events", 20, [0,1]],
              "Ak8_j1_mSD": ["Ak8_j1_mSD", "Events", 35, [0,400]],
              "N_Ak8": ["N_Ak8", "Events", 5, [0,5]],
              "N_Ak4": ["N_Ak4", "Events", 10, [0,10]],
              "TH_M": ["Mass_{t_{h}}", "Events", 40, [120,220]],
              "TL_M": ["Mass_{t_{l}}", "Events", 40, [120,220]], 
              "pT_had" : [ "p_{T}^{t_{had}} [GeV]", "Events",50, [0, 1200]],
              "pT_lep" : [ "p_{T}^{t_{#lep}} [GeV]", "Events",50, [0, 1200]],
              "pT_ttbar" : [ "p_{T}^{t#bar{t}} [GeV]", "Events",50, [0, 1200]],
#              "TH_M": ["TH_M", "Events", 40, [20,300]],
#              "TL_M": ["TL_M", "Events", 40, [20,300]],
#              "W_tagged_jet":["W_tagged_jet","Events",25,[40,200]], 
              "DeltaY": ["DeltaY", "Events", 2, [-2,2]], 
#	      "phi_Ak8Puppijets": ["#phi^{AK8Puppi jets}", "Events", 35, [-3.5, 3.5]], 
#              "phi_jet": ["#phi^{AK4 jets}", "Events", 35, [-3.5, 3.5]],
#              "dphi_jet1_MET":["#Delta#phi(jet1, MET)""Events", 40, [0, 4.0]],
#              "dphi_mu_MET":  ["#Delta#phi(#mu, MET)""Events", 40, [0, 4.0]],            
#              "st"    : ["S_{T} [GeV]", "Events", 50, [0,5000]],
#              "deepjetbscore_jet": ["DeepJet b-tag score all AK4 jets", "Events",20, [0, 1]],
#              "deepjetbscore_jet1": ["DeppJet b-tag score AK4 jet 1}", "Events",20, [0, 1]],
#              "pt_jet1": ["p_{T}^{jet 1} [GeV]", "Events", 20, [100, 900]],
#              "pt_jet": ["p_{T}^{jet} [GeV]","Events", 20, [100, 900]],
#              "pt_mu":[ "p_{T}^{#mu} [GeV]", "Events",50, [50, 500]],
#              "reliso_mu":[ "#mu rel. Iso", "Events", 20, [0, 0.5]],
#              "dR_mu_jet":[ "#DeltaR_{min}(#mu, jet)","Events", 60,[ 0, 3]],
#              "ptrel_mu_jet":["p_{T}^{rel}(#mu1, jet)", "Events",50, [0, 500]],
}

#sample_names = ["QCD","ST","DYJets","WJets","TTbar"]
#sample_names = ["TTToSemiLeptonic_2016_2"]
#sample_names = ["DYJetsToLL_M-50_HT_2016", "QCD_HT_2016", "WJetsToLNu_2016", "ST_2016", "WW_WZ_ZZ_2016", "TTToSemiLeptonic_2016","TTTo2L2Nu_2016", "TTToHadronic_2016"]
sample_names = ["TTToSemiLeptonic_2016","WJetsToLNu_2016","DYJetsToLL_M-50_HT_2016", "QCD_HT_2016", "ST_2016"]
#once ttbar is added:
#sample_names = ["TTToSemiLeptonic_2016","TToOthers","WJetsToLNu_2016","DYJetsToLL_M-50_HT_2016", "QCD_HT_2016", "ST_2016"]
for histName in histograms:
	for sample in sample_names:
        	print sample, histName
		_file[sample] = TFile("%s/uhh2.AnalysisModuleRunner.MC.%s.root"%(_fileDir,sample),"read")
        	print "%s/uhh2.AnalysisModuleRunner.MC.%s.root"%(_fileDir,sample)
		tree_MC[sample]=_file[sample].Get("AnalysisTree")
                tree_MC_up[sample]=_file[sample].Get("AnalysisTree")   
 	        tree_MC[sample].Draw("%s>>h2_%s(%i,%i,%f)"%(histName,sample,histograms[histName][2],histograms[histName][3][0],histograms[histName][3][1]),"(ttagN <= %s && wtagN <= %s && btagN >= 1 && rec_chi2<30)*weight*(weight_sfmu_HighPtID)*(weight_pu)*(weight_sfmu_MuonTrigger)*(weight_toptagSF_)*(weight_pt_rew_nolimit)*(weight_btag)*(muonrecSF_nominal)*0.95"%(argv[1],argv[2]))
                hist1_[sample] = tree_MC[sample].GetHistogram()

                if(sample == "DYJetsToLL_M-50_HT_2016" or sample == "QCD_HT_2016" or sample == "ST_2016"):
                    hist1_[sample].SetFillColor(stackList["Others"][0])
                    hist1_[sample].SetLineColor(stackList["Others"][0])
                    if(sample == "ST_2016"):
                        legendR.AddEntry(hist1_[sample],"Other bkg",'f')

                else:
                    hist1_[sample].SetFillColor(stackList[sample][0])
                    hist1_[sample].SetLineColor(stackList[sample][0])
                    if(sample == "TTToSemiLeptonic_2016"):
                        legendR.AddEntry(hist1_[sample],"t#bar{t} (l+jets)",'f')
                    if(sample == "WJetsToLNu_2016"):
                        legendR.AddEntry(hist1_[sample],"W+jets",'f')
                    if(sample == "TTToOthers"):
                        legendR.AddEntry(hist1_[sample],"t#bar{t} (others)",'f')



                    hist1_[sample].SetYTitle(histograms[histName][1])
                stack.Add(hist1_[sample])




_file["Data"] = TFile("%s/uhh2.AnalysisModuleRunner.DATA.DATA_SingleMuon_Run2016.root"%(_fileDir),"read")
print "%s/uhh2.AnalysisModuleRunner.DATA.DATA.root"%(_fileDir)

tree = _file["Data"].Get("AnalysisTree")
tree.Draw("%s>>dat_hist(%i,%i,%f)"%(histName,histograms[histName][2],histograms[histName][3][0],histograms[histName][3][1]),"ttagN <= %s && wtagN <= %s && btagN >= 1 && rec_chi2<30"%(argv[2],argv[3]))
dataHist=tree.GetHistogram()
print "total:",dataHist.Integral()
print "bins:",dataHist.GetNbinsX()
data=dataHist.GetBinContent(1)
print data
dataHist.SetMarkerColor(kBlack)
dataHist.SetYTitle(histograms[histName][1])     
dataHist.Draw("pe,x0")
stack.Draw("HIST,SAME")

errorban=stack.GetStack().Last().Clone("errorban")
errorban.Sumw2()
errorban.SetLineColor(kGray+2)
errorban.SetFillColor(kGray+2)
errorban.SetFillStyle(3245)
errorban.SetMarkerSize(0)
errorban.Draw("E2,SAME")

oneLine = TF1("oneline","1",-9e9,9e9)
oneLine.SetLineColor(kBlack)
oneLine.SetLineWidth(1)
oneLine.SetLineStyle(2)
	

maxVal =stack.GetMaximum()
minVal = 1
minVal = max(stack.GetStack()[0].GetMinimum(),1)
print minVal, maxVal

print "data, top:",dataHist.GetBinContent(1)

print "max stack:",1.5*maxVal
log=0
if log:
	stack.SetMaximum(10**(1.5*log10(maxVal) - 0.5*log10(minVal)))
else:
	stack.SetMaximum(1.7*maxVal)
stack.SetMinimum(minVal)

errorband=stack.GetStack().Last().Clone("error")
errorband.Sumw2()
errorband.SetLineColor(kBlack)
errorband.SetFillColor(kBlack)
errorband.SetFillStyle(3245)
errorband.SetMarkerSize(0)


canvas.Clear()
canvasRatio.cd()
canvasRatio.ResetDrawn()
canvasRatio.Draw()
canvasRatio.cd()

pad1.Draw()
pad2.Draw()

pad1.cd()

pad1.SetLogy(log)

y2 = pad1.GetY2()

stack.Draw("HIST")


stack.GetXaxis().SetTitle('')
stack.GetYaxis().SetTitle(dataHist.GetYaxis().GetTitle())

stack.SetTitle('')
stack.GetXaxis().SetLabelSize(0)
stack.GetYaxis().SetLabelSize(gStyle.GetLabelSize()/(1.-padRatio+padOverlap))
stack.GetYaxis().SetTitleSize(gStyle.GetTitleSize()/(1.-padRatio+padOverlap))
stack.GetYaxis().SetTitleOffset(gStyle.GetTitleYOffset()*(1.-padRatio+padOverlap))
stack.GetYaxis().SetMaxDigits(4)
stack.GetYaxis().SetTitle("Events")

dataHist.Draw("E,X0,SAME")

errorban.Draw("E2,SAME")

legendR.AddEntry(dataHist, "Data", 'pe')

ratio = dataHist.Clone("temp")
temp = stack.GetStack().Last().Clone("temp")

for i_bin in range(1,temp.GetNbinsX()+1):
       	temp.SetBinError(i_bin,0.)

ratio = dataHist.Clone("temp")
ratio.Divide(temp)
print ratio.GetNbinsX()
print ratio
print ratio.GetNbinsX()
ratio.SetTitle('')
histName=argv[1]
ratio.GetXaxis().SetLabelSize(gStyle.GetLabelSize()/(padRatio+padOverlap))
ratio.GetYaxis().SetLabelSize(gStyle.GetLabelSize()/(padRatio+padOverlap))
ratio.GetXaxis().SetTitleSize(gStyle.GetTitleSize()/(padRatio+padOverlap))
ratio.GetYaxis().SetTitleSize(gStyle.GetTitleSize()/(padRatio+padOverlap))
ratio.GetYaxis().SetTitleOffset(gStyle.GetTitleYOffset()*(padRatio+padOverlap-padGap))
ratio.GetYaxis().SetRangeUser(0.4,1.6)
ratio.GetYaxis().SetNdivisions(504)
ratio.GetXaxis().SetTitle(histograms[histName][0])
ratio.GetYaxis().SetTitle("Data/MC")
CMS_lumi.CMS_lumi(pad1, lum_num, 11)
legendR.SetTextSize(0.05)
legendR.Draw()
pad2.cd()
ratio.SetMarkerStyle(dataHist.GetMarkerStyle())
ratio.SetMarkerSize(dataHist.GetMarkerSize())
ratio.SetLineColor(dataHist.GetLineColor())
ratio.SetLineWidth(dataHist.GetLineWidth())
ratio.Draw('e,x0')
errorband.Divide(temp)



for i in range(1, errorband.GetNbinsX() + 1):
    if(errorban.GetBinContent(i) == 0):
        errorband.SetBinError(i,0)
    else:
        errorband.SetBinError(i, errorban.GetBinError(i)/errorban.GetBinContent(i))

errorband.Draw('e2,same')


oneLine.Draw("same")
	
canvasRatio.Update()
canvasRatio.RedrawAxis()

histName=argv[1]

if log:
	canvasRatio.SaveAs("%s_new_log.pdf"%(histName))
else:
	canvasRatio.SaveAs("%s_new.pdf"%(histName))




