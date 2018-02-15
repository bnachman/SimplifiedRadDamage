#!/usr/bin/env python

import sys

from ROOT import *
import numpy as np
from math import *
from RadHelpers import *

def myText(x,y,text, color = 1):
    l = TLatex()
    #l.SetTextSize(0.04)
    l.SetTextSize(0.025)
    l.SetNDC()
    l.SetTextColor(color)
    l.DrawLatex(x,y,text)
    pass

#Now, for the main program.
gStyle.SetOptStat(0)
gStyle.SetPadLeftMargin(0.2)

ChargeFrac_e = TH1F("","",200,0,200)
ChargeFrac_h = TH1F("","",200,0,200)
ChargeFrac = TH1F("","",200,0,200)

flvals = []
evals = []
hvals = []
vals = []

for fluence in [1e12,5e12,1e13,5e13,1e14,5e14,1e15,2e15,3e15,4e15,5e15,6e15,7e15,8e15,9e15,1e16]:
    for i in range(1,ChargeFrac.GetNbinsX()+1):
        z0 = ChargeFrac.GetXaxis().GetBinCenter(i) #in microns
        ChargeFrac_e.SetBinContent(i,Qe(z0,fluence,12))
        ChargeFrac_h.SetBinContent(i,Qh(z0,fluence,12))
        ChargeFrac.SetBinContent(i,Qe(z0,fluence,12)+Qh(z0,fluence,12))
        pass
    flvals += [fluence]
    evals += [ChargeFrac_e.Integral()/200.]
    hvals += [ChargeFrac_h.Integral()/200.]
    #vals += [ChargeFrac.Integral()/200.]
    vals += [Integrated(fluence,12)]
    pass

graph_e = TGraph(len(evals),np.array(flvals),np.array(evals))
graph_h = TGraph(len(evals),np.array(flvals),np.array(hvals))
graph_all = TGraph(len(evals),np.array(flvals),np.array(vals))

graph_e.SetLineColor(2)
graph_e.SetLineStyle(3)
graph_h.SetLineColor(4)
graph_h.SetLineStyle(7)
graph_all.SetLineColor(1)

c = TCanvas("a","a",650,500)
mydummy = TH1F("","",5,1e12,1e16)
gPad.SetLogx()
mydummy.GetYaxis().SetRangeUser(0,1.2)
mydummy.GetXaxis().SetTitle("Fluence [n_{eq}/cm^{2}]")
mydummy.GetYaxis().SetTitleOffset(2.2)
mydummy.GetXaxis().SetTitleOffset(1.2)
mydummy.GetYaxis().SetTitle("Induced Charge / electron charge")
mydummy.GetXaxis().SetNdivisions(505)
mydummy.GetYaxis().SetNdivisions(505)
mydummy.SetTitle("")
mydummy.Draw()

leg = TLegend(.6,.78,0.9,.9);
leg.SetTextFont(42);
leg.SetFillStyle(0);
leg.SetFillColor(0);
leg.SetBorderSize(0);
leg.SetMargin(0.1);
leg.AddEntry(graph_all,"e+h","L")
leg.AddEntry(graph_e,"e","L")
leg.AddEntry(graph_h,"h","L")
leg.Draw();

myText(0.19,0.92,"#it{#bf{#scale[1.2]{Simplified Radiation Damage Model}}}")
myText(0.25,0.85,"#bf{#scale[1.2]{IBL, 80V}}");

graph_e.Draw("samec")
graph_h.Draw("samec")
graph_all.Draw("samec")

#c.Update()
c.Print("chargefrac_integrated.pdf")


