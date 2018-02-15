#!/usr/bin/env python

import sys

if (len(sys.argv) < 2):
    print "Please enter a fluence."
    print "Run like python Simplified.py 1e14"
    exit(1)

fluence = float(sys.argv[1])
    
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

for i in range(1,ChargeFrac.GetNbinsX()+1):
    z0 = ChargeFrac.GetXaxis().GetBinCenter(i) #in microns
    ChargeFrac_e.SetBinContent(i,Qe(z0,fluence,12))
    ChargeFrac_h.SetBinContent(i,Qh(z0,fluence,12))
    ChargeFrac.SetBinContent(i,Qe(z0,fluence,12)+Qh(z0,fluence,12))
    pass

c = TCanvas("a","a",650,500)
mydummy = TH1F("","",200,0,200)
mydummy.GetYaxis().SetRangeUser(0,1.2)
mydummy.GetXaxis().SetTitle("Pixel Depth in Z [#mum]")
mydummy.GetYaxis().SetTitleOffset(2.2)
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
leg.AddEntry(ChargeFrac,"e+h","L")
leg.AddEntry(ChargeFrac_e,"e","L")
leg.AddEntry(ChargeFrac_h,"h","L")
leg.Draw();

myText(0.19,0.92,"#it{#bf{#scale[1.2]{Simplified Radiation Damage Model}}}")
myText(0.25,0.85,"#bf{#scale[1.2]{IBL, 80V, #Phi = "+sys.argv[1]+" n_{eq}/cm^{2}}}");

ChargeFrac.SetLineColor(1)
ChargeFrac.Draw("same")
ChargeFrac_e.SetLineStyle(7)
ChargeFrac_e.SetLineColor(2)
ChargeFrac_e.Draw("same")
ChargeFrac_h.SetLineStyle(3)
ChargeFrac_h.SetLineColor(4)
ChargeFrac_h.Draw("same")

#c.Update()
c.Print("chargefrac_"+sys.argv[1]+".pdf")


