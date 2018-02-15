#!/usr/bin/env python

from math import *

depth = 200. #um
bias = 80. #V
T=263.15+25 #K #Taken from PixelDistortionsSimpleSurfaceChargesGenerator.cxx
bField = 2.# Tesla = V*s/m^2
betaElectrons = 3.0E-16 #cm2/ns; 

#These parameterizations come from C. Jacoboni et al., Solid-State Electronics 20 (1977) 77-89. (see also https://cds.cern.ch/record/684187/files/indet-2001-004.pdf).

#Electrons
vs_e = 15.3*T**-0.87 #mm/ns
Ec_e = 1.01*10**-7*T**1.55 #MV/mm
beta_e = 2.57*10**-2*T**0.66
r_e = 1.13+0.0008*(T-273) #Hall factor

#Holes
vs_h = 1.62*T**-0.52 #mm/ns
Ec_h = 1.24*10**-7*T**1.68 #MV/mm
beta_h = 0.46*T**0.17
r_h = 0.72 - 0.0005*(T-273)

def mu(E,hole):
    #Returns the mobility in mm^2/(MV*ns)  
    #expects E in MV/mm
    Ec_h = 1.24*10**-7*T**1.68 
    if hole:
        num = vs_h/Ec_h
        den = (1+(E/Ec_h)**beta_h)**(1/beta_h)
        return r_h*num/den
    else:
        num = vs_e/Ec_e
        den = (1.+(E/Ec_e)**beta_e)**(1./beta_e)
        return r_e*num/den
    pass

def Ramo(z,a):
    L = depth
    den = 2.-exp(-a)-exp(-1.)
    num = exp(-a*z/L)+exp(-z/L)-exp(-a)-exp(-1.)
    return num/den

def Qe(z,fluence,a):
    myE = (bias/depth)*0.001 #MV/mm
    mymu = mu(myE,False)
    myv = r_e*mymu*myE*1000. #in microns/ns
    tau = 1./(fluence*betaElectrons) #in ns
    zbar = myv*tau
    lam = myv*tau/depth
    myexp = exp(-z/zbar)
    outval = myexp-Ramo(z,a)-(myexp/(2.-exp(-a)-exp(-1.)))*((1.-exp((1.-a*lam)*z/zbar))/(1.-a*lam) + (1.-exp((1.-lam)*z/zbar))/(1.-lam) - (exp(-a)+exp(-1.))*(1-myexp))
    return outval

def Qh(z,fluence,a):
    myE = (bias/depth)*0.001 #MV/mm
    mymu = mu(myE,True)
    myv = r_h*mymu*myE*1000. #in microns/ns
    tau = 1./(fluence*betaElectrons) #in ns
    zhat = myv*tau
    lamh = myv*tau/depth
    myexp = exp(z/zhat)
    outval = Ramo(z,a)-(myexp/(2.-exp(-a)-exp(-1.)))*((exp(-(1.+a*lamh)*z/zhat)-exp(-(1./lamh+a)))/(1.+a*lamh) + (exp(-(1.+lamh)*z/zhat)-exp(-(1./lamh+1)))/(1.+lamh) - (exp(-a)+exp(-1.))*(exp(-z/zhat)-exp(-depth/zhat)))
    return outval

def Integrated(fluence,a):
    myE = (bias/depth)*0.001 #MV/mm
    mymu = mu(myE,False)
    myv = r_e*mymu*myE*1000. #in microns/ns
    tau = 1./(fluence*betaElectrons) #in ns
    zbar = myv*tau
    lam = myv*tau/depth

    mymu = mu(myE,True)
    myv = r_h*mymu*myE*1000. #in microns/ns
    zhat = myv*tau
    lamh = myv*tau/depth

    x = (2.-exp(-a)-exp(-1.))**(-1)
    A = 1.-x/(1.-a*lam)-x/(1.-lam)+x*(exp(-a)+exp(-1.))
    B = x/(1.-a*lam)
    C = x/(1.-lam)
    D = x*exp(-(1./lamh+a))/(1.+a*lamh)+x*exp(-(1./lamh+1.))/(1.+lamh)-x*(exp(-a)+exp(-1.))*exp(-depth/zhat)
    E = -x/(1.+a*lamh)
    F = -x/(1.+lamh)
    G = -x*(exp(-a)+exp(-1.))

    return ((A*zbar)*(1.-exp(-depth/zbar)) + (B*zbar / (a*lam))*(1.-exp(-a*lam*depth/zbar)) + (C*zbar / (lam))*(1.-exp(-lam*depth/zbar)) + (D*zhat)*(exp(depth/zhat)-1.) + (E*zhat / (a*lamh))*(1.-exp(-a*lamh*depth/zhat)) + (F*zhat / (lamh))*(1.-exp(-lamh*depth/zhat)) + 0.5*G*(zbar*(1.-exp(-2.*depth/zbar))-2*depth)) / depth
