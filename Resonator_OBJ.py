import numpy as np

class Resonator_2DOF:
    def __init__(self,propS,force,freq_sweep):
        dim = propS.shape
        if dim == (5,3):
            self.mass = propS[0,0:2]
            self.spring_lin = propS[1,:]
            self.spring_cub = propS[2,:]
            self.spring_quad = propS[3,:]
            self.damp = propS[4,:]
            self.force = force
            self.freq_sweep = freq_sweep
            self.Amp_1F = []
            self.Amp_1B = []
            self.Amp_2F = []
            self.Amp_2B = []
        else:
            print('Check zise of props matrix ofr 2DOF')
    
    def add_Amp1F(self,ampData):
        dimA = ampData.shape
        nF = self.force.size
        nW = self.freq_sweep.size
        if dimA == (nF,nW):
            self.Amp_1F = ampData
        else:
            print('Wrong Amp_1F Size')
    
    def add_Amp1B(self,ampData):
        dimA = ampData.shape
        nF = self.force.size
        nW = self.freq_sweep.size
        if dimA == (nF,nW):
            self.Amp_1B = ampData
        else:
            print('Wrong Amp_1B Size')

    def add_Amp2F(self,ampData):
        dimA = ampData.shape
        nF = self.force.size
        nW = self.freq_sweep.size
        if dimA == (nF,nW):
            self.Amp_2F = ampData
        else:
            print('Wrong Amp_2F Size')
    
    def add_Amp2B(self,ampData):
        dimA = ampData.shape
        nF = self.force.size
        nW = self.freq_sweep.size
        if dimA == (nF,nW):
            self.Amp_2B = ampData
        else:
            print('Wrong Amp_2B Size')