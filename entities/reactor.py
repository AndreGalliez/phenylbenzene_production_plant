import numpy as np
from scipy.optimize import fsolve, newton
        
class GasPhaseReactor:        
    def __init__(self, Fin, Win, Vr, Kr, ReacCoefs, model, P, T):
            self.Fin=Fin
            self.Win=Win
            self.Vr=Vr
            self.Fout=None
            self.Wout=[None]*len(Win)
            self.Kr=Kr
            self.model = model
            self.ReacCoefs = ReacCoefs
            self.P=P
            self.T=T
            
    def formulate_equations(self,x_out):
        Wx=list()
        for i in range(len(self.Win)):
            Wx.append(x_out[i])
        self.Fout=x_out[-1]
        f=list()
        RateTotal = 0.0
        for i in range(len(Wx)):  ##Formulate Component Mass Balances for each Component i
            RateBySpecies=0.0
            for j in range(len(self.ReacCoefs)):
                if np.abs(self.ReacCoefs[j][i]) > 0.0:
                    ReacModel=ReactionModel(self.Kr[j],self.ReacCoefs[j],self.model)        ##Instantiante and calculate reaction rate for species i reaction j
                    RateBySpecies=RateBySpecies+ReacModel.get_reaction_rate(i,Wx,self.P,self.T)##*(np.abs(self.ReacCoefs[j][i])/np.abs(self.ReacCoefs[j][0]))                  ##Formulate reaction rate for species over all equations
            fx=(self.Fin*self.Win[i]/self.Fout)+RateBySpecies*(self.Vr/self.Fout)-Wx[i]  ##Formulate the equation
            f.append(fx)
            RateTotal = RateTotal + RateBySpecies
        fx=(self.Fin-self.Fout+self.Vr*RateTotal) ##Formulate Overall Mass Balance Equation
        f.append(fx)
        return tuple(f)
            
    def evaluate(self,x_in):
        aux = fsolve(self.formulate_equations,x_in)
        for i in range(len(self.Wout)):
            self.Wout[i]=aux[i]
        self.Fout=aux[-1]


class ReactionModel:
    def __init__(self, Kr, reactionCoefficients, model):
        self.Kr=Kr
        self.ReacCoefs=reactionCoefficients
        self.Model=model

    def get_reaction_rate(self, indx, Wx, P, T):
        Cx=list()
        for i in range(len(Wx)):
            Cx.append(Wx[i]*P)
        if self.Model == "DirectPowerlaw":
            ri = self.Kr
            for i in range(len(self.ReacCoefs)):
                if self.ReacCoefs[i] < 0:
                    ri = ri* (Cx[i]**(np.abs(self.ReacCoefs[i])))
            return ri 
        if self.Model == "TwoWaysPowerlaw":
            rdir = self.Kr[0]
            rinv = self.Kr[1]
            for i in range(len(self.ReacCoefs)):
                if self.ReacCoefs[i] < 0:
                    rdir = rdir* (Cx[i]**(np.abs(self.ReacCoefs[i])))
                elif self.ReacCoefs[i] > 0:
                    rinv = rinv* (Cx[i]**(np.abs(self.ReacCoefs[i])))
            ri = rdir -rinv
            return ri * self.ReacCoefs[indx]/np.abs(self.ReacCoefs[0])
    
class ReactionRateConstant:
    def __init__(self, Ko, E, T):
        self.Ko=Ko
        self.E=E
        self.T=T
        self.Kr=np.empty_like(Ko)
        self.R=1.9872
        
    def evaluate_K(self):
        for i in range(self.Kr.shape[0]):
            for j in range(self.Kr.shape[1]):
                self.Kr[i][j]=self.Ko[i][j]*np.exp(-self.E[i][j]/(self.T*self.R))