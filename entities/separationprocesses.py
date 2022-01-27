import numpy as np
from scipy.optimize import fsolve, newton

class Flash: ##Faz o cálculo de flash para quando as condições de equilibrio e composição na entrada são conhecidas
    def __init__(self, _type, Fin, z, p_sat, p_or_beta):
        if _type == 'PT':
            self.P = p_or_beta
            self.B = None
        elif _type == 'ZB':
            self.B = p_or_beta
            self.P = None
        self.Fin = Fin
        self.Z = z
        self.P_sat = p_sat
        self.K=list()
        self.V=None
        self.L=None
        self.X=[None]*len(z)
        self.Y=[None]*len(z)

    def evaluate_K(self):
        for Pi_sat in self.P_sat:
            self.K.append(Pi_sat/self.P)
               
    def formulate_equations_PT(self, x):
        f=0.0
        for i in range(len(self.Z)):
            f=f+self.Z[i]*self.K[i]/(1+x*(self.K[i]-1))
        return f-1

    def evaluate_flash_PT(self, x_in):
        self.evaluate_K()
        self.B = newton(self.formulate_equations_PT, x_in)
        for i in range(len(self.Z)):
            self.Y[i] = self.Z[i]*self.K[i]/(1+self.B*(self.K[i]-1))
        for i in range(len(self.Z)):
            self.X[i]= self.Y[i]/self.K[i]
        self.L=self.Fin*(1-self.B)
        self.V=self.Fin-self.L

    def formulate_equations_ZB(self, x):
        f=0.0
        for i in range(len(self.Z)):
            f=f+self.Z[i]*(self.P_sat[i]/x)/(1+self.B*((self.P_sat[i]/x)-1))
        return f-1

    def evaluate_flash_ZB(self, x_in):
        self.P = fsolve(self.formulate_equations_ZB, x_in,maxiter=500)
        self.evaluate_K()
        for i in range(len(self.Z)):
            self.Y[i] = self.Z[i]*self.K[i]/(1+self.B*(self.K[i]-1))
        for i in range(len(self.Z)):
            self.X[i]= self.Y[i]/self.K[i]
        self.L=self.Fin*(1-self.B)
        self.V=self.Fin-self.L
    
   
class LiquidVaporEquilibriumConstant:
    def __init__(self, model, model_input):
        self.model = model
        self.P_sat = None
        if model == 'Raoult and Antoine':
            self.T=model_input[0]
            self.elv_coefficients=model_input[1]
        elif model == 'Raoult and Other':
            self.T=model_input[0]
            self.elv_coefficients=model_input[1]
            

    ##O metodo implementado usa Pressao em psia e temperatura em farenheit
    ##A equação de Antoine utilizada usa a Pressão crítica do componente como pressão de referência
    @staticmethod
    def antoine_method(T,singleComponentCoefficients):
        pi_sat = singleComponentCoefficients[3]*(np.exp(singleComponentCoefficients[0]-singleComponentCoefficients[1]/(T+singleComponentCoefficients[2])))
        return pi_sat
    
    ##O metodo implementado usa Pressao em Pascal e temperatura em Kelvin
    @staticmethod 
    def forgot_name_method(T,singleComponentCoefficients):
        pi_sat= 1000*(np.exp(singleComponentCoefficients[0]*np.log(T)+(singleComponentCoefficients[1]/T) +singleComponentCoefficients[2]+singleComponentCoefficients[3]*(T**2)))
        return pi_sat

    def calc_psats(self):
        P_sat=list()
        if self.model == 'Raoult and Antoine':    
            for singleComponentCoefficients in self.elv_coefficients:
                P_sat.append(self.antoine_method(self.T,singleComponentCoefficients))
        elif self.model == 'Raoult and Other':
            for i in range(len(self.elv_coefficients)):
                singleComponentCoefficients = self.elv_coefficients[i]
                if i in [0,1,2]:
                    P_sat.append(6894.76*self.antoine_method(((self.T - 273.15) * (9/5) + 32),singleComponentCoefficients))
                else:
                    P_sat.append(self.forgot_name_method(self.T,singleComponentCoefficients))
        self.P_sat = P_sat
        return P_sat
        ##Print das Pressoes e Constantes de Equilibrio
        #print(np.array(P_sat)/(10**5))
        #print(self.K)