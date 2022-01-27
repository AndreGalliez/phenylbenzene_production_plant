import numpy as np

class Splitter: ## Divide a corrente de entrada em duas de saída de acordo com a proporção F1out/F2out =Cs
    def __init__(self,Fin,Win,Cs):
        self.Fin =Fin
        self.Win =Win
        self.Fout = dict()
        self.Wout = list()
        self.Cs=Cs
        
    def evaluate(self):
        self.Fout['F_recycle']=self.Fin*self.Cs
        self.Fout['F_purge']=self.Fin-self.Fout['F_recycle']
        for wi in self.Win: self.Wout.append(wi)
   
class Mixer: ##Faz o cálculo de mistura de correntes e seus componentes
    def __init__(self,Fin,Win):
        self.Fin =np.array(Fin) #Recebe uma lista normal e retorna um vetor linha np
        self.Win =np.array(Win).T #Recebe uma lista de listas(vetores) e transforma na matriz transposta np
        self.Fout = None 
        self.Wout= [None]*len(Win[0]) #Cria uma lista com o mesmo numero de entradas das listas composição de corrente
        
    def evaluate(self):
        self.Fout=np.sum(self.Fin)
        for i in range(len(self.Win[:,0])):
            self.Wout[i]=(self.Fin@self.Win[i,:])/self.Fout