import numpy as np

class Splitter: ## Divide a corrente de entrada em duas de saída de acordo com a proporção F1out/F2out =Cs
    """
        Responsável por realizar os cálculos referentes ao separados de correntes.
        Argumentos:
            Fin (float): Lista com a vazão da corrente de entrada.
            Win (list(float)): Lista com a composição da corrente de entrada.
            Cs (float): Razão de reciclo
        Atributos:
            Fout (dict(float)): Vazões da corrente de purga e reciclo a serem calculadas.
            Wout (list(float): Lista com a composição das correntes de saída a serem calculadas.
        Métodos:
            evaluate()
                : Realiza os calculos e atualiza os valores de Fout e Wout.
    """

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
    """
        Responsável por realizar os cálculos referentes ao misturador de correntes.
        Argumentos:
            Fin (list(float)): Lista com as vazões das correntes de entrada.
            Win (list(list(float))): Lista com as composições das correntes de entrada.
        Atributos:
            Fout (float): Vazões da corrente de saída a ser calculada.
            Wout (list(float): Lista com as composições da corrente de saída a serem calculadas.
        Métodos:
            evaluate()
                : Realiza os calculos e atualiza os valores de Fout e Wout.
    """

    def __init__(self,Fin,Win):
        self.Fin =np.array(Fin) #Recebe uma lista normal e retorna um vetor linha np
        self.Win =np.array(Win).T #Recebe uma lista de listas(vetores) e transforma na matriz transposta np
        self.Fout = None 
        self.Wout= [None]*len(Win[0]) #Cria uma lista com o mesmo numero de entradas das listas composição de corrente
        
    def evaluate(self):
        self.Fout=np.sum(self.Fin)
        for i in range(len(self.Win[:,0])):
            self.Wout[i]=(self.Fin@self.Win[i,:])/self.Fout