import numpy as np
from scipy.optimize import fsolve, newton
        
class GasPhaseReactor:      
    """
        Responsável por realizar os cálculos referentes ao reator.
        Argumentos:
            Fin (float): Vazão de entrada.
            Win (list(float)): Composições de entrada.
            Fin (float): Volume do reator.
            Kr (list(list((float))): Lista com os pares (k_direta,k_reversa) de constantes reacionais.
            P (float): Pressão no reator.
            P (float): Temperatura no reator.
            ReacCoefs (list(list(int))): Lista com os coeficientes reacionais para cada reação:
                                        [[R1 Coeficientes],[R2 Coeficientes],...]
        Atributos:
            Fout (float): Vazão da corrente de saída (a ser calculado).
            Wout (list(float)): Composições da corrente de saída (a ser calculado).
        Métodos:
            formulate_equations()
                : Formula o sistema de equações a ser resolvido.
            evaluate()
                : Resolve o sistema de equações e calcula os valores da vazão e composições de saída.
    """  
    def __init__(self, Fin, Win, Vr, Kr, ReacCoefs, P, T):
            self.Fin=Fin
            self.Win=Win
            self.Vr=Vr
            self.Kr=Kr
            self.ReacCoefs = ReacCoefs
            self.P=P
            self.T=T
            self.Fout=None
            self.Wout=[None]*len(Win)
            
    def formulate_equations(self,initial_guess):
        """
            Formula o sistema de equações a ser resolvido.
            Argumentos:
                initial_guess (list(float)): Chute inicial para as vazões e composições a serem calculados.
            Retorna:
                 f (list(float)): Valor do resíduo de cada equação do sistema. Será avaliado pelo solver do scipy .
        """
        W_initial=list()
        for i in range(len(self.Win)):
            W_initial.append(initial_guess[i])
        self.Fout=initial_guess[-1]
        f=list()
        RateTotal = 0.0
        for i in range(len(W_initial)):  
            reactionModel=ReactionModel(self.Kr,self.ReacCoefs,self.P,i,W_initial)
            RateBySpecies=reactionModel.get_global_reaction_rate()
            fx=(self.Fin*self.Win[i]/self.Fout)+RateBySpecies*(self.Vr/self.Fout)-W_initial[i]
            f.append(fx)
            RateTotal = RateTotal + RateBySpecies
        fx=(self.Fin-self.Fout+self.Vr*RateTotal)
        f.append(fx)
        print()
        return tuple(f)
            
    def evaluate(self,initial_guess):
        """
            Resolve o sistema de eqiações e calcula os valores da vazão e composições de saída.
            Chama a função formulate_equations, passa o valor dos resíduos para a função fsolve do scipy que resolve o sistema não-linear.
            Na sequência calcula os atributos da corrente de saida e atualiza os parâmetros.
            Argumentos:
                initial_guess (list(float)): Chute inicial para as vazões e composições a serem calculados.

        """
        aux = fsolve(self.formulate_equations,initial_guess)
        for i in range(len(self.Wout)):
            self.Wout[i]=aux[i]
        self.Fout=aux[-1]


class ReactionModel:
    """
        Responsável por calcular a taxa reacional (r_i) para um determinado componente i.
        Argumentos:
            Kr (list(float)): Lista com os pares (k_direta,k_reversa) de constantes reacionais.
            P (float): Pressão no reator.
            ReacCoefs (list(list(int))): Lista com os coeficientes reacionais para cada reação:
                                        [[R1 Coeficientes],[R2 Coeficientes],...]
            component_indx (int): Indice do componente cuja taxa reacional deve ser calculada.
            W_initial (list(float)): Lista com os valroes do chute inicial das composições de saída.
        Métodos:
            get_single_reaction_rate()
                : Calcula a taxa de reação para um componente i (component_indx) e uma ÚNICA reação j (reaction_number).
            get_global_reaction_rate()
                : Calcula a taxa global de reação para um componente i (component_indx)
                     iterando pelas j reações e chamando a função get_single_reaction_rate().
    """
    def __init__(self, Kr, reactionCoefficients,P,component_indx,W_initial):
        self.Kr=Kr
        self.ReacCoefs=reactionCoefficients
        self.P = P
        self.component_indx = component_indx
        self.W_initial = W_initial

    def get_single_reaction_rate(self,reaction_number):
        """
            Calcula a taxa de reação para um componente i (component_indx) e uma ÚNICA reação j (reaction_number).
            Argumentos:
                reaction_number (int): Indice j da reação cuja taxa será calculada.
            Retorna:
                (float) o valor da taxa da reação j para o componente i.
        """
        Pi=list()
        for wi in self.W_initial:
            Pi.append(wi*self.P)
        rdir = self.Kr[reaction_number][0]
        rinv = self.Kr[reaction_number][1]
        for i in range(len(self.ReacCoefs[reaction_number])):
            if self.ReacCoefs[reaction_number][i] < 0:
                rdir = rdir* (Pi[i]**(np.abs(self.ReacCoefs[reaction_number][i])))
            elif self.ReacCoefs[reaction_number][i] > 0:
                rinv = rinv* (Pi[i]**(np.abs(self.ReacCoefs[reaction_number][i])))
        ri = rdir -rinv
        return ri * self.ReacCoefs[reaction_number][self.component_indx]/np.abs(self.ReacCoefs[reaction_number][0])

    def get_global_reaction_rate(self):
        """
            Calcula a taxa global de reação para um componente i (component_indx) iterando pelas j reações e chamando a função get_single_reaction_rate().
            Retorna:
                (float) o valor global da taxa da reação i.
        """
        global_rate=0.0
        for reaction_number in range(len(self.ReacCoefs)):
            if np.abs(self.ReacCoefs[reaction_number][self.component_indx]) > 0.0:
                single_reaction_rate = self.get_single_reaction_rate(reaction_number)
                global_rate=global_rate+single_reaction_rate
        return global_rate


class ReactionRateConstant:
    """
        Responsável por calcular o valor das constantes reacionais usando a equação de Arrhenius.
        Argumentos:
            Ko (list(list(float))): Lista com os pares (ko_direta,ko_reversa) de constantes reacionais padrão.
            E (list(list(float))): Lista com os pares (Ea_direta,Ea_reversa) de energias de ativação das reações.
            T (float): Temperatura da reação.
        Atributos:
            Kr (list(list(float))): Lista com os pares (k_direta,k_reversa) de constantes reacionais a serem calculadas.
            R (float): Constante dos gases perfeitos (cal · K-1 · mol-1).
        Métodos:
            evaluate()
                : Realiza os calculos e atualiza os valores de Kr.
    """

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