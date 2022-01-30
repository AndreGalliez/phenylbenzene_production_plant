import numpy as np
from scipy.optimize import fsolve, newton

class Flash: ##Faz o cálculo de flash para quando as condições de equilibrio e composição na entrada são conhecidas
    """
        Responsável por realizar os cálculos referentes ao tanque de flash.
        Argumentos:
            P (float): Pressão de equilibrio no tanque de flash.
            Fin (float): Vazão de entrada.
            Z (numpy(float): Composições de entrada. Precisa ser um numpy, vai sofrer operações algébricas bidimensionais.
            P_sat (list(float)): Lista com o P_sat de cada componente. 
        Atributos:
            K (list(float)): Lista com as constantes de equilibrio líquido/vapor (P_sat/P) de cada componente.
            L (float): Vazão da corrente líquida de saída (a ser calculado).
            V (float): Vazão da corrente gasosa de saída (a ser calculado).
            X (list(float)): Composições da corrente líquida de saída (a ser calculado).
            Y (list(float)): Composições da corrente gasosa de saída (a ser calculado).
        Métodos:
            evaluate_K()
                : Faz as contas atualização o parâmetro K.
            formulate_equations_PT()
                : Formula a equação de Rashford-Rice a ser resolvida.
            evaluate_flash()
                : Resolve a equação de Rashford-Rice e calcula os valores das vazões e composições de saída.
    """

    def __init__(self, Fin, z, p_sat, P):
        self.P = P
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
        """
            Formula a equação de Rashford-Rice a ser resolvida que vai ser passada para um solver do scipy.
            Argumentos:
                x (float): Chute inicial para o parâmetro beta a ser encontrado.
            Retorna:
                Valor do resíduo da equação de Rashford-Rice a ser avaliado pelo solver do scipy .
        """
        f=0.0
        for i in range(len(self.Z)):
            f=f+self.Z[i]*self.K[i]/(1+x*(self.K[i]-1))
        return f-1

    def evaluate_flash_PT(self, x_in):
        """
            Resolve a equação de Rashford-Rice e calcula os valores das vazões e composições de saída.
            Chama a função formulate_equations_PT, passa o valor do resíduo para a função newton do scipy que resolve a equação não-linear.
            Na sequência calcula os atributos das correntes de saida e atualiza os parâmetros.
            Argumentos:
                x (float): Chute inicial para o parâmetro beta a ser encontrado.

        """
        self.evaluate_K()
        self.B = newton(self.formulate_equations_PT, x_in)
        for i in range(len(self.Z)):
            self.Y[i] = self.Z[i]*self.K[i]/(1+self.B*(self.K[i]-1))
        for i in range(len(self.Z)):
            self.X[i]= self.Y[i]/self.K[i]
        self.L=self.Fin*(1-self.B)
        self.V=self.Fin-self.L
    
   
class LiquidVaporEquilibriumConstant:
    """
        Responsável por calcular as pressoes de saturação de cada um dos componentes dada uma temperatura T.
        Argumentos:
            T (float): Temperatura de equilibrio no tanque de flash.
            elv_coefficients list((list(float))): Lista de listas com os coeficientes de cada componente a ser utilizado na correlação empírica usada para calcular P_sat.
                                                    Podem ser referentes à equação de Antoine ou Wagner. 
                                                    A equação de Antoine é usada para os componentes de A-C e a de Wagner para o componente D.
        Atributos:
            P_sat (list(float): Pressões de saturação de cada um dos componentes a serem calculadas.
        Métodos:
            antoine_method(singleComponentCoefficients)
                : Implementa a correlação de Antoine para calcular P_sat.
            wagner_method(singleComponentCoefficients)
                : Implementa a correlação de Wagner para calcular P_sat.
            calc_psats()
                : Realiza os calculos para todos os componentes e atualiza os valores de Fout e Wout.
    """
    
    def __init__(self, T, elv_coefficients):
        self.P_sat = None
        self.T=T
        self.elv_coefficients=elv_coefficients
            

    ##O metodo implementado usa Pressao em psia e temperatura em farenheit
    ##A equação de Antoine utilizada usa a Pressão crítica do componente como pressão de referência
    @staticmethod
    def antoine_method(T,singleComponentCoefficients):
        """
        Implementa a correlação de Antoine para calcular P_sat.
        Argumentos:
            singleComponentCoefficients (list(float)): Lista com os parametros da equação de um componente.
        Retorna:
            P_sat de um componente segundo os coeficientes fornecidos.
        """
        pi_sat = singleComponentCoefficients[3]*(np.exp(singleComponentCoefficients[0]-singleComponentCoefficients[1]/(T+singleComponentCoefficients[2])))
        return pi_sat
    
    ##O metodo implementado usa Pressao em Pascal e temperatura em Kelvin
    @staticmethod 
    def wagner_method(T,singleComponentCoefficients):
        """
        Implementa a correlação de Wagner para calcular P_sat.
        Argumentos:
            singleComponentCoefficients (list(float)): Lista com os parametros da equação de um componente.
        Retorna:
            P_sat de um componente segundo os coeficientes fornecidos.
        """
        pi_sat= 1000*(np.exp(singleComponentCoefficients[0]*np.log(T)+(singleComponentCoefficients[1]/T) +singleComponentCoefficients[2]+singleComponentCoefficients[3]*(T**2)))
        return pi_sat

    def calc_psats(self):
        P_sat=list()
        for i in range(len(self.elv_coefficients)):
            singleComponentCoefficients = self.elv_coefficients[i]
            if i in [0,1,2]:
                P_sat.append(6894.76*self.antoine_method(((self.T - 273.15) * (9/5) + 32),singleComponentCoefficients))
            else:
                P_sat.append(self.wagner_method(self.T,singleComponentCoefficients))
        self.P_sat = P_sat
        return P_sat
        ##Print das Pressoes e Constantes de Equilibrio
        #print(np.array(P_sat)/(10**5))
        #print(self.K)