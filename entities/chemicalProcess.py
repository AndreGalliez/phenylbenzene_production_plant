import numpy as np

class ChemicalProcess:

    def __init__(self,Frecycle_guess,Wrecycle_guess):
        """
        Responsável por determinar a ordem em que os equipamentos são calculados, 
        chamar o calculo e passar adiante os outputs de equipamentos que são inputs de outros.
        Argumentos:
            Frecycle_guess (float): Pressão de equilibrio no tanque de flash.
            Wrecycle_guess (float): Vazão de entrada.
        Atributos:
            F (list(float)): Lista de vazões de cada corrente do sistema (a ser calculado).
            W (list(list(float))): Lista de composições de cada corrente do sistema (a ser calculado).
            residual (float): Resíduo da iteração considerando a diferença entre valores iniciais e finais dos atributos da corrente de riclo.
        Métodos:
            calculate_mixer()
                : Instancia um objeto misturador com os parâmetros de entrada do processo, 
                    realiza os cálculos e incorpora sua corrente de saida nos atributos F e W.
            calculate_splitter()
                : Instancia um objeto separador com os parâmetros de entrada do processo,
                    realiza os cálculos e incorpora sua corrente de saida nos atributos F e W.
            get_reaction_constants()
                : Instancia um objeto ReactionRateConstant com os parâmetros de entrada do processo
                    e realiza os cálculos que define o valor das constantes reacionais a serem utilizadas no reator.
            calculate_reactor()
                : Instancia um objeto reator com os parâmetros de entrada do processo,
                    realiza os cálculos e incorpora sua corrente de saida nos atributos F e W.
            get_LVequilibrium_constant()
                : Instancia um objeto LiquidVaporEquilibriumConstant com os parâmetros de entrada do processo
                    e realiza os cálculos que define o valor das pressões de saturação a serem utilizadas no flash.
            calculate_flash()
                : Instancia um objeto flash com os parâmetros de entrada do processo,
                    realiza os cálculos e incorpora sua corrente de saida nos atributos F e W.
            evaluate_residual()
                : Calcula a norma do vetor diferença entre os atributos da corrente de riclo inicial e 
                a obtida após realizar os cálculos de processo. Atualiza o atrbuto residual.
            evaluate()
                : Chama todas as outras funções na ordem correta, organizando o passo a passo do processo.
                Funciona como a chamada para o cálculo.
    """
        self.Frecycle_guess = Frecycle_guess
        self.Wrecycle_guess = Wrecycle_guess
        self.F =[None] * 7
        self.W =[None] * 7
        self.residual = None  
    
    def calculate_mixer(self,Fin,Win,Frecycle_guess,Wrecycle_guess):
        from entities.connections import Mixer
        mixer=Mixer([Fin,Frecycle_guess],[Win,Wrecycle_guess])
        mixer.evaluate()
        self.F[1]=mixer.Fout
        self.W[1]=mixer.Wout

    def calculate_splitter(self,Fin,Win,Cs):
        from entities.connections import Splitter
        splitter=Splitter(Fin,Win,Cs)
        splitter.evaluate()
        self.F[6]=splitter.Fout['F_recycle']
        self.W[6]=splitter.Wout
        self.F[5]=splitter.Fout['F_purge']
        self.W[5]=splitter.Wout

    @staticmethod
    def get_reaction_constants(Ko,E,T):
        from entities.reactor import ReactionRateConstant
        reactionConstantSetter=ReactionRateConstant(Ko,E,T)
        reactionConstantSetter.evaluate_K()
        return reactionConstantSetter.Kr

    def calculate_reactor(self, Fin, Win, Vr, P, T, reactionCoefficients, Ko, E):
        from entities.reactor import GasPhaseReactor
        reactor=GasPhaseReactor(Fin, Win, Vr, self.get_reaction_constants(Ko,E,T), reactionCoefficients, P, T)
        reactor.evaluate((0.45,0.15,0.3,0.1,Fin)) ##initial guess for linear system 
        self.F[2]=reactor.Fout
        self.W[2]=reactor.Wout

    @staticmethod
    def get_LVequilibrium_constant(Tf, elv_coefficients):
        from entities.flash import LiquidVaporEquilibriumConstant
        equilibriumConstantSetter=LiquidVaporEquilibriumConstant(Tf, elv_coefficients)
        return equilibriumConstantSetter.calc_psats()

    def calculate_flash(self, Fin, Win, Tf, elv_coefficients, P):
        from entities.flash import Flash
        flash=Flash(Fin, Win, self.get_LVequilibrium_constant(Tf, elv_coefficients), P)
        flash.evaluate_flash_PT(0.6) ##initial guess for linear system 
        self.F[3]=flash.L
        self.W[3]=flash.X
        self.F[4]=flash.V
        self.W[4]=flash.Y

    def evaluate(self, Fin, Win,
                Vr, Pr, Tr, reactionCoefficients, Kor, Er,
                Pf, Tf, elv_coefficients, 
                Cs):
        self.F[0]=Fin
        self.W[0]=Win
        self.calculate_mixer(Fin,Win,self.Frecycle_guess,self.Wrecycle_guess)
        self.calculate_reactor(self.F[1], self.W[1], Vr, Pr, Tr, reactionCoefficients, Kor, Er)
        self.calculate_flash(self.F[2], self.W[2], Tf, elv_coefficients, Pf)
        self.calculate_splitter(self.F[4],self.W[4],Cs)
        self.evaluate_residual()

    def evaluate_residual(self):
        if self.F[6] == 0.0:
            self.residual = 0.0
        else:
            recycle_differences=list()
            for i in range(len(self.W[6])):
                recycle_differences.append((self.W[6][i]-self.Wrecycle_guess[i])/((self.W[6][i]+self.Wrecycle_guess[i])/2))
            recycle_differences.append((self.F[6]-self.Frecycle_guess)/((self.F[6]+self.Frecycle_guess)/2))
            self.residual = np.linalg.norm(recycle_differences)