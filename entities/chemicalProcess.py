import numpy as np

class ChemicalProcess:

    def __init__(self,Frecycle_guess,Wrecycle_guess):
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

    def calculate_reactor(self, Fin, Win, Vr, P, T, reactionCoefficients, reactionModel, Ko, E):
        from entities.reactor import GasPhaseReactor
        reactor=GasPhaseReactor(Fin, Win, Vr, self.get_reaction_constants(Ko,E,T), reactionCoefficients, reactionModel, P, T)
        reactor.evaluate((0.45,0.15,0.3,0.1,Fin)) ##initial guess for linear system 
        self.F[2]=reactor.Fout
        self.W[2]=reactor.Wout

    @staticmethod
    def get_LVequilibrium_constant(equilibrum_model, model_inputs):
        from entities.separationprocesses import LiquidVaporEquilibriumConstant
        equilibriumConstantSetter=LiquidVaporEquilibriumConstant(equilibrum_model, model_inputs)
        return equilibriumConstantSetter.calc_psats()

    def calculate_flash(self, Fin, Win, equilibrum_model, model_inputs, P):
        from entities.separationprocesses import Flash
        flash=Flash("PT", Fin, Win, self.get_LVequilibrium_constant(equilibrum_model, model_inputs), P)
        flash.evaluate_flash_PT(0.6) ##initial guess for linear system 
        self.F[3]=flash.L
        self.W[3]=flash.X
        self.F[4]=flash.V
        self.W[4]=flash.Y

    def evaluate(self, Fin, Win,
                Vr, Pr, Tr, reactionCoefficients, reactionModel, Kor, Er,
                lv_equilibrum_model, Pe, lv_model_inputs, 
                Cs):
        self.F[0]=Fin
        self.W[0]=Win
        self.calculate_mixer(Fin,Win,self.Frecycle_guess,self.Wrecycle_guess)
        self.calculate_reactor(self.F[1], self.W[1], Vr, Pr, Tr, reactionCoefficients, reactionModel, Kor, Er)
        self.calculate_flash(self.F[2], self.W[2], lv_equilibrum_model, lv_model_inputs, Pe)
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