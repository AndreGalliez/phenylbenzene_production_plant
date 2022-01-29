import unittest
import numpy as np
from connections import Splitter, Mixer
from reactor import GasPhaseReactor, ReactionRateConstant
from flash import Flash,LiquidVaporEquilibriumConstant

class TestConnections(unittest.TestCase):

    def test_splitter(self):
        splitter = Splitter(100,[0.25,0.25,0.25,0.25],0.6)
        splitter.evaluate()
        self.assertEqual(splitter.Fout['F_recycle'], 60.0)
        self.assertEqual(splitter.Fout['F_purge'], 40.0)
        self.assertEqual(splitter.Wout, [0.25,0.25,0.25,0.25])

    def test_mixer(self):
        mixer = Mixer([100,50],[[1.0,0.0,0.0,0.0],[0.0,0.5,0.5,0.0]])
        mixer.evaluate()
        self.assertEqual(mixer.Fout, 150.0)
        self.assertEqual(mixer.Wout, [4/6,1/6,1/6,0])

class TestReactor(unittest.TestCase):

    def test_reaction_rate(self):
        kmd = ReactionRateConstant([[3600*3.25*10**(-6), 3600*1.0205*10**(-5)], 
                                    [3600*3.7545*10**(-6), 3600*7.9544*10**(-6)]],[[30190, 30190],
                                                                                    [30190,30190]],1038.262085)
        kmd.evaluate_K()
        flatten_kmd = [item for sublist in kmd.Kr for item in sublist]
        expected_wi_result = [5.16928270*10**(-9), 1.62315477*10**(-8),5.97171444*10**(-9), 1.26518592*10**(-8)]
        for i in range(len(expected_wi_result)):
            self.assertAlmostEqual(flatten_kmd[i],expected_wi_result[i], places=5)

    def test_reactor(self):
        reactor = GasPhaseReactor(100,[1,0,0,0],1,[[5.16928270*10**(-9), 1.62315477*10**(-8)],[5.97171444*10**(-9), 1.26518592*10**(-8)]],[[-2,1,1,0],[-1,-1,1,1]],"TwoWaysPowerlaw",12*10.0**(5),1038.262085)
        reactor.evaluate((0.40,0.2,0.3,0.1,100))
        self.assertAlmostEqual(reactor.Fout,100)
        expected_wi_result = [0.41582426620362845, 0.15708837073377083, 0.3370876989529908, 0.08999966410960998]
        for i in range(len(expected_wi_result)):
            self.assertAlmostEqual(reactor.Wout[i],expected_wi_result[i], places=5)

class TestFlash(unittest.TestCase):

    def test_p_sat_calculation(self):
        coefficients = [[-8.43361,-6281.04,71.1072,6.1984*10**(-6)],[-11.1711,-10795.2,94.0474,3.8631*10**(-6)],[5.602657,418.1773,474.214,190.8],[-14.7697,-15484.2,122.524,3.7852*10**(-6)]]
        eq_LV = LiquidVaporEquilibriumConstant('Raoult and Other',[573,coefficients])
        eq_LV.calc_psats()
        expected_psat_result = [0.16275608612538364, 4.128048676381459, 239147540.2022536, 18974.149714045023]
        for i in range(len(expected_psat_result)):
            self.assertAlmostEqual(eq_LV.P_sat[i],expected_psat_result[i], places=5)

    def test_reactor(self):
        z= np.array([0.42206264789235637, 0.14409646369486045, 0.3078262669389916, 0.08382562981951494])
        flash=Flash("PT",250,z,[0.16275608612538364, 4.128048676381459, 239147540.2022536, 18974.149714045023],30*(10**5))
        flash.evaluate_flash_PT(0.3)
        self.assertAlmostEqual(flash.L,175.18288838424246, places=5)
        self.assertAlmostEqual(flash.V,74.81711161575754, places=5)
        expected_yi_result = [3.2676927767998985e-08, 2.8295993601598507e-07, 0.9992451227593246, 0.0007545616036956839]
        expected_xi_result = [0.6023171583794183, 0.20563706356099978, 0.012535087610529915, 0.11930362336138992] 
        for i in range(len(expected_yi_result)):
            self.assertAlmostEqual(flash.Y[i],expected_yi_result[i], places=5)
            self.assertAlmostEqual(flash.X[i],expected_xi_result[i], places=5)

if __name__ == '__main__':
    unittest.main()