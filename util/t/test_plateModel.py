import unittest
import plateModel as plate

class determineVelocities(unittest.TestCase):


    def test_one(self):
        sites = {}
        site_names = ('CEDU','KARR','YAR2')

        for name in site_names:
            sites[name] = {}
            sites[name]['X'] = []
            sites[name]['Y'] = []
            sites[name]['Z'] = []

        print(sites)
        sites['CEDU']['X'].append(-3753472.8476)
        sites['CEDU']['Y'].append(3912741.0128)
        sites['CEDU']['Z'].append(-3347960.1150) 

        sites['KARR']['X'].append(-2713832.9070)
        sites['KARR']['Y'].append(5303935.1072)
        sites['KARR']['Z'].append(-2269514.2088) 

        sites['YAR2']['X'].append(-2389026.2255)
        sites['YAR2']['Y'].append(5043316.9900)
        sites['YAR2']['Z'].append(-3078529.9783) 

        nuvel1A = plate.nuvel1APlateModel()
        print("plateModel:",nuvel1A['AUST'])
        asites = plate.determineVelocities(sites,nuvel1A['AUST'])
        
        self.assertAlmostEqual(asites['CEDU']['velX'][0],-0.04173,places=4)
        self.assertAlmostEqual(asites['CEDU']['velY'][0],0.00267,places=4)
        self.assertAlmostEqual(asites['CEDU']['velZ'][0],0.04990,places=4)
        #
        self.assertAlmostEqual(asites['KARR']['velX'][0],-0.04494,places=4)
        self.assertAlmostEqual(asites['KARR']['velY'][0],0.00074,places=4)
        self.assertAlmostEqual(asites['KARR']['velZ'][0],0.05547,places=4)
        #
        self.assertAlmostEqual(asites['YAR2']['velX'][0],-0.04745,places=4)
        self.assertAlmostEqual(asites['YAR2']['velY'][0],0.00912,places=4)
        self.assertAlmostEqual(asites['YAR2']['velZ'][0],0.05176,places=4)

# Positions input to Bernese
# CEDU -3753472.8476   3912741.0128  -3347960.1150 
# KARR -2713832.9070   5303935.1072  -2269514.2088 
# YAR2 -2389026.2255   5043316.9900  -3078529.9783 
# DARW -4091359.3114   4684606.5379  -1408579.6112 
# TID1 -4460996.6608   2682557.0851  -3674443.0169 
# DARM -4078219.3317   4709437.7089  -1363264.2995 
# 00NA -4073662.3047   4712064.7397  -1367874.4535 
# JAB2 -4236472.6817   4559859.4901  -1388764.5263 
# KUNU -3845767.2994   4789564.4930  -1712363.4225 
# DODA -4079102.4026   4661683.0535  -1515238.5009 
# LARR -4208011.3671   4479077.8760  -1701338.2019 
# KAT1 -4147413.5213   4581462.6974  -1573359.5572 
# ALIC -4052052.4096   4212836.0263  -2545105.0338 
# TOW2 -5054583.1345   3275504.2313  -2091538.9099 
# DWNI -4083215.3821   4704504.1465  -1365349.7260 

# Velocities obtained from Bernese
# NUvel1A-NNR velocities X Y Z (V M/Y)
# 00NA -0.03660 -0.01486  0.05781
# ALIC -0.03950 -0.00550  0.05378
# CEDU -0.04173  0.00267  0.04990
# DARM -0.03656 -0.01492  0.05781
# DARW -0.03663 -0.01465  0.05768
# DODA -0.03704 -0.01374  0.05744
# DWNI -0.03654 -0.01494  0.05780
# JAB2 -0.03575 -0.01572  0.05745
# KARR -0.04494  0.00074  0.05547
# KAT1 -0.03683 -0.01371  0.05716
# KUNU -0.03885 -0.01073  0.05724
# LARR -0.03684 -0.01309  0.05667
# TID1 -0.03568  0.00078  0.04389
# TOW2 -0.03129 -0.01535  0.05158
# YAR2 -0.04745  0.00912  0.05176


if __name__ == '__main__':
    unittest.main()

