import unittest
import geodetic as gd

#class TestStringMethods(unittest.TestCase):
#
#   def test_isupper(self):
#       self.assertTrue('FOO'.isupper())
#       self.assertFalse('Foo'.isupper())
#
#    def test_split(self):
#        s = 'hello world'
#        self.assertEqual(s.split(), ['hello', 'world'])
#        # check that s.split fails when the separator is not a string
#        with self.assertRaises(TypeError):
#            s.split(2)

class ell2xyzTEST(unittest.TestCase):

    def test_one(self):
        # Use default WGS84
        self.assertAlmostEqual(gd.llh2xyz(-45.,140.,0.)[0],-3460675.389,places=3)
        self.assertAlmostEqual(gd.llh2xyz(-45.,140.,0.)[1],2903851.443,places=3)
        self.assertAlmostEqual(gd.llh2xyz(-45.,140.,0.)[2],-4487348.409,places=3)
        #
        self.assertAlmostEqual(gd.llh2xyz(-41.,140.,0.)[0],-3692786.960,places=3)
        self.assertAlmostEqual(gd.llh2xyz(-41.,140.,0.)[1],3098616.176,places=3)
        self.assertAlmostEqual(gd.llh2xyz(-41.,140.,0.)[2],-4162423.201,places=3)
        #
        self.assertAlmostEqual(gd.llh2xyz(-25.,140.,0.)[0],-4430811.872,places=3)
        self.assertAlmostEqual(gd.llh2xyz(-25.,140.,0.)[1],3717892.607,places=3)
        self.assertAlmostEqual(gd.llh2xyz(-25.,140.,0.)[2],-2679074.463,places=3)
        #
        self.assertAlmostEqual(gd.llh2xyz(19.,140.,0.)[0],-4621383.516,places=3)
        self.assertAlmostEqual(gd.llh2xyz(19.,140.,0.)[1],3877801.204,places=3)
        self.assertAlmostEqual(gd.llh2xyz(19.,140.,0.)[2],2063349.463,places=3)
        #
        self.assertAlmostEqual(gd.llh2xyz(88.,140.,0.)[0],-171089.653,places=3)
        self.assertAlmostEqual(gd.llh2xyz(88.,140.,0.)[1],143561.265,places=3)
        self.assertAlmostEqual(gd.llh2xyz(88.,140.,0.)[2],6352853.879,places=3)

    def test_two(self):
        #print(gd.xyz2llh(-171089.653,143561.265,6352853.879))
        self.assertAlmostEqual(gd.xyz2llh(-171089.653,143561.265,6352853.879)[0],88.0,places=3)
        self.assertAlmostEqual(gd.xyz2llh(-171089.653,143561.265,6352853.879)[1],140.0,places=3)
        self.assertAlmostEqual(gd.xyz2llh(-171089.653,143561.265,6352853.879)[2],0.0,places=3)

    def test_three(self):
        xyz = [-171089.653, 143561.265, 6352853.879]
        neu = gd.geo2topo(88.,140.,xyz)
        XYZ = gd.topo2geo(88.,140.,neu)
        self.assertAlmostEqual((xyz[0] - XYZ[0]),0,places=5)
        self.assertAlmostEqual((xyz[1] - XYZ[1]),0,places=5)
        self.assertAlmostEqual((xyz[2] - XYZ[2]),0,places=5)

    # convert llh to neu and then subtract the two to get the difference in metres
    def test_four(self):

if __name__ == '__main__':
    unittest.main()

