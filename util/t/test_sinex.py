import unittest
import sinex as snx

class ell2xyzTEST(unittest.TestCase):

    def test_one(self):
        sinex,cova = snx.readSINEX('./t/ITRF2008-TRF-AUS.SNX')
        #
        #=====================================================================
        #
        auck = sinex['AUCK']
        self.assertAlmostEqual(auck['STAX'][0],-5105681.18110585,places=4)
        self.assertAlmostEqual(auck['STAY'][0],461564.034836750,places=4)
        self.assertAlmostEqual(auck['STAZ'][0],-3782181.48480215,places=4)
        #========
        self.assertAlmostEqual(auck['sSTAX'][0],0.79349E-03,places=4)
        self.assertAlmostEqual(auck['sSTAY'][0],0.55116E-03,places=4)
        self.assertAlmostEqual(auck['sSTAZ'][0],0.77486E-03,places=4)
        #========
        self.assertAlmostEqual(auck['VELX'][0],-.234727413865940E-01,places=4)
        self.assertAlmostEqual(auck['VELY'][0],-.240930121962410E-02,places=4)
        self.assertAlmostEqual(auck['VELZ'][0],0.324741037497387E-01,places=4)
        #========
        self.assertAlmostEqual(auck['sVELX'][0],0.78966E-04,places=4)
        self.assertAlmostEqual(auck['sVELY'][0],0.37008E-04,places=4)
        self.assertAlmostEqual(auck['sVELZ'][0],0.66142E-04,places=4)
        #========
        STAX = auck['STAZ_index']
        STAY = auck['STAY_index']
        STAZ = auck['STAZ_index']
        self.assertAlmostEqual(cova[STAX,STAX],0.629620028978910E-06,places=4)
        self.assertAlmostEqual(cova[STAY,STAY],0.303775016255092E-06,places=4)
        self.assertAlmostEqual(cova[STAZ,STAZ],0.600407541321165E-06,places=4)
        #========
        VELX = auck['VELX_index']
        VELY = auck['VELY_index']
        VELZ = auck['VELZ_index']
        self.assertAlmostEqual(cova[VELX,VELX],0.623557741031986E-08,places=4)
        self.assertAlmostEqual(cova[VELY,VELY],0.136957372931789E-08,places=4)
        self.assertAlmostEqual(cova[VELZ,VELZ],0.437481351831506E-08,places=4)
        #
        #=====================================================================
        #
        yar1 = sinex['YAR1']
        self.assertAlmostEqual(yar1['STAX'][0],-.238902590669338E+07,places=4)
        self.assertAlmostEqual(yar1['STAY'][0],0.504331693170434E+07,places=4)
        self.assertAlmostEqual(yar1['STAZ'][0],-.307853031466666E+07,places=4)
        #========
        self.assertAlmostEqual(yar1['sSTAX'][0],0.69473E-03,places=4)
        self.assertAlmostEqual(yar1['sSTAY'][0],0.86117E-03,places=4)
        self.assertAlmostEqual(yar1['sSTAZ'][0],0.78025E-03,places=4)
        #========
        self.assertAlmostEqual(yar1['VELX'][0],-.472518852152992E-01,places=4)
        self.assertAlmostEqual(yar1['VELY'][0],0.874632885656826E-02,places=4)
        self.assertAlmostEqual(yar1['VELZ'][0],0.504807997383417E-01,places=4)
        #========
        self.assertAlmostEqual(yar1['sVELX'][0],0.48406E-04,places=4)
        self.assertAlmostEqual(yar1['sVELY'][0],0.75481E-04,places=4)
        self.assertAlmostEqual(yar1['sVELZ'][0],0.57050E-04,places=4)
        #========
        STAX = yar1['STAZ_index']
        STAY = yar1['STAY_index']
        STAZ = yar1['STAZ_index']
        self.assertAlmostEqual(cova[STAX,STAX],0.482646230517592E-06,places=4)
        self.assertAlmostEqual(cova[STAY,STAY],0.741617781604912E-06,places=4)
        self.assertAlmostEqual(cova[STAZ,STAZ],0.608791002620522E-06,places=4)
        #========
        VELX = yar1['VELX_index']
        VELY = yar1['VELY_index']
        VELZ = yar1['VELZ_index']
        self.assertAlmostEqual(cova[VELX,VELX],0.234316187548364E-08,places=4)
        self.assertAlmostEqual(cova[VELY,VELY],0.569734362762740E-08,places=4)
        self.assertAlmostEqual(cova[VELZ,VELZ],0.325474278760286E-08,places=4)
        #=====================================================================
        #
        hob2 = sinex['HOB2']
        self.assertAlmostEqual(hob2['STAX'][0],-.395007167349903E+07,places=4)
        self.assertAlmostEqual(hob2['STAY'][0],0.252241525416228E+07,places=4)
        self.assertAlmostEqual(hob2['STAZ'][0],-.431163802558521E+07,places=4)
        #========
        self.assertAlmostEqual(hob2['sSTAX'][0],0.72176E-03,places=4)
        self.assertAlmostEqual(hob2['sSTAY'][0],0.60632E-03,places=4)
        self.assertAlmostEqual(hob2['sSTAZ'][0],0.80683E-03,places=4)
        #========
        self.assertAlmostEqual(hob2['VELX'][0],-.397396684594973E-01,places=4)
        self.assertAlmostEqual(hob2['VELY'][0],0.862000905404151E-02,places=4)
        self.assertAlmostEqual(hob2['VELZ'][0],0.407409725986825E-01,places=4)
        #========
        self.assertAlmostEqual(hob2['sVELX'][0],0.66860E-04,places=4)
        self.assertAlmostEqual(hob2['sVELY'][0],0.50939E-04,places=4)
        self.assertAlmostEqual(hob2['sVELZ'][0],0.73296E-04,places=4)
        #========
        STAX = hob2['STAZ_index']
        STAY = hob2['STAY_index']
        STAZ = hob2['STAZ_index']
        self.assertAlmostEqual(cova[STAX,STAX],0.520943082159015E-06,places=4)
        self.assertAlmostEqual(cova[STAY,STAY],0.367627183585689E-06,places=4)
        self.assertAlmostEqual(cova[STAZ,STAZ],0.650975735034499E-06,places=4)
        #========
        VELX = hob2['VELX_index']
        VELY = hob2['VELY_index']
        VELZ = hob2['VELZ_index']
        self.assertAlmostEqual(cova[VELX,VELX],0.447020891507436E-08,places=4)
        self.assertAlmostEqual(cova[VELY,VELY],0.259476595303455E-08,places=4)
        self.assertAlmostEqual(cova[VELZ,VELZ],0.537231066231415E-08,places=4)
        #========

    # remove a station, and make sure the coavriances are still correct
    def test_two(self):
        sinex,cova = snx.readSINEX('./t/ITRF2008-TRF-AUS.SNX')
        sinex,cova = snx.removeSolution(['STR1'],sinex,cova)
        #
        #=====================================================================
        #
        auck = sinex['AUCK']
        self.assertAlmostEqual(auck['STAX'][0],-5105681.18110585,places=4)
        self.assertAlmostEqual(auck['STAY'][0],461564.034836750,places=4)
        self.assertAlmostEqual(auck['STAZ'][0],-3782181.48480215,places=4)
        #========
        self.assertAlmostEqual(auck['sSTAX'][0],0.79349E-03,places=4)
        self.assertAlmostEqual(auck['sSTAY'][0],0.55116E-03,places=4)
        self.assertAlmostEqual(auck['sSTAZ'][0],0.77486E-03,places=4)
        #========
        self.assertAlmostEqual(auck['VELX'][0],-.234727413865940E-01,places=4)
        self.assertAlmostEqual(auck['VELY'][0],-.240930121962410E-02,places=4)
        self.assertAlmostEqual(auck['VELZ'][0],0.324741037497387E-01,places=4)
        #========
        self.assertAlmostEqual(auck['sVELX'][0],0.78966E-04,places=4)
        self.assertAlmostEqual(auck['sVELY'][0],0.37008E-04,places=4)
        self.assertAlmostEqual(auck['sVELZ'][0],0.66142E-04,places=4)
        #========
        STAX = auck['STAZ_index']
        STAY = auck['STAY_index']
        STAZ = auck['STAZ_index']
        self.assertAlmostEqual(cova[STAX,STAX],0.629620028978910E-06,places=4)
        self.assertAlmostEqual(cova[STAY,STAY],0.303775016255092E-06,places=4)
        self.assertAlmostEqual(cova[STAZ,STAZ],0.600407541321165E-06,places=4)
        #========
        VELX = auck['VELX_index']
        VELY = auck['VELY_index']
        VELZ = auck['VELZ_index']
        self.assertAlmostEqual(cova[VELX,VELX],0.623557741031986E-08,places=4)
        self.assertAlmostEqual(cova[VELY,VELY],0.136957372931789E-08,places=4)
        self.assertAlmostEqual(cova[VELZ,VELZ],0.437481351831506E-08,places=4)
        #
        #=====================================================================
        #
        yar1 = sinex['YAR1']
        self.assertAlmostEqual(yar1['STAX'][0],-.238902590669338E+07,places=4)
        self.assertAlmostEqual(yar1['STAY'][0],0.504331693170434E+07,places=4)
        self.assertAlmostEqual(yar1['STAZ'][0],-.307853031466666E+07,places=4)
        #========
        self.assertAlmostEqual(yar1['sSTAX'][0],0.69473E-03,places=4)
        self.assertAlmostEqual(yar1['sSTAY'][0],0.86117E-03,places=4)
        self.assertAlmostEqual(yar1['sSTAZ'][0],0.78025E-03,places=4)
        #========
        self.assertAlmostEqual(yar1['VELX'][0],-.472518852152992E-01,places=4)
        self.assertAlmostEqual(yar1['VELY'][0],0.874632885656826E-02,places=4)
        self.assertAlmostEqual(yar1['VELZ'][0],0.504807997383417E-01,places=4)
        #========
        self.assertAlmostEqual(yar1['sVELX'][0],0.48406E-04,places=4)
        self.assertAlmostEqual(yar1['sVELY'][0],0.75481E-04,places=4)
        self.assertAlmostEqual(yar1['sVELZ'][0],0.57050E-04,places=4)
        #========
        STAX = yar1['STAZ_index']
        STAY = yar1['STAY_index']
        STAZ = yar1['STAZ_index']
        self.assertAlmostEqual(cova[STAX,STAX],0.482646230517592E-06,places=4)
        self.assertAlmostEqual(cova[STAY,STAY],0.741617781604912E-06,places=4)
        self.assertAlmostEqual(cova[STAZ,STAZ],0.608791002620522E-06,places=4)
        #========
        VELX = yar1['VELX_index']
        VELY = yar1['VELY_index']
        VELZ = yar1['VELZ_index']
        self.assertAlmostEqual(cova[VELX,VELX],0.234316187548364E-08,places=4)
        self.assertAlmostEqual(cova[VELY,VELY],0.569734362762740E-08,places=4)
        self.assertAlmostEqual(cova[VELZ,VELZ],0.325474278760286E-08,places=4)
        #=====================================================================
        #
        hob2 = sinex['HOB2']
        self.assertAlmostEqual(hob2['STAX'][0],-.395007167349903E+07,places=4)
        self.assertAlmostEqual(hob2['STAY'][0],0.252241525416228E+07,places=4)
        self.assertAlmostEqual(hob2['STAZ'][0],-.431163802558521E+07,places=4)
        #========
        self.assertAlmostEqual(hob2['sSTAX'][0],0.72176E-03,places=4)
        self.assertAlmostEqual(hob2['sSTAY'][0],0.60632E-03,places=4)
        self.assertAlmostEqual(hob2['sSTAZ'][0],0.80683E-03,places=4)
        #========
        self.assertAlmostEqual(hob2['VELX'][0],-.397396684594973E-01,places=4)
        self.assertAlmostEqual(hob2['VELY'][0],0.862000905404151E-02,places=4)
        self.assertAlmostEqual(hob2['VELZ'][0],0.407409725986825E-01,places=4)
        #========
        self.assertAlmostEqual(hob2['sVELX'][0],0.66860E-04,places=4)
        self.assertAlmostEqual(hob2['sVELY'][0],0.50939E-04,places=4)
        self.assertAlmostEqual(hob2['sVELZ'][0],0.73296E-04,places=4)
        #========
        STAX = hob2['STAZ_index']
        STAY = hob2['STAY_index']
        STAZ = hob2['STAZ_index']
        self.assertAlmostEqual(cova[STAX,STAX],0.520943082159015E-06,places=4)
        self.assertAlmostEqual(cova[STAY,STAY],0.367627183585689E-06,places=4)
        self.assertAlmostEqual(cova[STAZ,STAZ],0.650975735034499E-06,places=4)
        #========
        VELX = hob2['VELX_index']
        VELY = hob2['VELY_index']
        VELZ = hob2['VELZ_index']
        self.assertAlmostEqual(cova[VELX,VELX],0.447020891507436E-08,places=4)
        self.assertAlmostEqual(cova[VELY,VELY],0.259476595303455E-08,places=4)
        self.assertAlmostEqual(cova[VELZ,VELZ],0.537231066231415E-08,places=4)
        #========

if __name__ == '__main__':
    unittest.main()

