import unittest

import os
import sys
sys.path.insert(0,os.path.abspath('..'))

import gpsTime as gt

class test_yyyy2yy(unittest.TestCase):

    def test_one(self):
        self.assertEqual(gt.yyyy2yy(1984),84)
        self.assertEqual(gt.yyyy2yy(1999),99)
        self.assertEqual(gt.yyyy2yy(2000),00)
        self.assertEqual(gt.yyyy2yy(2015),15)
        self.assertEqual(gt.yyyy2yy(2025),25)

class test_yy2yyyy(unittest.TestCase):
    def test_one(self):
        self.assertEqual(gt.yy2yyyy(84),1984)
        self.assertEqual(gt.yy2yyyy(99),1999)
        self.assertEqual(gt.yy2yyyy(00),2000)
        self.assertEqual(gt.yy2yyyy(15),2015)
        self.assertEqual(gt.yy2yyyy(25),2025)

class test_wwww2doys(unittest.TestCase):
    def test_one(self):
        self.assertEqual(gt.wwww2doys(1772),['356','357','358','359','360','361','362'])

class test_yyyydoy2cal(unittest.TestCase):
    def test_one(self):
        self.assertEqual(gt.yyyydoy2cal(2006,228),(2006,8,16))
        self.assertEqual(gt.yyyydoy2cal(2016,256),(2016,9,12))

class test_ymd2yyyyddd(unittest.TestCase):
    def test_one(self):
        self.assertEqual(gt.ymd2yyyyddd(2016,9,18),(2016,262))
        self.assertEqual(gt.ymd2yyyyddd(1999,3,4),(1999,63))

if __name__ == '__main__':
    unittest.main()
