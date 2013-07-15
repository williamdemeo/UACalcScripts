'''
Created on Jul 14, 2013

@author: williamdemeo
'''
import unittest
from aljebra.closure import Closure

class TestAljebra(unittest.TestCase):


    def setUp(self):
        pass


    def tearDown(self):
        pass


    def testsd_embedding(self,pars):
        closure = Closure()
        sd_embedding = closure.sd_embedding()
        pass


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_sd_embedding']
    unittest.main()