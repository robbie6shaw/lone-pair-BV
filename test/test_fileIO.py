import unittest
import numpy as np
import pandas as pd
import calculate, fileIO
import logging, math
from test_BVStructure import alteredTestCase

class TestBVProdDatabase (alteredTestCase):

    def setUp(self):
        self.db = fileIO.BVDatabase("soft-bv-params.sqlite3")
    
    def test_get_bv_params(self):

        self.db.get_bv_params("Sn2+", "F1-")
        self.db.get_bv_params(('Sn', 2), ('F', -1))

    def test_calc_bv_params(self):
        
        params = self.db.get_bv_params("Na+", "O2-", True)
        self.assertAlmostEqual(params.d0, 0.57523)
        self.assertAlmostEqual(params.rmin, 2.37433)

    def test_d0_rmin_calc(self):

        precalc = pd.read_excel("results/BVSE-params.xlsx")

        for row in precalc.itertuples():
            params = self.db.get_bv_params(row.cation, "O2-", True)
            self.assertWithinBounds(params.rmin, row.rmin, 0.5)
            self.assertWithinBounds(params.d0, row.d0, 0.5)