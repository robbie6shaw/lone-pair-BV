import unittest
import numpy as np
import pandas as pd
import bv2, fileIO
import logging, math
from test_BVStructure import alteredTestCase

class TestBVProdDatabase (alteredTestCase):

    def setUp(self):
        self.db = fileIO.BVDatabase("soft-bv-params.sqlite3")
    
    def test_get_bv_params(self):

        self.db.getParams("Sn2+", "F1-")
        self.db.getParams(('Sn', 2), ('F', -1))