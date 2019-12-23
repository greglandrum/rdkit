#
#  Copyright (C) 2019 Greg Landrum
#         All Rights Reserved
#
from rdkit import RDConfig, rdBase
from rdkit import Chem
from io import BytesIO
import unittest


class TestCase(unittest.TestCase):
    def test1Basics(self):
        m = Chem.MolFromSmiles('CC(C)[C@](F)(O)Cl')
        rnks = Chem.ChiralRankMolAtoms(m)
        self.assertEqual(list(rnks),[0, 2, 0, 3, 5, 4, 6])
        
        m = Chem.MolFromSmiles('CC[C@](F)(O)C=C')
        rnks = Chem.ChiralRankMolAtoms(m)
        #print(list(rnks))
        self.assertEqual(list(rnks),[0, 2, 4, 6, 5, 3, 1])


if __name__ == '__main__':
  unittest.main()
