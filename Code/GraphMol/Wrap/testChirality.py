#
#  Copyright (C) 2019 Greg Landrum
#         All Rights Reserved
#
from rdkit import RDConfig, rdBase
from rdkit import Chem
from io import BytesIO
import unittest


class TestCase(unittest.TestCase):
    def _compareRankings(self,rnks,ornks):
        tmp = sorted(set(rnks))
        rnks = [tmp.index(x) for x in rnks]
        tmp = sorted(set(ornks))
        ornks = [tmp.index(x) for x in ornks]        
        self.assertEqual(list(rnks),list(ornks))

    def test1Basics(self):
        m = Chem.MolFromSmiles('CC(C)[C@](F)(O)Cl')
        rnks = Chem.ChiralRankMolAtoms(m)
        self.assertEqual(list(rnks),[0, 2, 0, 3, 5, 4, 6])
        
        m = Chem.MolFromSmiles('CC[C@](F)(O)C=C')
        rnks = Chem.ChiralRankMolAtoms(m)
        #print(list(rnks))
        self.assertEqual(list(rnks),[0, 2, 4, 6, 5, 3, 1])

    def test2DependantChirality(self):
        m = Chem.MolFromSmiles('C[C@](O)(F)[C@H](C)[C@](C)(O)F')
        rnks = Chem.ChiralRankMolAtoms(m)
        #print(list(rnks))
        self.assertGreater(rnks[6],rnks[1])

    def test3CompareToOld(self):
        m = Chem.MolFromSmiles('CC(C)[C@](F)(O)Cl')
        rnks = list(Chem.ChiralRankMolAtoms(m))
        ornks = list(Chem.CIPRankMolAtoms(m))
        self._compareRankings(rnks,ornks)

        m = Chem.MolFromSmiles('CC[C@](F)(O)C=C')
        rnks = list(Chem.ChiralRankMolAtoms(m))
        ornks = list(Chem.CIPRankMolAtoms(m))
        self._compareRankings(rnks,ornks)


if __name__ == '__main__':
  unittest.main()
