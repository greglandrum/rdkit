import doctest
import gzip
import os
import sys
import unittest

from rdkit import Chem, RDConfig, __version__, rdBase

class TestCase(unittest.TestCase):
    def testDetermineFormat(self):
       testd = (
          ("foo.sdf","sdf",""),
          ("foo.sdf.gz","sdf","gz"),
          ("foo.sdfgz","sdf","gz"),
          ("foo.txt","txt",""),
          ("foo.smi","smi",""),
       )
       for tpl in testd:
          fmt,cfmt = Chem.DetermineFileFormat(tpl[0])
          self.assertEqual(fmt,tpl[1])
          self.assertEqual(cfmt,tpl[2])
          
       fails = ( "foo.sdq", "foo.tot", "foo.sdf.az")
       for fail in fails:
          with self.assertRaises(OSError):
             Chem.DetermineFileFormat(fail)


if __name__ == '__main__':
  unittest.main()
