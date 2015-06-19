# Copyright 2014-2015 Novartis Institutes for Biomedical Research

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import unittest
import tempfile, shutil, gzip

from railroadtracks import core
import railroadtracks.model.simulate

try:
    import ngs_plumbing
    has_ngsp = True
except ImportError:
    has_ngsp = False



class AssetsTestCase(unittest.TestCase):

    def setUp(self):
        self.tempdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tempdir)

    def test_assetfactory(self):
        Foo = core.assetfactory('Foo', [core.AssetAttr('bar', core.File, '')])
        # check that the asset bar is defined
        foo = Foo(core.File(''))
        # check that an incorrect type raises an error
        self.assertRaises(AssertionError, Foo, 123)
        self.assertTrue(foo.bar._defined)
        # check that incorrect parameters make it fail
        self.assertRaises(AssertionError,
                          core.assetfactory, 
                          'Foo', [core.AssetAttr('bar', core.File, None)])
        # check that trying to modify an asset raises an error
        self.assertRaises(AttributeError, setattr, foo, 'bar', 123)


    def test_AssetStep_NoTarget(self):
        Foo = core.assetfactory('Foo', [core.AssetAttr('bar', core.File, '')])
        # no target
        assetset = core.AssetsStep(Foo)
        self.assertTrue(isinstance(assetset.__str__(), str))

    def test_AssetStep_Target(self):
        Foo = core.assetfactory('Foo', [core.AssetAttr('bar', core.File, '')])
        Bar = core.assetfactory('Bar', [core.AssetAttr('foo', core.File, '')])
        assetset = core.AssetsStep(Foo, Bar)
        self.assertTrue(isinstance(assetset.__str__(), str))

    def test_assetfactory_allownone(self):
        # check that allownone=True allows unspecified assets
        Foo = core.assetfactory('Foo', [core.AssetAttr('bar', core.File, '', allownone=True)])
        # check that the asset bar is defined
        foo = Foo(core.File(''))
        self.assertTrue(foo.bar._defined)
        # check that an incorrect type raises an error
        self.assertRaises(AssertionError, Foo, 123)
        foo = Foo(None)
        self.assertTrue(foo.bar is None)
        # check that trying to modify an asset raises an error
        self.assertRaises(AttributeError, setattr, foo, 'bar', 123)

    def test_createundefined(self):
        # -- createundefined-begin
        Foo = core.assetfactory('Foo', [core.AssetAttr('bar', core.File, '')])
        # create an undefined set of assets of type Foo
        undefoo = Foo.createundefined()
        # -- createundefined-end
        self.assertTrue(isinstance(undefoo, Foo))
        # check that the asset bar is an undefined "saved entity"
        self.assertFalse(undefoo.bar._defined)
        # check that trying to modify an asset raises an error
        # (modifying an asset value will be possible though)
        self.assertRaises(AttributeError, setattr, undefoo, 'bar', 123)
        undefoo.bar.name = '123'
        self.assertEqual('123', undefoo.bar.name)

    @unittest.skipIf(not has_ngsp,
                     'The Python package ngs-plumbing is missing.')
    def test_GzipFastqFilePair(self):
        NFRAGMENTS_MATCH = 300
        read1_fh = tempfile.NamedTemporaryFile(prefix='read1', suffix='.fq.gz', dir=self.tempdir, delete=False)
        read1_fh.close()
        read2_fh = tempfile.NamedTemporaryFile(prefix='read2', suffix='.fq.gz', dir=self.tempdir, delete=False)
        read2_fh.close()
        with open(railroadtracks.model.simulate.PHAGEFASTA) as fasta_fh:
            reference = next(railroadtracks.model.simulate.readfasta_iter(fasta_fh))
            read1_io = gzip.GzipFile(read1_fh.name, mode='w')
            read2_io = gzip.GzipFile(read2_fh.name, mode='w')
            # reads from the target genome
            read1_fh, read2_fh=railroadtracks.model.simulate.randomPEreads(read1_io,
                                                                           read2_io,
                                                                           reference,
                                                                           n = NFRAGMENTS_MATCH)
        read1_fh.close()
        read2_fh.close()
        fqp = rnaseq.GzipFastqFilePair(read1_fh.name,
                                       read2_fh.name)
        readpairs = tuple(fqp)
        for i, (r1,r2) in enumerate(readpairs):
            self.assertTrue(hasattr(r1, 'header'))
            self.assertTrue(hasattr(r1, 'sequence'))
            self.assertTrue(hasattr(r1, 'quality'))
            self.assertTrue(hasattr(r2, 'header'))
            self.assertTrue(hasattr(r2, 'sequence'))
            self.assertTrue(hasattr(r2, 'quality'))
        self.assertEqual(NFRAGMENTS_MATCH, i+1)

if __name__ == '__main__':
    unittest.main()
