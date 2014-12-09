# Copyright 2014 Novartis Institutes for Biomedical Research

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
import os, tempfile
import subprocess
import collections
import itertools
import csv
import sqlite3
import railroadtracks
import railroadtracks.model.simulate
import railroadtracks.model.aligners
import railroadtracks.model.diffexp
import railroadtracks.model.files
import shutil
import gzip
import re
import random
from . import core
from . import hortator
from . import rnaseq
from . import environment
from . import unifex
from . import easy

has_R = environment.Executable.ispresent('R')

_DEBUGGER = True
PHAGEGFF = railroadtracks.model.simulate.PHAGEGFF
PHAGEGTF = railroadtracks.model.simulate.PHAGEGTF
PHAGEFASTA = railroadtracks.model.simulate.PHAGEFASTA


class EnvironmentTestCase(unittest.TestCase):

    def test_executable(self):
        ls_exec = environment.Executable('ls')
        #FIXME: technically no test here

    @unittest.skipIf(not has_R,
                     'R is missing')
    def test_R(self):
        r_exec = environment.R('R')
        r_version = r_exec.version

        # missing package
        self.assertRaises(ValueError, r_exec.packageversion, 'foobarbaz')

        version = r_exec.packageversion('stats')
        self.assertTrue(r_version.startswith(version))

    @unittest.skipIf(not (has_R and \
                          environment.R('R').packageversion_or_none('rjson') is not None),
                     'R and its package "rjson" must be present.')
    def test_R_run_snippet(self):
        r_exec = environment.R('R')
        # run a snippet
        magicvariable = 'railroadtracks_import'
        code = """        
        run <- function(p) {
          res <- 1 + p$x
          # output file
          fh_out <- file(p$output_fn, "w")
          write(res, file = fh_out)
          flush(fh_out)
        }
        run(%s)
        """ % magicvariable
        with tempfile.NamedTemporaryFile() as fh_out:
            var_in = {'x': 3,
                      'output_fn': fh_out.name}
            r_exec.run_snippet(code, var_in)
            res = fh_out.read()
        self.assertEquals('4', res.rstrip())


class ModelSimulateTestCase(unittest.TestCase):

    def setUp(self):
        self._dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self._dir)

    def _test_FluxsimulatorExpression(self, fse):
        A = fse.Assets
        annotation_fn = PHAGEGTF
        genome_dir = os.path.dirname(annotation_fn)

        params_fn = os.path.join(self._dir, 'original.par')
        with open(params_fn, 'w') as fh_out:
            fh_out.write("""
FRAG_SUBSTRATE  RNA
FRAG_METHOD     UR
""")

        paramsandpro_prefix = os.path.join(self._dir, 'modified')
        assets = A(A.Source(rnaseq.SavedGFF(annotation_fn),
                            rnaseq.FilePattern(genome_dir),
                            rnaseq.FluxSimulatorParameters(params_fn)), 
                   A.Target(rnaseq.FluxSimulatorParametersAndPro(paramsandpro_prefix)))
        return assets

    @unittest.skipIf(not environment.Executable.ispresent('flux-simulator'),
                     'flux-simulator is not in the PATH')
    def test_FluxsimulatorExpression(self):
        fse = rnaseq.FluxsimulatorExpression()
        assets = self._test_FluxsimulatorExpression(fse)
        cmd, returncode = fse.run(assets, parameters=('--NB-MOLECULES', '1000'))
        self.assertEquals(0, returncode)
        # check that the target file is not empty
        self.assertGreater(os.stat(assets.target.parameters_and_pro.name + '.pro').st_size, 0)

    @unittest.skipIf(not environment.Executable.ispresent('flux-simulator'),
                     'flux-simulator is not in the PATH')
    def test_FluxsimulatorSequencing(self):
        # expression
        fse = rnaseq.FluxsimulatorExpression()
        fse_assets = self._test_FluxsimulatorExpression(fse)
        cmd, returncode = fse.run(fse_assets, parameters=('--NB-MOLECULES', '1000'))
        read1_fn = os.path.join(self._dir, 'read1.fq.gz')
        read2_fn = os.path.join(self._dir, 'read2.fq.gz')

        fss = rnaseq.FluxsimulatorSequencing()
        A = fss.Assets
        paramsandpro_prefix = os.path.join(self._dir, 'expr_modified')
        lib_fn = os.path.join(self._dir, 'mapping.lib')
        #annotation_fn = fse_assets.source.annotation_gtf
        assets = A(A.Source(fse_assets.target.parameters_and_pro),
                   A.Target(rnaseq.FluxSimulatorParametersAndPro(paramsandpro_prefix),
                            rnaseq.FASTQPossiblyGzipCompressed(read1_fn),
                            rnaseq.FASTQPossiblyGzipCompressed(read2_fn),
                            core.File(lib_fn)))
        cmd, returncode = fss.run(assets, parameters=('--read-number', '1000'))
        self.assertEquals(0, returncode)
        # check that the target files are not empty
        self.assertGreater(os.stat(assets.target.parameters_and_pro.name + '.pro').st_size, 0)
        self.assertGreater(os.stat(read1_fn).st_size, 0)
        self.assertGreater(os.stat(read2_fn).st_size, 0)
        self.assertGreater(os.stat(lib_fn).st_size, 0)


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
        self.assertEquals('123', undefoo.bar.name)

    def test_GzipFastqFilePair(self):
        NFRAGMENTS_MATCH = 300
        read1_fh = tempfile.NamedTemporaryFile(prefix='read1', suffix='.fq.gz', dir=self.tempdir, delete=False)
        read1_fh.close()
        read2_fh = tempfile.NamedTemporaryFile(prefix='read2', suffix='.fq.gz', dir=self.tempdir, delete=False)
        read2_fh.close()
        with open(PHAGEFASTA) as fasta_fh:
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


class ModelUtilsTestCase(unittest.TestCase):

    def setUp(self):
        cls_index = rnaseq.Bowtie2Build
        executable_index = 'bowtie2-build'
        cls_align = rnaseq.Bowtie2
        executable_align = 'bowtie2'
        read1_fh = tempfile.NamedTemporaryFile(prefix='read1', suffix='.fq')
        read2_fh = tempfile.NamedTemporaryFile(prefix='read2', suffix='.fq')
        with open(PHAGEFASTA) as fasta_fh:
            reference = next(railroadtracks.model.simulate.readfasta_iter(fasta_fh))
        self._read1_fh, self._read2_fh = railroadtracks.model.simulate.randomPEreads(read1_fh, read2_fh, reference)
        self.tempdir = tempfile.mkdtemp()
        assets_align, cmd, returncode, fh = _build_UpToAlign(cls_index, executable_index,
                                                             cls_align, executable_align, self._read1_fh,
                                                             os.path.join(self.tempdir, 'reference'))
        self._assets_align = assets_align
        self._fh = fh # protect temp file from deletion

    def tearDown(self):
        self._read1_fh.close()
        self._read2_fh.close()
        shutil.rmtree(self.tempdir)

    @unittest.skipIf(not environment.Executable.ispresent('bowtie2-build'),
                     'bowtie2-build is not in the PATH')
    def test_ensuresam(self):
        assets_align = self._assets_align
        fn = assets_align.target.alignment.name
        self.assertTrue(fn.endswith('.bam'))
        fh_sam = railroadtracks.model.files.ensure_sam(fn)
        self.assertTrue(fh_sam.name.endswith('.sam'))
        fh_sam_again = railroadtracks.model.files.ensure_sam(fh_sam.name)
        self.assertEquals(fh_sam.name, fh_sam_again.name)

    @unittest.skipIf(not environment.Executable.ispresent('bowtie2-build'),
                     'bowtie2-build is not in the PATH')
    def test_ensurebam(self):
        assets_align = self._assets_align
        fn = assets_align.target.alignment.name
        self.assertTrue(fn.endswith('.bam'))
        fh_sam = railroadtracks.model.files.ensure_sam(fn)
        fh_bam = railroadtracks.model.files.ensure_bam(fh_sam.name)
        fh_bam_again = railroadtracks.model.files.ensure_bam(fh_bam.name)
        self.assertEquals(fh_bam.name, fh_bam_again.name)
        self.assertTrue(fh_bam.name.endswith('.bam'))
        # this is round trip for the SAM file
        with open(fn) as fh_bam_orig:
            self.assertEqual(fh_bam_orig.read(), fh_bam.read())

    @unittest.skipIf(not environment.Executable.ispresent('bowtie2-build'),
                     'bowtie2-build is not in the PATH')
    def test_SamtoolsBamToSam(self):
        executable = 'samtools'
        bamtosam = railroadtracks.model.files.SamtoolsBamToSam(executable)
        version = bamtosam.version # not a test - just here to check that the accessor is working
        Assets = bamtosam.Assets
        assets_align = self._assets_align
        source = Assets.Source(assets_align.target.alignment)
        sam_fh = tempfile.NamedTemporaryFile(suffix='.sam')
        target = Assets.Target(rnaseq.SavedSAM(sam_fh.name))
        assets = Assets(source, target)
        cmd, returncode = bamtosam.run(assets)
        self.assertEqual(0, returncode)
        # check the infamous '.bam.bam' issue
        self.assertFalse(os.path.exists(sam_fh.name + '.sam')) # if this fails, a newer version of samtools is required
        # check that the target file is not empty
        self.assertGreater(os.stat(sam_fh.name).st_size, 0)

    @unittest.skipIf(not environment.Executable.ispresent('bowtie2-build'),
                     'bowtie2-build is not in the PATH')
    def test_SamtoolsSamToBam(self):
        executable = 'samtools'
        bamtosam = railroadtracks.model.files.SamtoolsBamToSam(executable)
        Assets = bamtosam.Assets
        assets_align = self._assets_align
        sam_fh = tempfile.NamedTemporaryFile(suffix='.sam')
        assets_bamtosam = Assets(Assets.Source(assets_align.target.alignment), 
                                 Assets.Target(rnaseq.SavedSAM(sam_fh.name)))
        cmd, returncode = bamtosam.run(assets_bamtosam)

        samtobam = railroadtracks.model.files.SamtoolsSamToBam(executable)
        version = samtobam.version # not a test - just here to check that the accessor is working
        Assets = samtobam.Assets
        bam_fh = tempfile.NamedTemporaryFile(suffix='.bam')
        assets_samtobam = Assets(Assets.Source(assets_bamtosam.target.samfile), 
                                 Assets.Target(rnaseq.BAMFile(bam_fh.name)))
        cmd, returncode = samtobam.run(assets_samtobam)
        self.assertEqual(0, returncode)
        # check that the target file is not empty
        self.assertGreater(os.stat(bam_fh.name).st_size, 0)


class UnifexTestCase(unittest.TestCase):
    def test__extract_argdict(self):
        # invalid one arg
        args = ('/bar/baz.fq', )
        self.assertRaises(ValueError, unifex._extract_argdict, args)

        # one arg
        arglist = ('foo=/bar/baz.fq', )
        d = unifex._extract_argdict(arglist)
        self.assertEqual(set(('foo', )), set(d.keys()))

        # two args
        arglist = ('foo=/bar/baz.fq', 'bar=/baz/foo.bam')
        d = unifex._extract_argdict(arglist)
        self.assertEqual(set(('foo', 'bar')), set(d.keys()))

        # three args
        arglist = (('a=123', 'b=abc', 'a=456'))
        argdict = unifex._extract_argdict(arglist)
        # 
        self.assertEquals(set(argdict.keys()), set(('a','b')))
        self.assertEquals(argdict['a'], ['123','456'])
        self.assertEquals(argdict['b'], ['abc',])
        # test round trip
        arglist2 = unifex._extract_arglist(argdict)
        self.assertEqual(list(arglist).sort(), list(arglist2).sort())

    def test_build_AssetSet(self):
        values = (('abc',), ('123','456'))
        AssetSet = core.assetfactory('AssetSet', [core.AssetAttr('afile', core.File,''),
                                                  core.AssetAttr('asequenceoffiles', core.FileSequence, '')])
        assetset = unifex.build_AssetSet(AssetSet, values)
        #FIXME: where is the test ?
        #FIXME: how are allownone=True assets handled ?

        # values not allowed to be None cannot be None
        values = (None, ('123','456'))
        self.assertRaises(AssertionError, unifex.build_AssetSet, AssetSet, values)

        values = (('abc',), ('123','456'))
        AssetSet = core.assetfactory('AssetSet', [core.AssetAttr('afile', core.File,'', allownone=True),
                                                  core.AssetAttr('asequenceoffiles', core.FileSequence, '')])
        assetset = unifex.build_AssetSet(AssetSet, values)
        values = (None, ('123','456'))
        assetset = unifex.build_AssetSet(AssetSet, values)


    def test_AssetSet_repr(self):
        values = (('abc',), ('123','456'))
        AssetSet = core.assetfactory('AssetSet', [core.AssetAttr('afile', core.File,''),
                                                  core.AssetAttr('asequenceoffiles', core.FileSequence, '')])
        assetset = unifex.build_AssetSet(AssetSet, values)
        res = str(assetset)
        self.assertTrue(isinstance(res, str))


class ModelTestCase(unittest.TestCase):

    def test__make_stepdict(self):
        class MyStep(core.StepAbstract):
            _name = 'My Step'
            _default_execpath = None
            activities = (core.DEFAULT_ACTIVITY.MISC, )
            Assets = None
            version = '0.1'
            def run(self):
                pass
        class MyModule(object):
            pass
        m = MyModule()
        m.Mystep = MyStep
        sl = unifex._make_stepdict(m)
        self.assertEqual(1, len(sl))
        class MyOtherStep(core.StepAbstract):
            _name = 'My Other Step'
            _default_execpath = None
            activities = (core.DEFAULT_ACTIVITY.MISC, )
            Assets = None
            version = '0.1'
            def run(self):
                pass

        m = MyModule()
        m.Mystep = MyStep
        m.MyOtherStep = MyOtherStep
        sl = unifex._make_stepdict(m)
        self.assertEqual(2, len(sl))
        class MyYetAnotherStep(core.StepAbstract):
            _name = 'My Other Step'
            _default_execpath = None
            activities = (core.DEFAULT_ACTIVITY.MISC, )
            Assets = None
            version = '0.1'
            def run(self):
                pass
        m = MyModule()
        m.Mystep = MyStep
        m.MyOtherStep = MyOtherStep
        m.MyYetAnotherStep = MyYetAnotherStep
        self.assertRaises(ValueError, unifex._make_stepdict, m)

    # note: testing core.StepAbstract
    def test_newStepAbstract(self):
        # dummy step that uses the UNIX command 'cat' to concatenate files
        class CatStep(core.StepAbstract):
            _name = 'cat'
            _default_execpath = 'cat'
            Assets = None
            activities = (core.DEFAULT_ACTIVITY.MISC, )
            def getversion(self):
                # dummy version that is the time stamp for the executable
                path = subprocess.check_output(('which', 'cat'))
                ctime = os.path.ctime(path)
                return str(ctime)
            version = property(getversion)
            def run(self, sources, targets, parameters):
                output_filename = targets['output'][0]
                with open(output_filename, 'w') as output:
                    args = ['cat', ]
                    args.extend(sources['input'])
                    res = subprocess.check_call(args, stdout=output)
                return res

        cs = CatStep()
        inputs = (tempfile.NamedTemporaryFile(), tempfile.NamedTemporaryFile())
        inputs[0].write('abc')
        inputs[0].flush()
        inputs[1].write('def')
        inputs[1].flush()
        outputs = (tempfile.NamedTemporaryFile(),)
        sources = {'input': tuple(x.name for x in inputs, )}
        targets = {'output': tuple(x.name for x in outputs)}
        res = cs.run(sources, targets, '')
        out_res = outputs[0].readlines()
        self.assertEqual(['abcdef'], out_res)



def _build_StepIndex(cls, executable, 
                     tempdir, parameters=()):
    """ 
    :param cls: class used to model the step
    :param executable: name (or path) of executable 
    :param tempdir: temporary directory 

    :return (assets, cmd, returncode):"""
    runner = cls(executable)
    AssetsIndexer = railroadtracks.model.aligners.AssetsIndexer
    referencegenome = rnaseq.SavedFASTA(PHAGEFASTA)
    source = AssetsIndexer.Source(referencegenome)
    indexthing = rnaseq.FilePattern(tempdir)
    target = AssetsIndexer.Target(indexthing)
    assets = railroadtracks.model.aligners.AssetsIndexer(source,
                                                         target)
    cmd, returncode = runner.run(assets, parameters = parameters)
    return (assets, cmd, returncode)


def _build_StepAlign(assets_index,
                     cls_align, executable_align,
                     read1_fh):
    runner = cls_align(executable_align)
    AssetsAligner = railroadtracks.model.aligners.AssetsAligner
    # single reads
    read1 = rnaseq.FASTQPossiblyGzipCompressed(read1_fh.name)
    source = AssetsAligner.Source(assets_index.target.indexfilepattern,
                                  read1, None)
    field_i = AssetsAligner.Target._fields.index('alignment')
    cls_alignedreads = getattr(AssetsAligner.Target, core.AssetMetaReserved.SOURCES.value)[field_i].cls
    if cls_alignedreads is rnaseq.SavedSAM:
        suffix='.sam'
    elif cls_alignedreads is rnaseq.BAMFile:
        suffix='.bam'
    else:
        raise ValueError('Invalid class for alignments')
    alignment_fh = tempfile.NamedTemporaryFile(suffix=suffix)
    alignment_fh.close()
    alignthing = cls_alignedreads(alignment_fh.name)
    target = AssetsAligner.Target(alignthing)
    assets = AssetsAligner(source,
                           target)
    cmd, returncode = runner.run(assets)
    return (runner.version, 
            assets, cmd, returncode, alignment_fh)

def _build_UpToAlign(cls_index, executable_index,
                     cls_align, executable_align,
                     read1_fh, tempdir):
    assets_index, cmd, returncode = _build_StepIndex(cls_index, executable_index,
                                                     tempdir)
    version, assets, cmd, returncode, fh = _build_StepAlign(assets_index,
                                                            cls_align, executable_align,
                                                            read1_fh)
    return (assets, cmd, returncode, fh)

def _build_StepQuantify(assets_align,
                        cls_count, executable_count,
                        parameters_count):
    #FIXME: should be sorted by ID !
    runner = cls_count(executable_count)
    AssetsQuantifier = rnaseq.AssetsQuantifier
    alignmentfiles = tuple(assets_align.target.alignment)
    #FIXME: len == 1 (tested elsewhere)
    f_align = rnaseq.BAMFile(os.path.join(os.path.dirname(assets_align.target.alignment.name),
                                          alignmentfiles[0]))
    referenceannotation = rnaseq.SavedGFF(PHAGEGFF)
    source = AssetsQuantifier.Source(f_align,
                                     referenceannotation)
    fh = tempfile.NamedTemporaryFile(suffix='.csv')
    counts = rnaseq.SavedCSV(fh.name)
    target = AssetsQuantifier.Target(counts)
    assets = AssetsQuantifier(source, target)
    cmd, returncode = runner.run(assets, parameters = parameters_count)
    return (assets, cmd, returncode, fh) # fh returned to be protected from premature destruction
    

class ModelIndexTestCase(unittest.TestCase):

    def setUp(self):
        self.tempdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tempdir)

    @unittest.skipIf(not environment.Executable.ispresent('bowtie2-build'),
                     'bowtie2-build is not in the PATH')
    def test_StepBowtie2Build(self):
        execname = 'bowtie2-build'
        runner = rnaseq.Bowtie2Build(execname)
        #FIXME: test version
        version = runner.version
        self.assertTrue(isinstance(version, str))
        self.assertEqual(set(runner.activities), set((rnaseq.ACTIVITY.INDEX, )))
        assets, cmd, returncode = _build_StepIndex(rnaseq.Bowtie2Build, execname,
                                                   os.path.join(self.tempdir, 'reference'))
        self.assertEqual(0, returncode)

    @unittest.skipIf(not environment.Executable.ispresent('bowtie-build'),
                     'bowtie-build is not in the PATH')
    def test_StepBowtieBuild(self):
        execname = 'bowtie-build'
        runner = rnaseq.BowtieBuild(execname)
        #FIXME: test version
        version = runner.version
        self.assertTrue(isinstance(version, str))
        self.assertEqual(set(runner.activities), set((rnaseq.ACTIVITY.INDEX, )))
        assets, cmd, returncode = _build_StepIndex(rnaseq.BowtieBuild, execname,
                                                   os.path.join(self.tempdir, 'reference'))
        self.assertEqual(0, returncode)
        # check that files actually present
        files = tuple(assets.target.indexfilepattern.iterlistfiles())
        self.assertGreater(len(files), 0)

    @unittest.skipIf(not environment.Executable.ispresent('bwa'),
                     'bwa is not in the PATH')
    def test_StepBWAIndex(self):
        execname = 'bwa'
        runner = rnaseq.BWAIndex(execname)
        #FIXME: test version
        version = runner.version
        self.assertTrue(isinstance(version, str))
        self.assertEqual(set(runner.activities), set((rnaseq.ACTIVITY.INDEX, )))
        assets, cmd, returncode = _build_StepIndex(rnaseq.BWAIndex, execname,
                                                   os.path.join(self.tempdir, 'reference'))
        self.assertEqual(0, returncode)
        # check that files actually present
        files = tuple(assets.target.indexfilepattern.iterlistfiles())
        self.assertGreater(len(files), 0)

    @unittest.skipIf(not environment.Executable.ispresent('STAR'),
                     'STAR is not in the PATH')
    def test_StepStarIndex(self):
        execname = 'STAR'
        runner = rnaseq.StarIndex(execname)
        #FIXME: test version
        version = runner.version
        self.assertTrue(isinstance(version, str))
        self.assertEqual(set(runner.activities), set((rnaseq.ACTIVITY.INDEX, )))
        assets, cmd, returncode = _build_StepIndex(rnaseq.StarIndex, execname, 
                                                   os.path.join(self.tempdir, 'reference'),
                                                   parameters=('--genomeChrBinNbits', '12'))
        self.assertEqual(0, returncode)

    @unittest.skipIf(not environment.Executable.ispresent('gmap_build'),
                     'gmap_build is not in the PATH')
    def test_StepGsnapIndex(self):
        execname = 'gmap_build'
        runner = rnaseq.GsnapIndex(execname)
        #FIXME: test version
        version = runner.version
        self.assertTrue(isinstance(version, str))
        self.assertEqual(set(runner.activities), set((rnaseq.ACTIVITY.INDEX, )))
        assets, cmd, returncode = _build_StepIndex(rnaseq.GsnapIndex, execname,
                                                   os.path.join(self.tempdir, 'reference'))
        self.assertEqual(0, returncode)

    @unittest.skipIf(not environment.Executable.ispresent('sailfish'),
                     'sailfish is not in the PATH')
    def test_StepSailfishIndex(self):
        execname = 'sailfish'
        runner = rnaseq.SailfishIndex(execname)
        #FIXME: test version
        version = runner.version
        self.assertTrue(isinstance(version, str))
        self.assertEqual(set(runner.activities), set((rnaseq.ACTIVITY.INDEX, )))
        assets, cmd, returncode = _build_StepIndex(rnaseq.SailfishIndex, execname, 
                                                   os.path.join(self.tempdir, 'reference'),
                                                   parameters = rnaseq.SailfishIndex.PARAMETERS_DEFAULT)
        self.assertEqual(0, returncode)


class ModelAlignTestCase(unittest.TestCase):

    def setUp(self):
        NFRAGMENTS_NOMATCH = 10
        read1_fh = tempfile.NamedTemporaryFile(prefix='read1', suffix='.fq')
        read2_fh = tempfile.NamedTemporaryFile(prefix='read2', suffix='.fq')
        with open(PHAGEFASTA) as fasta_fh:
            reference = next(railroadtracks.model.simulate.readfasta_iter(fasta_fh))
        self._read1_fh, self._read2_fh = railroadtracks.model.simulate.randomPEreads(read1_fh, read2_fh, reference)
        # reads from a random reference (should not match)
        randint = railroadtracks.model.simulate.random.randint
        rand_reference = railroadtracks.model.simulate.Entry('> Random DNA',
                                                             ''.join('ATCG'[randint(0, 3)] for x in range(500)))
        self._read1_fh, self._read2_fh = railroadtracks.model.simulate.randomPEreads(read1_fh, read2_fh,
                                                                                     rand_reference,
                                                                                     n = NFRAGMENTS_NOMATCH)
        self.tempdir = tempfile.mkdtemp()

    def tearDown(self):
        self._read1_fh.close()
        self._read2_fh.close()
        shutil.rmtree(self.tempdir)

    def _test_StepAlign_noreads(self, 
                                cls_index, executable_index,
                                cls_align, executable_align):
        assets_index, cmd, index_res = _build_StepIndex(cls_index, executable_index,
                                                        os.path.join(self.tempdir, 'reference'))
        runner = cls_align(executable_align)
        self.assertEqual(set(runner.activities), set((rnaseq.ACTIVITY.ALIGN, )))
        AssetsAligner = railroadtracks.model.aligners.AssetsAligner
        # no reads
        self.assertRaises(AssertionError,
                          AssetsAligner.Source,
                          assets_index.target.indexfilepattern)

    def _test_StepAlign_singlereads(self, 
                                    cls_index, executable_index,
                                    cls_align, executable_align):
        runner = cls_align(executable_align)
        self.assertEqual(set(runner.activities), set((rnaseq.ACTIVITY.ALIGN, )))
        # single reads
        assets, cmd, returncode, fh = _build_UpToAlign(cls_index, executable_index,
                                                       cls_align, executable_align, self._read1_fh,
                                                       os.path.join(self.tempdir, 'reference'))
        self.assertTrue(isinstance(runner.version, str))
        self.assertEquals(0, returncode)
        # FIXME: check that the alignment file contains what it should
        self.assertTrue(os.path.exists(fh.name))
        # check that the target file is not empty
        self.assertGreater(os.stat(fh.name).st_size, 0)

    def _test_StepAlign_pairedreads(self, 
                                    cls_index, executable_index,
                                    cls_align, executable_align):
        assets_index, cmd, returncode = _build_StepIndex(cls_index, executable_index,
                                                         os.path.join(self.tempdir, 'reference'))
        runner = cls_align(executable_align)
        self.assertEqual(set(runner.activities), set((rnaseq.ACTIVITY.ALIGN, )))
        AssetsAligner = railroadtracks.model.aligners.AssetsAligner
        # paired-end reads
        read1 = rnaseq.FASTQPossiblyGzipCompressed(self._read1_fh.name)
        read2 = rnaseq.FASTQPossiblyGzipCompressed(self._read2_fh.name)
        source = AssetsAligner.Source(assets_index.target.indexfilepattern,
                                      read1, read2)
        field_i = AssetsAligner.Target._fields.index('alignment')
        cls_alignedreads = getattr(AssetsAligner.Target, core.AssetMetaReserved.SOURCES.value)[field_i].cls
        if cls_alignedreads is rnaseq.SavedSAM:
            suffix='.sam'
        elif cls_alignedreads is rnaseq.BAMFile:
            suffix='.bam'
        else:
            raise ValueError('Invalid class for alignments')
        fh = tempfile.NamedTemporaryFile(suffix=suffix)
        fh.close()
        try:
            alignthing = cls_alignedreads(fh.name)
            target = AssetsAligner.Target(alignthing)
            assets = AssetsAligner(source,
                                   target)
            cmd, returncode = runner.run(assets)
            # FIXME: check that the alignment file contains what it should
            self.assertEqual(0, returncode)
            self.assertTrue(os.path.exists(fh.name))
            # check that the target file is not empty
            self.assertGreater(os.stat(fh.name).st_size, 0)
        finally:
            if os.path.exists(fh.name):
                os.unlink(fh.name)


    @unittest.skipIf(not environment.Executable.ispresent('bwa'),
                     'bwa is not in the PATH')
    def test_StepBWA_noread(self):
        self._test_StepAlign_noreads(rnaseq.BWAIndex, 'bwa',
                                     rnaseq.BWA, 'bwa')
    @unittest.skipIf(not environment.Executable.ispresent('bwa'),
                     'bwa is not in the PATH')
    def test_StepBWA_singlereads(self):
        self._test_StepAlign_singlereads(rnaseq.BWAIndex, 'bwa',
                                         rnaseq.BWA, 'bwa')
    @unittest.skipIf(not environment.Executable.ispresent('bwa'),
                     'bwa is not in the PATH')
    def test_StepBWA_pairedreads(self):
        self._test_StepAlign_pairedreads(rnaseq.BWAIndex, 'bwa',
                                         rnaseq.BWA, 'bwa')

       
    @unittest.skipIf(not environment.Executable.ispresent('bowtie'),
                     'bowtie is not in the PATH')
    def test_StepBowtie_noread(self):
        self._test_StepAlign_noreads(rnaseq.Bowtie2Build, 'bowtie-build',
                                     rnaseq.Bowtie2, 'bowtie')
    @unittest.skipIf(not environment.Executable.ispresent('bowtie'),
                     'bowtie is not in the PATH')
    def test_StepBowtie_singlereads(self):
        self._test_StepAlign_singlereads(rnaseq.BowtieBuild, 'bowtie-build',
                                         rnaseq.Bowtie, 'bowtie')
    @unittest.skipIf(not environment.Executable.ispresent('bowtie'),
                     'bowtie is not in the PATH')
    def test_StepBowtie_pairedreads(self):
        self._test_StepAlign_pairedreads(rnaseq.BowtieBuild, 'bowtie-build',
                                         rnaseq.Bowtie, 'bowtie')

    @unittest.skipIf(not environment.Executable.ispresent('bowtie2'),
                     'bowtie2 is not in the PATH')
    def test_StepBowtie2_noread(self):
        self._test_StepAlign_noreads(rnaseq.Bowtie2Build, 'bowtie2-build',
                                     rnaseq.Bowtie2, 'bowtie2')
    @unittest.skipIf(not environment.Executable.ispresent('bowtie2'),
                     'bowtie2 is not in the PATH')
    def test_StepBowtie2_singlereads(self):
        self._test_StepAlign_singlereads(rnaseq.Bowtie2Build, 'bowtie2-build',
                                         rnaseq.Bowtie2, 'bowtie2')
    @unittest.skipIf(not environment.Executable.ispresent('bowtie2'),
                     'bowtie2 is not in the PATH')
    def test_StepBowtie2_pairedreads(self):
        self._test_StepAlign_pairedreads(rnaseq.Bowtie2Build, 'bowtie2-build',
                                         rnaseq.Bowtie2, 'bowtie2')


    @unittest.skipIf(not environment.Executable.ispresent('tophat'),
                     'tophat is not in the PATH')
    def test_StepTopHat_noread(self):
        self._test_StepAlign_noreads(rnaseq.Bowtie2Build, 'bowtie2-build',
                                     rnaseq.TopHat, 'tophat')
    @unittest.skipIf(not environment.Executable.ispresent('tophat'),
                     'tophat is not in the PATH')
    def test_StepTopHat_singlereads(self):
        self._test_StepAlign_singlereads(rnaseq.Bowtie2Build, 'bowtie2-build',
                                         rnaseq.TopHat, 'tophat')
    @unittest.skipIf(not environment.Executable.ispresent('tophat'),
                     'tophat is not in the PATH')
    def test_StepTopHat_pairedreads(self):
        self._test_StepAlign_pairedreads(rnaseq.BowtieBuild, 'bowtie2-build',
                                         rnaseq.TopHat, 'tophat')

    @unittest.skipIf(not environment.Executable.ispresent('tophat2'),
                     'tophat is not in the PATH')
    def test_StepTopHat2_noread(self):
        self._test_StepAlign_noreads(rnaseq.Bowtie2Build, 'bowtie2-build',
                                     rnaseq.TopHat2, 'tophat2')
    @unittest.skipIf(not environment.Executable.ispresent('tophat2'),
                     'tophat2 is not in the PATH')
    def test_StepTopHat2_singlereads(self):
        self._test_StepAlign_singlereads(rnaseq.Bowtie2Build, 'bowtie2-build',
                                         rnaseq.TopHat2, 'tophat2')
    @unittest.skipIf(not environment.Executable.ispresent('tophat2'),
                     'tophat2 is not in the PATH')
    def test_StepTopHat2_pairedreads(self):
        self._test_StepAlign_pairedreads(rnaseq.BowtieBuild, 'bowtie2-build',
                                         rnaseq.TopHat, 'tophat2')

    @unittest.skipIf(not environment.Executable.ispresent('gsnap'),
                     'gsnap is not in the PATH')
    @unittest.skipIf(not environment.Executable.ispresent('gmap_build'),
                     'gmap_build is not in the PATH')
    def test_StepGsnap_singlereads(self):
        self._test_StepAlign_singlereads(rnaseq.GsnapIndex, 'gmap_build',
                                         rnaseq.GsnapAlign, 'gsnap')
    @unittest.skipIf(not environment.Executable.ispresent('gsnap'),
                     'gsnap is not in the PATH')
    @unittest.skipIf(not environment.Executable.ispresent('gmap_build'),
                     'gmap_build is not in the PATH')
    def test_StepGsnapAlign_pairedreads(self):
        self._test_StepAlign_pairedreads(rnaseq.GsnapIndex, 'gmap_build',
                                         rnaseq.GsnapAlign, 'gsnap')

    @unittest.skipIf(not environment.Executable.ispresent('STAR'),
                     'STAR is not in the PATH')
    def test_StepStarAlign_singlereads(self):
        self._test_StepAlign_singlereads(rnaseq.StarIndex, 'STAR',
                                         rnaseq.StarAlign, 'STAR')
    @unittest.skipIf(not environment.Executable.ispresent('STAR'),
                     'STAR is not in the PATH')
    def test_StepStarAlign_pairedreads(self):
        self._test_StepAlign_pairedreads(rnaseq.StarIndex, 'STAR',
                                         rnaseq.StarAlign, 'STAR')


class ModelHandleBAMTestCase(unittest.TestCase):
    """
    Tests for the handling of BAM files (Steps and associated).
    """

    def setUp(self):
        NFRAGMENTS_NOMATCH = 10
        cls_index = rnaseq.Bowtie2Build
        executable_index = 'bowtie2-build'
        cls_align = rnaseq.Bowtie2
        executable_align = 'bowtie2'
        executable_samtools = 'samtools'
        read1_fh = tempfile.NamedTemporaryFile(prefix='read1', suffix='.fq')
        read2_fh = tempfile.NamedTemporaryFile(prefix='read2', suffix='.fq')
        with open(PHAGEFASTA) as fasta_fh:
            reference = next(railroadtracks.model.simulate.readfasta_iter(fasta_fh))
        self._read1_fh, self._read2_fh = railroadtracks.model.simulate.randomPEreads(read1_fh, read2_fh, reference)

        # reads from a random reference (should not match)
        randint = railroadtracks.model.simulate.random.randint
        rand_reference = railroadtracks.model.simulate.Entry('> Random DNA',
                                                             ''.join('ATCG'[randint(0, 3)] for x in range(500)))
        self._read1_fh, self._read2_fh = railroadtracks.model.simulate.randomPEreads(read1_fh, read2_fh,
                                                                                     rand_reference,
                                                                                     n = NFRAGMENTS_NOMATCH)

        self.tempdir = tempfile.mkdtemp()
        assets_align, cmd, returncode, fh = _build_UpToAlign(cls_index, executable_index,
                                                             cls_align, executable_align, self._read1_fh,
                                                             os.path.join(self.tempdir, 'reference'))
        self._assets_align = assets_align

    def tearDown(self):
        self._read1_fh.close()
        self._read2_fh.close()
        shutil.rmtree(self.tempdir)


    @unittest.skipIf(not (environment.Executable.ispresent('samtools') and \
                          environment.Executable.ispresent('bowtie2-build')),
                     'samtools and bowtie2-build should be in the PATH')        
    def test_SamtoolsSortByID(self):
        AssetsSorter = rnaseq.AssetsSorter
        SamtoolsSorterByID = rnaseq.SamtoolsSorterByID
        executable_sort = 'samtools'
        samsort_byID = SamtoolsSorterByID(executable_sort)
        self.assertEqual(set((rnaseq.ACTIVITY.SORT,)), set(samsort_byID.activities))
        # mostly to check that version is working
        self.assertTrue(isinstance(samsort_byID.version, str))
        source = AssetsSorter.Source(self._assets_align.target.alignment)
        fh = tempfile.NamedTemporaryFile(suffix='.bam')
        sortedreads = rnaseq.BAMFile(fh.name)
        target = AssetsSorter.Target(sortedreads)
        assets = AssetsSorter(source, target)
        cmd, returncode = samsort_byID.run(assets)
        self.assertEqual(0, returncode)
        # check that the sorted BAM is not an empty file
        self.assertGreater(os.stat(fh.name).st_size, 0)

    @unittest.skipIf(not (environment.Executable.ispresent('samtools') and \
                          environment.Executable.ispresent('bowtie2-build')),
                     'samtools and bowtie2-build should be in the PATH')        
    def test_SamtoolsExtractUnaligned(self):
        SamtoolsExtractUnaligned = rnaseq.SamtoolsExtractUnaligned
        executable_sort = 'samtools'
        extractunaligned = rnaseq.SamtoolsExtractUnaligned(executable_sort)
        self.assertEqual(set((rnaseq.ACTIVITY.UTILITY,)), set(extractunaligned.activities))
        # mostly to check that version is working
        self.assertTrue(isinstance(extractunaligned.version, str))
        Assets = SamtoolsExtractUnaligned.Assets
        source = Assets.Source(self._assets_align.target.alignment)
        fh = tempfile.NamedTemporaryFile(suffix='.bam')
        sortedreads = rnaseq.BAMFile(fh.name)
        target = Assets.Target(sortedreads)
        assets = Assets(source, target)
        cmd, returncode = extractunaligned.run(assets)
        self.assertEqual(0, returncode)
        #FIXME: other tests on the results


    @unittest.skipIf(not (environment.Executable.ispresent('samtools') and \
                          environment.Executable.ispresent('bowtie2-build')),
                     'samtools and bowtie2-build should be in the PATH')        
    def test_SamtoolsFilter(self):
        SamtoolsFilter = railroadtracks.model.files.SamtoolsFilter
        executable_sort = 'samtools'
        readfilter = SamtoolsFilter(executable_sort)
        self.assertEqual(set((railroadtracks.model.files.ACTIVITY.FILTERREADS,)), set(readfilter.activities))
        # mostly to check that version is working
        self.assertTrue(isinstance(readfilter.version, str))
        Assets = SamtoolsFilter.Assets
        source = Assets.Source(self._assets_align.target.alignment)
        # 
        fh_f0x4 = tempfile.NamedTemporaryFile(suffix='.bam')
        sortedreads = rnaseq.BAMFile(fh_f0x4.name)
        target = Assets.Target(sortedreads)
        assets = Assets(source, target)
        cmd, returncode = readfilter.run(assets, ('--filter-include', '0x0004'))
        self.assertEqual(0, returncode)
        in_stat = os.stat(assets.source.bamfile.name)
        out_stat_inc = os.stat(assets.target.bamfile.name)
        self.assertGreater(in_stat.st_size, out_stat_inc.st_size)
        #
        fh_F0x4 = tempfile.NamedTemporaryFile(suffix='.bam')
        sortedreads = rnaseq.BAMFile(fh_F0x4.name)
        target = Assets.Target(sortedreads)
        assets = Assets(source, target)
        cmd, returncode = readfilter.run(assets, ('--filter-exclude', '0x0004'))
        self.assertEqual(0, returncode)
        out_stat_exc = os.stat(assets.target.bamfile.name)
        self.assertGreater(out_stat_exc.st_size, out_stat_inc.st_size)


class ModelAlignQuantificationTestCase(unittest.TestCase):
    """ At the time of writing, this is for Sailfish """
    def setUp(self):
        self.tempdir = tempfile.mkdtemp()
        read1_fh = tempfile.NamedTemporaryFile(prefix='read1', suffix='.fq', dir=self.tempdir)
        read2_fh = tempfile.NamedTemporaryFile(prefix='read2', suffix='.fq', dir=self.tempdir)
        with open(PHAGEFASTA) as fasta_fh:
            reference = next(railroadtracks.model.simulate.readfasta_iter(fasta_fh))
        self._read1_fh, self._read2_fh = railroadtracks.model.simulate.randomPEreads(read1_fh, read2_fh, reference)

    def tearDown(self):
        shutil.rmtree(self.tempdir)

    @unittest.skipIf(not environment.Executable.ispresent('sailfish'),
                     'sailfish is not in the PATH')        
    def test_pairedend(self):
        execname = 'sailfish'
        runner = rnaseq.SailfishIndex(execname)
        assets, cmd, returncode = _build_StepIndex(rnaseq.SailfishIndex, execname, 
                                                   os.path.join(self.tempdir, 'reference'),
                                                   parameters = runner.PARAMETERS_DEFAULT)
        if returncode != 0:
            raise Exception("Index not built")
        self.assets = assets
        sf = rnaseq.SailfishQuant(execname)
        A = sf.Assets
        count_csv = os.path.join(self.tempdir, 'counts.csv')
        output_dir = os.path.join(self.tempdir, 'sailfish_output_dir')
        assets = A(A.Source(self.assets.target.indexfilepattern,
                            rnaseq.FASTQPossiblyGzipCompressed(self._read1_fh.name),
                            rnaseq.FASTQPossiblyGzipCompressed(self._read2_fh.name)),
                   A.Target(rnaseq.SavedCSV(count_csv),
                            rnaseq.FilePattern(output_dir)))
        cmd, returncode = sf.run(assets, parameters = ('--libtype', sf.LIBRARY_PE))

class ModelQuantificationTestCase(unittest.TestCase):
    """
    """

    def setUp(self):
        cls_index = rnaseq.Bowtie2Build
        executable_index = 'bowtie2-build'
        cls_align = rnaseq.Bowtie2
        executable_align = 'bowtie2'
        self.tempdir = tempfile.mkdtemp()
        read1_fh = tempfile.NamedTemporaryFile(prefix='read1', suffix='.fq')
        read2_fh = tempfile.NamedTemporaryFile(prefix='read2', suffix='.fq')
        with open(PHAGEFASTA) as fasta_fh:
            reference = next(railroadtracks.model.simulate.readfasta_iter(fasta_fh))
        self._read1_fh, self._read2_fh = railroadtracks.model.simulate.randomPEreads(read1_fh, read2_fh, reference)
        assets_align, cmd, returncode, fh = _build_UpToAlign(cls_index, executable_index,
                                                             cls_align, executable_align,
                                                             self._read1_fh,
                                                             os.path.join(self.tempdir, 'reference'))
        self._assets_align = assets_align
        self._fh = fh # protect temp file from deletion

    def tearDown(self):
        self._read1_fh.close()
        self._read2_fh.close()
        shutil.rmtree(self.tempdir)

    @unittest.skipIf(not (environment.Executable.ispresent('htseq-count') and \
                          environment.Executable.ispresent('bowtie2-build')),
                     'htseq-count is not in the PATH')        
    def test_HTSeqCount(self):
        AssetsQuantifier = rnaseq.AssetsQuantifier
        HTSeqCount = rnaseq.HTSeqCount
        referenceannotation = rnaseq.SavedGFF(PHAGEGFF)
        executable_count = 'htseq-count'
        # non-default parameters to fit our demo GFF
        htseqcount = HTSeqCount(executable_count)
        self.assertEqual(set((rnaseq.ACTIVITY.QUANTIFY,)), set(htseqcount.activities))

        version = htseqcount.version

        assets_align = self._assets_align
        # sort to remove the warning
        stsort = rnaseq.SamtoolsSorterByID()
        with tempfile.NamedTemporaryFile(suffix='.bam') as sortedbam_fh, \
             tempfile.NamedTemporaryFile(suffix='.csv') as fh:
            assets_sort = stsort.Assets(stsort.Assets.Source(assets_align.target.alignment),
                                        stsort.Assets.Target(railroadtracks.model.files.BAMFileSortedByID(sortedbam_fh.name)))
            stsort.run(assets_sort)
            source = AssetsQuantifier.Source(assets_sort.target.sortedbam,
                                             referenceannotation)
            
            counts = rnaseq.SavedCSV(fh.name)
            target = AssetsQuantifier.Target(counts)
            assets = AssetsQuantifier(source, target)
            cmd, returncode = htseqcount.run(assets,
                                             parameters=HTSeqCount._noexons_parameters)
            self.assertEqual(0, returncode)
            cmd, returncode = htseqcount.run(assets,
                                             parameters=itertools.chain(HTSeqCount._noexons_parameters,
                                                                        ('-m', 'union', '-s', 'no', '-a', '10')))
            self.assertEqual(0, returncode)

    @unittest.skipIf(not (environment.Executable.ispresent('htseq-count') and \
                          environment.Executable.ispresent('R') and \
                          environment.R('R').packageversion_or_none('Rsubread') is not None),
                     'Htseq-count is not in the PATH')        
    def test_FeatureCount(self):
        AssetsQuantifier = rnaseq.AssetsQuantifier
        FeatureCount = rnaseq.FeatureCount
        referenceannotation = rnaseq.SavedGFF(PHAGEGFF)
        # wrapped R script
        executable_count = 'R'
        # non-default parameters to fit our demo GFF
        featurecount = FeatureCount(executable_count)
        self.assertEqual(set((rnaseq.ACTIVITY.QUANTIFY,)), set(featurecount.activities))

        version = featurecount.version

        assets_align = self._assets_align
        source = AssetsQuantifier.Source(assets_align.target.alignment,
                                         referenceannotation)
        fh = tempfile.NamedTemporaryFile(suffix='.csv')
        counts = rnaseq.SavedCSV(fh.name)
        target = AssetsQuantifier.Target(counts)
        assets = AssetsQuantifier(source, target)
        cmd, returncode = featurecount.run(assets, 
                                           parameters=featurecount._noexons_parameters)
        self.assertEqual(0, returncode)
        # try with paired-end (the total of counts should now be approx. half)
        cmd, returncode = featurecount.run(assets, parameters=tuple(itertools.chain(featurecount._noexons_parameters, featurecount._pe_parameters)))
        self.assertEqual(0, returncode)

class ModelColumnMergerTestCase(unittest.TestCase):

    def test_merge(self):
        with tempfile.NamedTemporaryFile(suffix=".csv") as fh_col1, \
             tempfile.NamedTemporaryFile(suffix=".csv") as fh_col2, \
             tempfile.NamedTemporaryFile(suffix=".csv") as fh_merged:
            fh_col1_csv = csv.writer(fh_col1)
            fh_col2_csv = csv.writer(fh_col2)
            for x in (1,2,3):
                fh_col1_csv.writerow([x])
            fh_col1.seek(0)
            for x in (4,5,6):
                fh_col2_csv.writerow([x])
            fh_col2.seek(0)
            ColumnMerger = rnaseq.ColumnMerger
            Assets = ColumnMerger.Assets
            merger = ColumnMerger()
            filestomerge = tuple(rnaseq.SavedCSV(x.name) for x in (fh_col1, fh_col2))
            source = Assets.Source(rnaseq.SavedCSVSequence(filestomerge))
            assets = Assets(source,
                            Assets.Target(rnaseq.SavedCSV(fh_merged.name)))
            cmd, returncode = merger.run(assets,
                                         parameters=('None', '0'))
            self.assertEqual(0, returncode)
            # just in case... rewind to 0
            fh_merged.seek(0)
            self.assertEqual(('1,4',
                              '2,5',
                              '3,6'), tuple(x.rstrip() for x in fh_merged))

    def test_merge_withIDs(self):
        with tempfile.NamedTemporaryFile(suffix=".csv") as fh_col1, \
             tempfile.NamedTemporaryFile(suffix=".csv") as fh_col2, \
             tempfile.NamedTemporaryFile(suffix=".csv") as fh_merged:
            fh_col1_csv = csv.writer(fh_col1)
            fh_col2_csv = csv.writer(fh_col2)
            for row in ((1,'a'),(2,'b'),(3,'c')):
                fh_col1_csv.writerow(row)
            fh_col1.flush()
            fh_col1.seek(0)
            for row in ((1,'d'),(2,'e'),(3,'f')):
                fh_col2_csv.writerow(row)
            fh_col2.flush()
            fh_col2.seek(0)
            ColumnMerger = rnaseq.ColumnMerger
            Assets = ColumnMerger.Assets
            merger = ColumnMerger()
            filestomerge = tuple(rnaseq.SavedCSV(x.name) for x in (fh_col1, fh_col2))
            source = Assets.Source(rnaseq.SavedCSVSequence(filestomerge))
            assets = Assets(source,
                            Assets.Target(rnaseq.SavedCSV(fh_merged.name)))
            cmd, returncode = merger.run(assets, parameters=('0', '1'))
            self.assertEqual(0, returncode)
            # just in case... rewind to 0
            fh_merged.seek(0)
            self.assertEqual(('1,a,d',
                              '2,b,e',
                              '3,c,f'),
                             tuple(x.rstrip() for x in fh_merged))


class ActivitiesSplit(core.Enum):
    HEADTAIL = 'Split in 2'
    MERGE = 'Merge'
    FOO = 'bar'
# Asset for a step that splits a file into head and tail
class AssetsSplit(core.AssetsStep):
    Source = core.assetfactory('Source', [core.AssetAttr('file', core.File, '')])
    Target = core.assetfactory('Target', [core.AssetAttr('head', core.File, ''),
                                          core.AssetAttr('tail', core.File, '')])
# Step that splits a file into head and tail
class Split(core.StepAbstract):
    _name = 'split'
    _default_execpath = None
    Assets=AssetsSplit
    activities = (ActivitiesSplit.HEADTAIL, )
    def __init__(self, executable):
        self._execpath = executable
    @property
    def version(self):
        return '0.1'
    def run(self, assets, parameters=()):
        with open(assets.source.file, 'rb') as fh_in:
            with open(assets.target.head.name, 'wb') as fh_out:
                head = fh_in.read(1)
                fh_out.write(head)
            with open(assets.target.tail.name, 'wb') as fh_out:
                fh_out.write(fh_in)
        cmd = None
        returncode = 0
        return cmd, returncode

class AssetsCat(core.AssetsStep):
    Source = core.assetfactory('Source', [core.AssetAttr('heads', core.FileSequence, '')])
    Target = core.assetfactory('Target', [core.AssetAttr('result', core.File, '')])

class Cat(core.StepAbstract):
    _name = 'cat'
    Assets=AssetsCat
    activities = (ActivitiesSplit.MERGE, )
    _default_execpath = None
    def __init__(self, executable):
        self._execpath = executable
    @property
    def version(self):
        return '0.1'
    def run(self, assets, parameters=()):
        with open(assets.target.result.name, 'wb') as fh_out:
            for x in assets.source.heads:
                with open(x.name, 'rb') as fh_in:
                    fh_out.write(fh_in)
        cmd = None
        returncode = 0
        return cmd, returncode



class ModelDExpressionTestCase(unittest.TestCase):
    """
    """
    
    def setUp(self):
        self.tempdir = tempfile.mkdtemp()
        sampleinfo_fh = tempfile.NamedTemporaryFile(dir=self.tempdir, suffix='.csv', delete=False)
        csv_w = csv.writer(sampleinfo_fh)
        csv_w.writerow(['sample_id', 'group'])
        for i in range(6):
            csv_w.writerow([str(i), ('A','B')[i%2]])
        sampleinfo_fh.flush()
        self._sampleinfo_fh = sampleinfo_fh
        with open(PHAGEFASTA) as fasta_fh:
            reference = next(railroadtracks.model.simulate.readfasta_iter(fasta_fh))
        random.seed(123)    
        with tempfile.NamedTemporaryFile(dir=self.tempdir, suffix='.csv', delete=False) as fh_merged:
            csv_w = csv.writer(fh_merged)
            self.asset_merge = rnaseq.SavedCSV(fh_merged.name)
            self._r_exec = environment.R('R')
            csv_w.writerow(['',] + ['V%i'%x for x in range(1,7)])
            with open(PHAGEGFF) as gff_fh:
                csv_r = csv.reader(gff_fh, delimiter="\t")
                for row in csv_r:
                    if row[0].startswith('#'):
                        continue
                    if row[2] != 'CDS':
                        continue
                    entry_id = re.sub('ID="(.+?)".+', '\\1', row[8])
                    count_row = [entry_id, ] + [int(random.lognormvariate(4, 2)) for x in range(6)]
                    csv_w.writerow(count_row)
            
                    
    def tearDown(self):
        shutil.rmtree(self.tempdir)


    @unittest.skipIf(not (environment.Executable.ispresent('R') and \
                          environment.R('R').packageversion_or_none('DESeq') is not None),
                     'R (with package "DESeq") must be in the PATH')
    def test_StepDESeq(self):
        deseq = rnaseq.DESeq(self._r_exec)
        self.assertEqual(set((rnaseq.ACTIVITY.DIFFEXP,)), set(deseq.activities))
        AssetsDE = railroadtracks.model.diffexp.AssetsDifferentialExpression
        fh = tempfile.NamedTemporaryFile()
        source = AssetsDE.Source(self.asset_merge, 
                                 core.File(self._sampleinfo_fh.name))
        target = AssetsDE.Target(core.File(fh.name))
        assets = AssetsDE(source,
                          target)
        cmd, returncode = deseq.run(assets, parameters=('--dispersion-fittype=local', ))
        self.assertEquals(0, returncode)

    @unittest.skipIf(not (environment.Executable.ispresent('R') and \
                          environment.R('R').packageversion_or_none('DESeq2') is not None),
                     'R (with package "DESeq2") must be in the PATH')
    def test_StepDESeq2(self):
        deseq2 = rnaseq.DESeq2(self._r_exec)
        self.assertEqual(set((rnaseq.ACTIVITY.DIFFEXP,)), set(deseq2.activities))
        AssetsDE = railroadtracks.model.diffexp.AssetsDifferentialExpression
        fh = tempfile.NamedTemporaryFile()
        source = AssetsDE.Source(self.asset_merge,
                                 core.File(self._sampleinfo_fh.name))
        target = AssetsDE.Target(core.File(fh.name))
        assets = AssetsDE(source,
                          target)
        cmd, returncode = deseq2.run(assets)
        self.assertEquals(0, returncode)


    @unittest.skipIf(not (environment.Executable.ispresent('R') and \
                          environment.R('R').packageversion_or_none('DESeq2') is not None),
                     'R (with package "edgeR") must be in the PATH')
    def test_StepEdgeR(self):
        edger = rnaseq.EdgeR(self._r_exec)
        self.assertEqual(set((rnaseq.ACTIVITY.DIFFEXP,)), set(edger.activities))
        AssetsDE = railroadtracks.model.diffexp.AssetsDifferentialExpression
        fh = tempfile.NamedTemporaryFile()
        source = AssetsDE.Source(self.asset_merge,
                                 core.File(self._sampleinfo_fh.name))
        target = AssetsDE.Target(core.File(fh.name))
        assets = AssetsDE(source,
                          target)
        cmd, returncode = edger.run(assets)
        self.assertEquals(0, returncode)


    @unittest.skipIf(not (environment.Executable.ispresent('R') and \
                          environment.R('R').packageversion_or_none('limma') is not None),
                     'R (with package "limma") must be in the PATH')
    def test_StepLimmaVoom(self):
        voom = rnaseq.LimmaVoom(self._r_exec)
        self.assertEqual(set((rnaseq.ACTIVITY.DIFFEXP,)), set(voom.activities))
        AssetsDE = railroadtracks.model.diffexp.AssetsDifferentialExpression
        fh = tempfile.NamedTemporaryFile()
        source = AssetsDE.Source(self.asset_merge,
                                 core.File(self._sampleinfo_fh.name))
        target = AssetsDE.Target(core.File(fh.name))
        assets = AssetsDE(source,
                          target)
        cmd, returncode = voom.run(assets)
        self.assertEquals(0, returncode)


class PersistanceTestCase(unittest.TestCase):

    def setUp(self):
        self.cache_file = tempfile.NamedTemporaryFile()
        class PythonTime(core.StepAbstract):
            """ Dummy step - print the local date/time using Python. """
            _name = 'time'
            _default_execpath = 'time'
            activities = (core.DEFAULT_ACTIVITY.MISC, )
            Assets = None
            def __init__(self, executable=None):
                if executable is None:
                    self._execpath = self._default_execpath
                else:
                    self._execpath = executable
            @property
            def version(self):
                res = subprocess.check_output([self._execpath, '--version'],
                                              stderr=subprocess.STDOUT)
                version = res[:res.find(os.linesep)]
                return version
            def run(self, assets, parameters=('%D', )):
                # assets not used
                res = subprocess.check_output([self._execpath, '-c',
                                               "import %s; print(time.strftime('%D', time.localtime()))" % self._default_execpath])
                print(res)
        self.PythonTime = PythonTime

    def tearDown(self):
        self.cache_file.close()

    def test_PersistentTaskList(self):
        model = core.Model(tuple()) # getting away with an empty model for this step
        # test the initialization
        cache = hortator.PersistentTaskList(self.cache_file.name, model, force_create=True)
        # check that it was created
        self.assertTrue(cache.created)
        # empty set of steps
        s = tuple(cache.iter_steps())
        self.assertEqual(0, len(s))

    def test_statuslist(self):
        model = core.Model(tuple()) # getting away with an empty model for this step
        cache = hortator.PersistentTaskList(self.cache_file.name, model, force_create=True)
        statuslist = cache.statuslist
        self.assertEqual(set((y,x) for x,y in hortator._TASK_STATUS_LIST.items()), set(statuslist))

    def test_reopen(self):
        model = core.Model(tuple()) # getting away with an empty model for this step
        #create
        cache = hortator.PersistentTaskList(self.cache_file.name, model, force_create=True)
        #reopen
        cache_2 = hortator.PersistentTaskList(self.cache_file.name, model)
        # check that it was not created
        self.assertFalse(cache_2.created)
        # check that the statuslist is matching (as it also means that the inserts
        # to set up the DB were committed.
        self.assertEqual(cache.statuslist, cache_2.statuslist)

    def test_id_stored_entity(self):
        model = core.Model(tuple()) # getting away with an empty model for this step
        cache = hortator.PersistentTaskList(self.cache_file.name, rnaseq, force_create=True)
        db_id = cache.id_stored_entity(rnaseq.FASTQPossiblyGzipCompressed, 'hohoho.fastq')
        self.assertTrue(db_id.new)
        self.assertEqual(1, db_id.id)
        # same will give the same id.
        db_id_same = cache.id_stored_entity(rnaseq.FASTQPossiblyGzipCompressed, 'hohoho.fastq')
        self.assertFalse(db_id_same.new)
        self.assertEqual(db_id.id, db_id_same.id)
        # different will give a different id.
        db_id_differ = cache.id_stored_entity(rnaseq.FASTQPossiblyGzipCompressed, 'hahaha.fastq')
        self.assertTrue(db_id_differ.new)
        self.assertNotEqual(db_id.id, db_id_differ.id)

    def test_id_stored_sequence_len1(self):
        model = core.Model(tuple()) # getting away with an empty model for this step
        cache = hortator.PersistentTaskList(self.cache_file.name, rnaseq, force_create=True)
        db_id = cache.id_stored_sequence(core.FileSequence, ((core.File, 'hohoho.fastq'),))
        self.assertTrue(db_id.new)
        self.assertEqual(1, db_id.id)
        # same will give the same id.
        db_id_same = cache.id_stored_sequence(core.FileSequence, ((core.File, 'hohoho.fastq'),))
        self.assertFalse(db_id_same.new)
        self.assertEqual(db_id.id, db_id_same.id)
        # different will give a different id.
        db_id_differ = cache.id_stored_sequence(core.FileSequence, ((core.File, 'hahaha.fastq'),))
        self.assertTrue(db_id_differ.new)
        self.assertNotEqual(db_id.id, db_id_differ.id)

    def test_id_stored_sequence_len2(self):
        model = core.Model(tuple()) # getting away with an empty model for this step
        cache = hortator.PersistentTaskList(self.cache_file.name, rnaseq, force_create=True)
        db_id = cache.id_stored_sequence(core.FileSequence, ((core.File, 'hohoho.fastq'),
                                                             (core.File, 'hahaha.fastq')
                                                         ))
        self.assertTrue(db_id.new)
        self.assertEqual(1, db_id.id)
        # same will give the same id.
        db_id_same = cache.id_stored_sequence(core.FileSequence, ((core.File, 'hohoho.fastq'),
                                                                  (core.File, 'hahaha.fastq')
                                                            ))
        self.assertFalse(db_id_same.new)
        self.assertEqual(db_id.id, db_id_same.id)
        # different will give a different id.
        db_id_differ = cache.id_stored_sequence(core.FileSequence, ((core.File, 'hahaha.fastq'),
                                                                    (core.File, 'hohoho.fastq')
                                                                ))
        self.assertTrue(db_id_differ.new)
        self.assertNotEqual(db_id.id, db_id_differ.id)


    def test_id_step_activity(self):
        cache = hortator.PersistentTaskList(self.cache_file.name, rnaseq, force_create=True)
        # new activity
        db_id = cache.id_step_activity(rnaseq.ACTIVITY.ALIGN)
        self.assertTrue(db_id.new)
        self.assertEqual(1, db_id.id)
        # other activity - different ID
        db_id_other = cache.id_step_activity(rnaseq.ACTIVITY.QUANTIFY)
        self.assertTrue(db_id_other.new)
        self.assertNotEqual(db_id.id, db_id_other.id)
        db_id_another = cache.id_step_activity(rnaseq.ACTIVITY.DIFFEXP)
        self.assertTrue(db_id_another.new)
        self.assertNotEqual(db_id.id, db_id_another.id)
        # already know activity - do not create a new one
        db_id_same = cache.id_step_activity(rnaseq.ACTIVITY.ALIGN)
        self.assertFalse(db_id_same.new)
        self.assertEqual(db_id.id, db_id_same.id)

    def test_id_step_type_oneactivity(self):
        # step_type qualifies the activites carried out by a step
        cache = hortator.PersistentTaskList(self.cache_file.name, rnaseq, force_create=True)
        # new step_type that does only alignments
        db_id = cache.id_step_type((rnaseq.ACTIVITY.ALIGN, ))
        self.assertTrue(db_id.new)
        self.assertEqual(1, db_id.id)
        # step_type that does only alignments
        db_id = cache.id_step_type((rnaseq.ACTIVITY.ALIGN, ))
        self.assertFalse(db_id.new)
        self.assertEqual(1, db_id.id)

    def test_id_step_type_severalactivities(self):
        # step_type qualifies the activites carried out by a step
        cache = hortator.PersistentTaskList(self.cache_file.name, rnaseq, force_create=True)
        # new step_type that does pretty much anything
        db_id_all = cache.id_step_type((rnaseq.ACTIVITY.ALIGN,
                                        rnaseq.ACTIVITY.QUANTIFY,
                                        rnaseq.ACTIVITY.DIFFEXP))
        self.assertTrue(db_id_all.new)
        self.assertEqual(1, db_id_all.id)
        # same step_type - indedependent from the order in which activities are specified
        # (this is by design)
        db_id_all_same = cache.id_step_type((rnaseq.ACTIVITY.ALIGN,
                                             rnaseq.ACTIVITY.DIFFEXP,
                                             rnaseq.ACTIVITY.QUANTIFY))
        self.assertFalse(db_id_all_same.new)
        self.assertEqual(db_id_all.id, db_id_all_same.id)

    def test_id_step_variant(self):
        model = core.Model(tuple()) # getting away with an empty model for this step
        cache = hortator.PersistentTaskList(self.cache_file.name, model, force_create=True)
        class Activities(core.Enum):
            COMPRESS = 'Compress'
        class GZip(core.StepAbstract):
            _name = 'gzip'
            _default_execpath = 'gzip'
            activities = (core.DEFAULT_ACTIVITY.MISC, )
            Assets = None
            def __init__(self, executable=None):
                if executable is None:
                    self._execpath = _default_execpath
                else:
                    self._execpath = executable
            @property
            def version(self):
                res = subprocess.check_output([self._execpath, '--version'])
                version = res[:res.find(os.linesep)]
                return version
            def run(self):
                pass
        gzip = GZip('gzip')
        db_id = cache.id_step_variant(gzip,
                                      (Activities.COMPRESS,))
        self.assertTrue(isinstance(db_id.id, int))
        #FIXME: dummy executable (no sources or targets)

    def test_id_step_variant_noexecutable(self):
        model = core.Model(tuple()) # getting away with an empty model for this step
        cache = hortator.PersistentTaskList(self.cache_file.name, model, force_create=True)
        class Activities(core.Enum):
            NOODLE = 'Noodle'
        class Foo(core.StepAbstract):
            _name = 'foo'
            _default_execpath = None
            activities = (core.DEFAULT_ACTIVITY.MISC, )
            Assets = None
            def __init__(self, executable):
                if executable is not None:
                    raise ValueError('No executable needed. Should be None.')
                self._execpath = None
            @property
            def version(self):
                return '0.1.0'
            def run(self):
                pass
        foo = Foo(None)
        db_id = cache.id_step_variant(foo,
                                      (Activities.NOODLE,))
        self.assertTrue(db_id.new)
        self.assertTrue(isinstance(db_id.id, int))
        db_id_same = cache.id_step_variant(foo,
                                           (Activities.NOODLE,))
        self.assertFalse(db_id_same.new)
        self.assertEqual(db_id.id, db_id_same.id)
        
    def test_id_stepconcrete(self):
        model = core.Model(tuple()) # getting away with an empty model for this step
        cache = hortator.PersistentTaskList(self.cache_file.name, model, force_create=True)
        class Activities(core.Enum):
            DATETIME = 'Give date/time'
        PythonTime = self.PythonTime

        python = PythonTime('python')
        stepvariant_db_id = cache.id_step_variant(python,
                                                  (Activities.DATETIME,))
        # empty sources is a special case
        sources = core.AssetSet() # source
        targets = core.AssetSet() # targets
        parameters = tuple()
        db_id = cache.id_stepconcrete(stepvariant_db_id.id,
                                      sources, targets, parameters)
        db_id_same = cache.id_stepconcrete(stepvariant_db_id.id,
                                           sources, targets, parameters)
        self.assertEqual(db_id.id, db_id_same.id)
        db_id_notthesame = cache.id_stepconcrete(stepvariant_db_id.id,
                                                 sources, targets, parameters, tag = 2)
        self.assertNotEqual(db_id.id, db_id_notthesame.id)

        # 1-element sources
        sources = railroadtracks.model.aligners.AssetsIndexer.Source(rnaseq.SavedFASTA('foo.fasta'))

        db_id_nothesame = cache.id_stepconcrete(stepvariant_db_id.id,
                                                sources, targets, parameters)
        self.assertNotEqual(db_id.id, db_id_nothesame.id)
        db_id_sameagain = cache.id_stepconcrete(stepvariant_db_id.id,
                                                sources, targets, parameters)
        
        self.assertEqual(db_id_sameagain.id, db_id_nothesame.id)
        db_id_nothesameagain = cache.id_stepconcrete(stepvariant_db_id.id,
                                                     sources, targets, ("%Y",))
        self.assertNotEqual(db_id.id, db_id_nothesameagain.id)
        self.assertNotEqual(db_id_sameagain.id, db_id_nothesameagain.id)


        # 2-elements sources
        SrcCls = core.assetfactory('Source', 
                                   [core.AssetAttr('reference', rnaseq.SavedFASTA, ''),
                                    core.AssetAttr('otherreference', rnaseq.SavedFASTA, '')])
        sources = SrcCls(rnaseq.SavedFASTA('foo.fasta'),
                         rnaseq.SavedFASTA('bar.fasta'))
        db_id_notthesame = cache.id_stepconcrete(stepvariant_db_id.id,
                                                 sources, targets, parameters)
        self.assertNotEqual(db_id.id, db_id_notthesame.id)
        db_id_sameagain = cache.id_stepconcrete(stepvariant_db_id.id,
                                                sources, targets, parameters)
        self.assertEqual(db_id_sameagain.id, db_id_notthesame.id)


        # 1-element source / 1-element target
        sources = railroadtracks.model.aligners.AssetsIndexer.Source(rnaseq.SavedFASTA('foo.fasta'))
        targets = railroadtracks.model.aligners.AssetsIndexer.Target(rnaseq.FilePattern('foo_idx'))

        foo_sh = rnaseq.Anyscript()
        stepvariant_db_id = cache.id_step_variant(foo_sh,
                                                  (Activities.DATETIME,))

        db_id_nothesame = cache.id_stepconcrete(stepvariant_db_id.id,
                                                sources, targets, parameters)
        self.assertNotEqual(db_id.id, db_id_nothesame.id)
        db_id_sameagain = cache.id_stepconcrete(stepvariant_db_id.id,
                                                sources, targets, parameters)
        self.assertEqual(db_id_sameagain.id, db_id_nothesame.id)

        targets_bar = railroadtracks.model.aligners.AssetsIndexer.Target(rnaseq.FilePattern('bar_idx'))
        self.assertRaises(ValueError, cache.id_stepconcrete, stepvariant_db_id.id,
                          sources, targets_bar, parameters)

    def test_get_assets(self):
        model = core.Model(tuple()) # getting away with an empty model for this step
        cache = hortator.PersistentTaskList(self.cache_file.name, model, force_create=True)
        class Activities(core.Enum):
            DATETIME = 'Give date/time'
        PythonTime = self.PythonTime
        python = PythonTime('python')
        stepvariant_db_id = cache.id_step_variant(python,
                                                  (Activities.DATETIME,))
        # 2-elements sources
        SrcCls = core.assetfactory('Source', 
                                   [core.AssetAttr('reference', rnaseq.SavedFASTA, ''),
                                    core.AssetAttr('otherreference', rnaseq.SavedFASTA, ''),
                                    core.AssetAttr('listoffiles', rnaseq.SavedCSVSequence, '')])
        sources = SrcCls(rnaseq.SavedFASTA('foo.fasta'),
                         rnaseq.SavedFASTA('bar.fasta'),
                         rnaseq.SavedCSVSequence((rnaseq.SavedCSV('baz.csv'),rnaseq.SavedCSV('baz2.csv'))))
        targets = core.AssetSet() # targets
        parameters = tuple()
        db_id = cache.id_stepconcrete(stepvariant_db_id.id,
                                      sources, targets, parameters)
        for storedthing in cache.get_srcassets(db_id.id):
            thing = storedthing.resurrect(rnaseq)

    def test_nconcrete_steps(self):
        model = core.Model(tuple()) # getting away with an empty model for this step
        cache = hortator.PersistentTaskList(self.cache_file.name, model, force_create=True)
        class Activities(core.Enum):
            DATETIME = 'Give date/time'
        PythonTime = self.PythonTime
        python = PythonTime('python')
        stepvariant_db_id = cache.id_step_variant(python,
                                                  (Activities.DATETIME,))
        sources = core.AssetSet() # source
        targets = core.AssetSet() # targets
        parameters = tuple()
        self.assertEqual(0, cache.nconcrete_steps)
        taskstatuscount = cache.nconcrete_steps_status
        self.assertEqual(0, len(taskstatuscount))
        db_id = cache.id_stepconcrete(stepvariant_db_id.id,
                                      sources, targets, parameters)
        self.assertEqual(1, cache.nconcrete_steps)
        taskstatuscount = cache.nconcrete_steps_status
        self.assertEqual(taskstatuscount[0].label, hortator._TASK_TODO)
        self.assertEqual(taskstatuscount[0].count, 1)


    def test_gettargetsofactivity(self):
        model = core.Model(tuple()) # getting away with an empty model for this step
        cache = hortator.PersistentTaskList(self.cache_file.name, model, force_create=True)
        split = Split(None)
        parameters = tuple()
        input_file = tempfile.NamedTemporaryFile()
        input_file.write(b'123')
        input_file.flush()
        head_file = tempfile.NamedTemporaryFile()
        tail_file = tempfile.NamedTemporaryFile()
        sources = split.Assets.Source(core.File(input_file.name))
        targets = split.Assets.Target(core.File(head_file.name),
                                      core.File(tail_file.name))
        self.assertRaises(ValueError, cache.get_targetsofactivity, ActivitiesSplit.HEADTAIL)
        # ensure that step variant is tracked
        stepvariant_db_id = cache.id_step_variant(split,
                                                  split.activities)
        self.assertRaises(ValueError, cache.get_targetsofactivity, ActivitiesSplit.FOO)
        # get the targets of the activity (obviously there are not any yet)
        res = cache.get_targetsofactivity(ActivitiesSplit.HEADTAIL)
        self.assertEquals(0, len(res))
        # create a new task
        task_id = cache.id_stepconcrete(stepvariant_db_id.id,
                                        sources, targets, parameters)
        # there is only one such 
        res = cache.get_targetsofactivity(ActivitiesSplit.HEADTAIL)
        self.assertEquals(2, len(res))
        #
        head_file = tempfile.NamedTemporaryFile()
        tail_file = tempfile.NamedTemporaryFile()
        targets_other = split.Assets.Target(core.File(head_file.name),
                                      core.File(tail_file.name))
        task_id_other = cache.id_stepconcrete(stepvariant_db_id.id,
                                              sources, targets_other, parameters=(1,))
        res = cache.get_targetsofactivity(ActivitiesSplit.HEADTAIL)
        self.assertEquals(4, len(res))
        # --
        cat = Cat(None)
        concat = tempfile.NamedTemporaryFile()
        sources = cat.Assets.Source(core.FileSequence((targets.head, targets_other.head)))
        targets = cat.Assets.Target(core.File(concat.name))
        # ensure that step variant is tracked
        stepvariant_db_id = cache.id_step_variant(cat,
                                                  cat.activities)
        task_id_cat = cache.id_stepconcrete(stepvariant_db_id.id,
                                            sources, targets, parameters)
        res = cache.get_targetsofactivity(ActivitiesSplit.MERGE)
        self.assertEquals(1, len(res))
                                        
    def test_get_parenttask_of_storedentity(self):
        model = core.Model(tuple()) # getting away with an empty model for this step
        cache = hortator.PersistentTaskList(self.cache_file.name, model, force_create=True)
        split = Split(None)
        parameters = tuple()
        input_file = tempfile.NamedTemporaryFile()
        input_file.write(b'123')
        input_file.flush()
        head_file = tempfile.NamedTemporaryFile()
        tail_file = tempfile.NamedTemporaryFile()
        sources = split.Assets.Source(core.File(input_file.name))
        targets = split.Assets.Target(core.File(head_file.name),
                                      core.File(tail_file.name))
        # ensure that step variant is tracked
        stepvariant_db_id = cache.id_step_variant(split,
                                                  split.activities)
        # create a new task
        task_id = cache.id_stepconcrete(stepvariant_db_id.id,
                                        sources, targets, parameters)
        for sa in cache.get_srcassets(task_id):
            # all are root nodes
            res = cache.get_parenttask_of_storedentity(sa)
            self.assertTrue(res is None)
        for sa in cache.get_targetassets(task_id):
            res = cache.get_parenttask_of_storedentity(sa)
            self.assertEquals(task_id.id, res.id)
        #
        head_file_other = tempfile.NamedTemporaryFile()
        tail_file_other = tempfile.NamedTemporaryFile()
        targets_other = split.Assets.Target(core.File(head_file_other.name),
                                            core.File(tail_file_other.name))
        task_id_other = cache.id_stepconcrete(stepvariant_db_id.id,
                                              split.Assets.Source(targets.tail), targets_other, parameters=(1,))
        for sa in cache.get_srcassets(task_id_other):
            res = cache.get_parenttask_of_storedentity(sa)
            self.assertEquals(task_id.id, res.id)
        for sa in cache.get_targetassets(task_id_other):
            res = cache.get_parenttask_of_storedentity(sa)
            self.assertEquals(task_id_other.id, res.id)

        # --
        cat = Cat(None)
        sources_cat = cat.Assets.Source(core.FileSequence((targets.head, targets_other.head)))
        concat = tempfile.NamedTemporaryFile()
        targets_cat = cat.Assets.Target(core.File(concat.name))
        # ensure that step variant is tracked
        stepvariant_db_id_cat = cache.id_step_variant(cat,
                                                      cat.activities)
        task_id_cat = cache.id_stepconcrete(stepvariant_db_id_cat.id,
                                            sources_cat, targets_cat, parameters)
        for sa in cache.get_srcassets(task_id_cat):
            if hasattr(sa, 'iter_storedentities'):
                for sa_sub, t in zip(sa.iter_storedentities(), (task_id, task_id_other)):
                    res = cache.get_parenttask_of_storedentity(sa_sub)
                    self.assertEquals(t.id, res.id)                    
            else:
                res = cache.get_parenttask_of_storedentity(sa)
                self.assertEquals(task_id_cat.id, res.id)

        head_file_other = tempfile.NamedTemporaryFile()
        tail_file_other = tempfile.NamedTemporaryFile()
        targets_other = split.Assets.Target(core.File(head_file_other.name),
                                            core.File(tail_file_other.name))
        task_id = cache.id_stepconcrete(stepvariant_db_id.id,
                                        split.Assets.Source(targets_cat.result), targets_other, tuple())
        for sa in cache.get_srcassets(task_id):
            res = cache.get_parenttask_of_storedentity(sa)
            self.assertEquals(task_id_cat.id, res.id)
        


    def test_gettargetsoftype(self):
        model = core.Model(tuple()) # getting away with an empty model for this step
        cache = hortator.PersistentTaskList(self.cache_file.name, model, force_create=True)
        class SplitCSV(Split):
            _name = 'split'
            activities = (core.DEFAULT_ACTIVITY.MISC, )
            version = '0.1'
            _default_execpath = 'None'
            class Assets(core.AssetsStep):
                Source = core.assetfactory('Source', [core.AssetAttr('file', railroadtracks.model.files.SavedCSV, '')])
                Target = core.assetfactory('Target', [core.AssetAttr('head', railroadtracks.model.files.SavedCSV, ''),
                                                      core.AssetAttr('tail', railroadtracks.model.files.SavedCSV, '')])
            def run(self, assets, parameters=()):
                with open(assets.source.file, 'rb') as fh_in:
                    csv_r = csv.reader(fh_in)
                    with open(assets.target.head, 'wb') as fh_out:
                        csv_w = csv.writer(fh_out)
                        head = next(csv_r)
                        csv_w.writerow(head)
                    with open(assets.target.tail, 'wb') as fh_out:
                        csv_w = csv.writer(fh_out)
                        for row in csv_r:
                            csv_w.writerow(row)
                cmd = None
                returncode = 1
                return cmd, returncode

        split = Split(None)
        parameters = tuple()
        input_file = tempfile.NamedTemporaryFile()
        input_file.write(b'123')
        input_file.flush()
        head_file = tempfile.NamedTemporaryFile()
        tail_file = tempfile.NamedTemporaryFile()
        sources = split.Assets.Source(core.File(input_file.name))
        targets = split.Assets.Target(core.File(head_file.name),
                                      core.File(tail_file.name))
        stepvariant_db_id = cache.id_step_variant(split,
                                                  split.activities)
        res = cache.get_targetsoftype(core.File.__name__)
        self.assertEquals(0, len(res))
        task_id = cache.id_stepconcrete(stepvariant_db_id.id,
                                        sources, targets, parameters)
        res = cache.get_targetsoftype(core.File.__name__)
        self.assertEquals(2, len(res))

        splitcsv = SplitCSV(None)
        parameters = tuple()
        input_file = tempfile.NamedTemporaryFile(suffix=".csv")
        input_file.write(b'123')
        input_file.flush()
        head_file = tempfile.NamedTemporaryFile(suffix=".csv")
        tail_file = tempfile.NamedTemporaryFile(suffix=".csv")
        sources = splitcsv.Assets.Source(railroadtracks.model.files.SavedCSV(input_file.name))
        targets = splitcsv.Assets.Target(railroadtracks.model.files.SavedCSV(head_file.name),
                                         railroadtracks.model.files.SavedCSV(tail_file.name))
        stepvariant_db_id = cache.id_step_variant(splitcsv,
                                                  splitcsv.activities)
        res = cache.get_targetsoftype(railroadtracks.model.files.SavedCSV.__name__)
        self.assertEquals(0, len(res))
        task_id = cache.id_stepconcrete(stepvariant_db_id.id,
                                        sources, targets, parameters)
        res = cache.get_targetsoftype(railroadtracks.model.files.SavedCSV.__name__)
        self.assertEquals(2, len(res))

        head_file = tempfile.NamedTemporaryFile()
        tail_file = tempfile.NamedTemporaryFile()
        targets = split.Assets.Target(core.File(head_file.name),
                                      core.File(tail_file.name))
        task_id_other = cache.id_stepconcrete(stepvariant_db_id.id,
                                              sources, targets, parameters=(1,))
        res = cache.get_targetsoftype(core.File.__name__)
        self.assertEquals(4, len(res))

class ModelCRCHeadTailTestCase(unittest.TestCase):

    def _test(self, data, crc):
        fh = tempfile.NamedTemporaryFile()
        fh.write(data)
        fh.flush()
        out_fh = tempfile.NamedTemporaryFile(suffix='.csv')
        assets = crc.Assets(crc.Assets.Source(core.File(fh.name)),
                            crc.Assets.Target(rnaseq.SavedCSV(out_fh.name)))
        cmd, returncode = crc.run(assets)
        self.assertEqual(0, returncode)
        csv_r = csv.reader(out_fh)
        header = next(csv_r)
        filename, crc = next(csv_r)
        self.assertEqual(fh.name, filename)
        return crc

    def test_init(self):
        CRC = rnaseq.CRCHeadTail
        crc = CRC(None)
        data = 'foobarbaz'
        crc_a = self._test(data, crc)
        #
        crc_b = self._test(data, crc)
        #
        self.assertEqual(crc_a, crc_b)
        #
        data = '123456789'
        crc_c = self._test(data, crc)
        self.assertNotEqual(crc_a, crc_c)
        
        
class StepTestCase(unittest.TestCase):

    def test_Step(self):
        #'id label classname entityname'
        fh = tempfile.NamedTemporaryFile()
        fh.write('foobarbaz')
        fh.flush()
        out_fh = tempfile.NamedTemporaryFile(suffix='.csv')
        src = (hortator.StoredEntity(1,'file', 'File', fh.name), )
        targets = (hortator.StoredEntity(1,'crc', 'SavedCSV', out_fh.name), )
        #id status steptype_id executable clsname version
        sc = hortator.StepConcrete_DbEntry(1,
                                           hortator._TASK_TODO,
                                           1,
                                           1, None,
                                           'CRCHeadTail', None, ())
        parameters = tuple()
        step = hortator.Step(sc, src, targets, parameters,
                             rnaseq)
        self.assertFalse(step.iscomplete())
        returncode = step.run()
        self.assertEqual(0, returncode)
        sc = hortator.StepConcrete_DbEntry(1,
                                           hortator._TASK_DONE,
                                           1,
                                           1, None,
                                           'CRCHeadTail', None, ())
        model = core.Model(tuple()) # getting away with an empty model for this step
        step = hortator.Step(sc, src, targets, parameters, model)
        self.assertTrue(step.iscomplete())

class StepGraphTestCase(unittest.TestCase):
    
    def setUp(self):
        self.cache_file = tempfile.NamedTemporaryFile()
        model = core.Model(tuple()) # getting away with an empty model for this step
        self.cache = hortator.PersistentTaskList(self.cache_file.name, model, 
                                                 force_create=True)

    @unittest.skip('test not implemented')
    def test_StepGraph(self):
        raise NotImplementedError()
        #FIXME: where are the tests ?

    @unittest.skip('test not implemented')
    def test_addstep(self):
        raise NotImplementedError()
        #FIXME: where are the tests ?


class EasyTestCase(unittest.TestCase):

    def setUp(self):
        wd = tempfile.mkdtemp()
        self.wd = wd

    def tearDown(self):
        shutil.rmtree(self.wd)

    def test_FrozenNameSpace(self):
        kv = (('a', 1), ('b', 2))
        fz = easy.FrozenNamespace(kv)
        fz_keys = list(fz._fields)
        fz_keys.sort()
        for x,y in zip(kv, fz_keys):
            self.assertEqual(x[0],y)

    @unittest.skipIf(not environment.Executable.ispresent('bowtie2-build'),
                     'bowtie2-build is not in the PATH')            
    def test_command_line(self):
        packagedir = os.path.dirname(railroadtracks.__file__)
        PHAGEFASTA = os.path.join(packagedir, 'EF204940.FASTA')
        bowtie2index_exec = environment.Executable('bowtie2-build')
        project = easy.Project(rnaseq, self.wd)
        bowtie2 = rnaseq.Bowtie2Build()
        assets = bowtie2.Assets(bowtie2.Assets.Source(rnaseq.SavedFASTA(PHAGEFASTA)),
                                bowtie2.Assets.Target(rnaseq.SavedBowtie2Index(os.path.join(self.wd, 'foo'))))
        parameters = tuple()
        step_concrete_id = project.todo.add(bowtie2, assets, parameters = parameters)
        cmd = easy.command_line(project, step_concrete_id, bowtie2, assets, parameters = parameters)
        returncode = subprocess.check_call(cmd)
        self.assertEqual(0, returncode)

    @unittest.skipIf(not environment.Executable.ispresent('bowtie2-build'),
                     'bowtie2-build is not in the PATH')
    def test_EasyProject(self):
        packagedir = os.path.dirname(railroadtracks.__file__)
        PHAGEFASTA = os.path.join(packagedir, 'EF204940.FASTA')
        bowtie2index_exec = environment.Executable('bowtie2-build')
        project = easy.Project(rnaseq, self.wd)
        # Initial number of steps is zero
        self.assertEqual(0, project.todo._cache.nconcrete_steps)
        bowtie2index = rnaseq.Bowtie2Build(bowtie2index_exec.path)
        reference_fn = PHAGEFASTA
        task_index = project.add_task(bowtie2index, 
                                      bowtie2index.Assets(bowtie2index.Assets.Source(rnaseq.SavedFASTA(reference_fn))))
        # Number of steps is one
        self.assertEqual(1, project.todo._cache.nconcrete_steps)
        task_index_same = project.add_task(bowtie2index, 
                                           bowtie2index.Assets(bowtie2index.Assets.Source(rnaseq.SavedFASTA(reference_fn))))
        # Number of steps should remain one (same assets and paremeters)
        self.assertEqual(1, project.todo._cache.nconcrete_steps)

    @unittest.skipIf(not environment.Executable.ispresent('bowtie2-build'),
                     'bowtie2-build is not in the PATH')
    def test_Task(self):
        project = easy.Project(rnaseq, self.wd)
        #
        packagedir = os.path.dirname(railroadtracks.__file__)
        PHAGEFASTA = os.path.join(packagedir, 'EF204940.FASTA')
        bowtie2index_exec = environment.Executable('bowtie2-build')
        bowtie2align_exec = environment.Executable('bowtie2')

        bowtie2index = rnaseq.Bowtie2Build(bowtie2index_exec.path)
        reference_fn = PHAGEFASTA
        task_index = project.add_task(bowtie2index, 
                                      bowtie2index.Assets(bowtie2index.Assets.Source(rnaseq.SavedFASTA(reference_fn))))

        bowtie2align = rnaseq.Bowtie2(bowtie2align_exec.path)
        Assets = bowtie2align.Assets
        alignment_tasks = list()
        for read1_name, read2_name in (('f1_1.fq', 'f1_2.fq'),
                                       ('f2_1.fq', 'f2_2.fq')):
            assets = Assets(Assets.Source(task_index.call.assets.target.indexfilepattern, 
                                          rnaseq.FASTQPossiblyGzipCompressed(read1_name),
                                          rnaseq.FASTQPossiblyGzipCompressed(read2_name)),
                            None)
            alignment_tasks.append(project.add_task(bowtie2align, assets))

        # find all child
        res = task_index.all_child_tasks()
        s1 = set(x.task_id.id for x in alignment_tasks)
        s2 = set(x.task_id.id for x in res)
        self.assertEqual(0, len(s1.difference(s2)))
        self.assertEqual(0, len(s2.difference(s1)))

        # find the parents of an alignment
        task = next(iter(alignment_tasks))
        parent_taskset = task.parent_tasks()
        # only one parent (the indexing step)
        self.assertEqual(1, len(parent_taskset))
        parent_task = next(iter(parent_taskset))
        self.assertEqual(bowtie2index.activities, parent_task.call.step.activities)
        self.assertEqual(bowtie2index._name, parent_task.call.step._name)

    @unittest.skipIf(not environment.Executable.ispresent('bowtie2-build'),
                     'bowtie2-build is not in the PATH')
    def test_Task_withNone(self):
        project = easy.Project(rnaseq, self.wd)
        #
        packagedir = os.path.dirname(railroadtracks.__file__)
        PHAGEFASTA = os.path.join(packagedir, 'EF204940.FASTA')
        bowtie2index_exec = environment.Executable('bowtie2-build')
        bowtie2align_exec = environment.Executable('bowtie2')

        bowtie2index = rnaseq.Bowtie2Build(bowtie2index_exec.path)
        reference_fn = PHAGEFASTA
        task_index = project.add_task(bowtie2index, 
                                      bowtie2index.Assets(bowtie2index.Assets.Source(rnaseq.SavedFASTA(reference_fn))))

        bowtie2align = rnaseq.Bowtie2(bowtie2align_exec.path)
        Assets = bowtie2align.Assets
        alignment_tasks = list()
        for read1_name, read2_name in (('f1_1.fq', None),
                                       ('f2_1.fq', None)):
            assets = Assets(Assets.Source(task_index.call.assets.target.indexfilepattern, 
                                          rnaseq.FASTQPossiblyGzipCompressed(read1_name),
                                          read2_name),
                            None)
            alignment_tasks.append(project.add_task(bowtie2align, assets))

        # find all child
        res = task_index.all_child_tasks()
        s1 = set(x.task_id.id for x in alignment_tasks)
        s2 = set(x.task_id.id for x in res)
        self.assertEqual(0, len(s1.difference(s2)))
        self.assertEqual(0, len(s2.difference(s1)))

        # find the parents of an alignment
        task = next(iter(alignment_tasks))
        parent_taskset = task.parent_tasks()
        # only one parent (the indexing step)
        self.assertEqual(1, len(parent_taskset))
        parent_task = next(iter(parent_taskset))
        self.assertEqual(bowtie2index.activities, parent_task.call.step.activities)
        self.assertEqual(bowtie2index._name, parent_task.call.step._name)

# Test the writing of recipes
from railroadtracks import easy

class RecipeTestCase(unittest.TestCase):

    def setUp(self):
        # -- recipe-init-begin
        # -- initialization boiler plate code
        import tempfile
        from railroadtracks import hortator, rnaseq, easy
        from environment import Executable

        wd = tempfile.mkdtemp()
        project = easy.Project(rnaseq, wd=wd)

        # declare the 3rd-party command-line tools we will use
        env = easy.Environment(rnaseq)
        # -- recipe-init-end

        # -- recipe-data-begin
        # Phage genome shipped with the package for testing purposes
        import railroadtracks.model.simulate
        PHAGEFASTA = railroadtracks.model.simulate.PHAGEFASTA
        PHAGEGFF = railroadtracks.model.simulate.PHAGEGFF

        # create random data for 6 samples (just testing here)
        nsamples = 6
        samplereads = list()
        with open(PHAGEFASTA) as fasta_fh:
            reference = next(railroadtracks.model.simulate.readfasta_iter(fasta_fh))
        for sample_i in range(nsamples):
            read1_fh = tempfile.NamedTemporaryFile(prefix='read1', suffix='.fq')
            read2_fh = tempfile.NamedTemporaryFile(prefix='read2', suffix='.fq')
            read1_fh, read2_fh = railroadtracks.model.simulate.randomPEreads(read1_fh,
                                                                             read2_fh, 
                                                                             reference)
            samplereads.append((read1_fh, read2_fh))

        sampleinfo_fh = tempfile.NamedTemporaryFile(suffix='.csv')
        csv_w = csv.writer(sampleinfo_fh)
        csv_w.writerow(['sample_id', 'group'])
        for i in range(6):
            csv_w.writerow([str(i), ('A','B')[i%2]])
        sampleinfo_fh.flush()
        referenceannotation = rnaseq.SavedGFF(PHAGEGFF)
        # -- recipe-data-end

        self._wd = wd
        self.project = project
        self.reference_fn = PHAGEFASTA
        self.env = env
        self.nsamples = nsamples
        self.samplereads = samplereads
        self.sampleinfo_fh = sampleinfo_fh
        self.referenceannotation = referenceannotation
        self._PHAGEFASTA = PHAGEFASTA
        self._PHAGEGFF = PHAGEGFF

    def tearDown(self):
        samplereads = self.samplereads
        # -- recipe-teardown-begin
        for read1_fh, read2_fh in self.samplereads:
            read1_fh.close()
            read2_fh.close()
        # FIXME: delete the temporary directory
        shutil.rmtree(self.project.wd)
        # -- recipe-teardown-end

    def test_File(self):
        #FIXME: rather test it in the model ?
        reference = core.File(self.reference_fn)

    @unittest.skipIf(not (environment.Executable.ispresent('bowtie2-build') and \
                          environment.Executable.ispresent('htseq-count') and \
                          environment.Executable.ispresent('R') and \
                          environment.R('R').packageversion_or_none('edgeR') is not None),
                     'bowtie2 and/or htseq-count is not in the PATH')
    def test_RecipeSimpleIncremental(self):
        project = self.project
        env = self.env
        nsamples = self.nsamples
        samplereads = self.samplereads
        sampleinfo_fh = self.sampleinfo_fh
        reference_fn = self.reference_fn
        referenceannotation = self.referenceannotation
        PHAGEFASTA = self._PHAGEFASTA
        PHAGEGFF = self._PHAGEGFF

        # steps used
        bowtie2index = env.activities.INDEX.bowtie2build
        bowtie2align = env.activities.ALIGN.bowtie2
        htseqcount = env.activities.QUANTIFY.htseqcount
        merge = env.activities.UTILITY.columnmerger
        edger = env.activities.DIFFEXP.edger

        import easy

        # sequence of tasks to run
        torun = list()

        # index for alignment
        Assets = bowtie2index.Assets
        assets = Assets(Assets.Source(rnaseq.SavedFASTA(reference_fn)),
                        Assets.Target.createundefined())
        task_index = project.add_task(bowtie2index, 
                                      assets)
        # the step is not done
        self.assertEqual(hortator._TASK_TODO, task_index.info[1])
        torun.append(task_index)
        # run the tasks
        for task in torun:
            # run only if not done
            if task.info[1] != hortator._TASK_DONE:
                task.execute()
                task.status = hortator._TASK_DONE

        self.assertEqual(1, project.todo._cache.nconcrete_steps)
        # now that the tasks have run let's open the same project
        project_same = easy.Project(project.model, wd=project.wd)

        # index for alignment
        Assets = bowtie2index.Assets
        assets = Assets(Assets.Source(rnaseq.SavedFASTA(reference_fn)),
                        Assets.Target.createundefined())
        task_index_same = project_same.add_task(bowtie2index, 
                                                assets)

        self.assertNotEqual(task_index, task_index_same)
        self.assertNotEqual(task_index.call.assets, task_index_same.call.assets)
        self.assertListEqual(list(task_index.call.assets.source.reference), 
                             list(task_index_same.call.assets.source.reference))
        self.assertListEqual(list(task_index.call.assets.target.indexfilepattern), 
                             list(task_index_same.call.assets.target.indexfilepattern))
        self.assertEqual(hortator._TASK_DONE, task_index_same.info[1])
        self.assertEqual(1, project.todo._cache.nconcrete_steps)


    def _recipesimpleincremental(self, runtasks):
        project = self.project
        env = self.env
        nsamples = self.nsamples
        samplereads = self.samplereads
        sampleinfo_fh = self.sampleinfo_fh
        reference_fn = self.reference_fn
        referenceannotation = self.referenceannotation
        PHAGEFASTA = self._PHAGEFASTA
        PHAGEGFF = self._PHAGEGFF
        
        # steps used
        bowtie2index = env.activities.INDEX.bowtie2build
        bowtie2align = env.activities.ALIGN.bowtie2
        htseqcount = env.activities.QUANTIFY.htseqcount
        merge = env.activities.UTILITY.columnmerger
        edger = env.activities.DIFFEXP.edger

        for iteration in range(5):
            nextiteration = False
            # sequence of tasks to run
            torun = list()

            # index for alignment
            Assets = bowtie2index.Assets
            assets = Assets(Assets.Source(rnaseq.SavedFASTA(reference_fn)),
                            Assets.Target.createundefined())
            task_index = project.add_task(bowtie2index, assets)
            torun.append(task_index)
            if iteration < 1:
                nextiteration = True
                runtasks(torun)
                self.assertEqual(1, project.todo._cache.nconcrete_steps)
                continue
            # process all samples
            sample_counts = list()
            for sample_i, (read1_fh, read2_fh) in enumerate(samplereads):
                # align
                Assets = bowtie2align.Assets
                assets = Assets(Assets.Source(task_index.call.assets.target.indexfilepattern, 
                                              rnaseq.FASTQPossiblyGzipCompressed(read1_fh.name),
                                              rnaseq.FASTQPossiblyGzipCompressed(read2_fh.name)),
                                Assets.Target.createundefined())
                task_align = project.add_task(bowtie2align, assets)
                torun.append(task_align)
                if iteration < 2:
                    nextiteration = True
                    runtasks(torun)
                    self.assertEqual(1+(sample_i+1), project.todo._cache.nconcrete_steps)
                    continue

                # quantify
                # (non-default parameters to fit our demo GFF)
                params = rnaseq.HTSeqCount._noexons_parameters
                Assets = htseqcount.Assets
                assets = Assets(Assets.Source(task_align.call.assets.target.alignment,
                                              rnaseq.SavedGFF(referenceannotation)),
                                Assets.Target.createundefined())
                task_quantify = project.add_task(htseqcount,
                                                 assets,
                                                 parameters=params)
                torun.append(task_quantify)
                if iteration < 3:
                    nextiteration = True
                    runtasks(torun)
                    self.assertEqual(1+len(samplereads)+(sample_i+1), 
                                     project.todo._cache.nconcrete_steps)
                    continue

                # keep a pointer to the counts, as we will use it in the merge step
                sample_counts.append(task_quantify.call.assets)

            if nextiteration:
                continue
            # merge the sample data into a table (so differential expression can be computed)
            Assets = merge.Assets
            counts = tuple(x.target.counts for x in sample_counts)
            assets = Assets(Assets.Source(rnaseq.SavedCSVSequence(counts)),
                            merge.Assets.Target.createundefined())

            task_merge = project.add_task(merge,
                                          assets,
                                          parameters=("0", "1"))
            torun.append(task_merge)
            if iteration < 4:
                nextiteration = True
                runtasks(torun)
                self.assertEqual(1+2*len(samplereads)+1, 
                                 project.todo._cache.nconcrete_steps)
                continue

            # differential expression with edgeR
            Assets = edger.Assets
            assets = Assets(Assets.Source(task_merge.call.assets.target.counts,
                                          rnaseq.SavedCSV(sampleinfo_fh.name)),
                            Assets.Target.createundefined())
            task_de = project.add_task(edger,
                                       assets)
            if iteration < 5:
                nextiteration = True
                runtasks(torun)
                self.assertEqual(1+2*len(samplereads)+2, # 1 index + 2 FASTQ per sample + 1 merge + 1 differential expression
                                 project.todo._cache.nconcrete_steps)
                continue

    @unittest.skipIf(not (environment.Executable.ispresent('bowtie2-build') and \
                          environment.Executable.ispresent('htseq-count') and \
                          environment.Executable.ispresent('R') and \
                          environment.R('R').packageversion_or_none('edgeR') is not None),
                     'bowtie2 and/or htseq-count is not in the PATH')
    def test_RecipeSimpleIncrementalComplete(self):
        def runtasks(torun):
            # run the tasks
            for task in torun:
                # run only if not done
                if task.info[1] != hortator._TASK_DONE:
                    task.execute()
        self._recipesimpleincremental(runtasks)

    @unittest.skipIf(not (environment.Executable.ispresent('bowtie2-build') and \
                          environment.Executable.ispresent('htseq-count') and \
                          environment.Executable.ispresent('R') and \
                          environment.R('R').packageversion_or_none('edgeR') is not None),
                     'bowtie2, htseq-count, R (with package "edgeR") must be in the PATH')
    def test_RecipeSimpleIncrementalCompleteNoRun(self):
        def runtasks(torun):
            # do nothing
            pass
        self._recipesimpleincremental(runtasks)


    @unittest.skipIf(not (environment.Executable.ispresent('bowtie2-build') and \
                          environment.Executable.ispresent('htseq-count') and \
                          environment.Executable.ispresent('R') and \
                          environment.R('R').packageversion_or_none('edgeR') is not None),
                     'bowtie2, htseq-count, R (with package "edgeR") must be in the PATH')
    def test_RecipeSimple(self):
        project = self.project
        env = self.env
        nsamples = self.nsamples
        samplereads = self.samplereads
        sampleinfo_fh = self.sampleinfo_fh
        reference_fn = self.reference_fn
        referenceannotation = self.referenceannotation
        PHAGEFASTA = self._PHAGEFASTA
        PHAGEGFF = self._PHAGEGFF
        
        # -- recipesimple-test-begin

        # steps used
        bowtie2index = env.activities.INDEX.bowtie2build
        bowtie2align = env.activities.ALIGN.bowtie2
        htseqcount = env.activities.QUANTIFY.htseqcount
        merge = env.activities.UTILITY.columnmerger
        edger = env.activities.DIFFEXP.edger

        import easy

        # sequence of tasks to run
        torun = list()
                            
        # index for alignment
        Assets = bowtie2index.Assets
        assets = Assets(Assets.Source(rnaseq.SavedFASTA(reference_fn)),
                        Assets.Target.createundefined())
        task_index = project.add_task(bowtie2index, assets)
        torun.append(task_index)

        # process all samples
        sample_counts = list()
        for read1_fh, read2_fh in samplereads:
            # align
            Assets = bowtie2align.Assets
            assets = Assets(Assets.Source(task_index.call.assets.target.indexfilepattern, 
                                          rnaseq.FASTQPossiblyGzipCompressed(read1_fh.name),
                                          rnaseq.FASTQPossiblyGzipCompressed(read2_fh.name)),
                            Assets.Target.createundefined())
            task_align = project.add_task(bowtie2align, assets)
            torun.append(task_align)

            # quantify
            # (non-default parameters to fit our demo GFF)
            params = rnaseq.HTSeqCount._noexons_parameters
            Assets = htseqcount.Assets
            assets = Assets(Assets.Source(task_align.call.assets.target.alignment,
                                          rnaseq.SavedGFF(referenceannotation)),
                            Assets.Target.createundefined())
            task_quantify = project.add_task(htseqcount,
                                             assets,
                                             parameters=params)
            torun.append(task_quantify)
            # keep a pointer to the counts,
            # as we will use them in the merge step
            sample_counts.append(task_quantify.call.assets)

        # merge the sample data into a table
        # (so differential expression can be computed)
        Assets = merge.Assets
        counts = tuple(x.target.counts for x in sample_counts)
        assets = Assets(Assets.Source(rnaseq.SavedCSVSequence(counts)),
                        merge.Assets.Target.createundefined())
        task_merge = project.add_task(merge,
                                      assets,
                                      parameters=("0","1"))
        torun.append(task_merge)

        # differential expression with edgeR
        Assets = edger.Assets
        assets = Assets(Assets.Source(task_merge.call.assets.target.counts,
                                      rnaseq.SavedCSV(sampleinfo_fh.name)),
                        Assets.Target.createundefined())
        task_de = project.add_task(edger,
                                   assets)

        # run the tasks
        for task in torun:
            # run only if not done
            if task.info[1] != hortator._TASK_DONE:
                task.execute()

        # get results
        final_storedentities = project.get_targetsofactivity(rnaseq.ACTIVITY.DIFFEXP)

        # get the step that created the results files
        final_steps = list()
        for stored_entity in final_storedentities:
            final_steps.append(project.todo._cache.get_parenttask_of_storedentity(stored_entity))
        
        # -- recipesimple-test-end
        
        self.assertEqual(1, len(final_storedentities))
        self.assertEqual(core.File.__name__, final_storedentities[0].clsname)
        self.assertEqual('railroadtracks.model.diffexp.EdgeR', final_steps[0].clsname)

        # FIXME: not yet implemented
        # now that we have all steps, we "only" have to run them
        #steps = todo.stepcrawler()
        #for s in steps:
        #    print('%s' % (s.unifiedname))
        #    s.run()

    @unittest.skipIf(not (environment.Executable.ispresent('bowtie2-build') and \
                          environment.Executable.ispresent('bowtie-build') and \
                          environment.Executable.ispresent('STAR') and \
                          environment.Executable.ispresent('R') and \
                          environment.R('R').packageversion_or_none('edgeR') is not None and \
                          environment.R('R').packageversion_or_none('DESeq') is not None and \
                          environment.R('R').packageversion_or_none('DESeq2') is not None and \
                          environment.R('R').packageversion_or_none('limma') is not None),
                     'bowtie2, bowtie, STAR, TopHat2, and R (with packages "edgeR", "DESeq", "DESeq2", "limma") must be in the PATH')
    def test_RecipeLoop(self):
        project = self.project
        env = self.env
        nsamples = self.nsamples
        samplereads = self.samplereads
        sampleinfo_fh = self.sampleinfo_fh
        reference_fn = self.reference_fn
        referenceannotation = self.referenceannotation
        PHAGEFASTA = self._PHAGEFASTA
        PHAGEGFF = self._PHAGEGFF

        # -- recipeloop-test-begin
        import easy

        torun = list()

        # bowtie
        bowtie1index = env.activities.INDEX.bowtiebuild
        bowtie1align = env.activities.ALIGN.bowtie
        Assets = bowtie1index.Assets
        fa_file = rnaseq.SavedFASTA(reference_fn)
        task_index_bowtie1 = project.add_task(bowtie1index, 
                                              Assets(Assets.Source(fa_file),
                                                     None))
        torun.append(task_index_bowtie1)

        # bowtie2
        bowtie2index = env.activities.INDEX.bowtie2build
        bowtie2align = env.activities.ALIGN.bowtie2
        Assets = bowtie2index.Assets
        fa_file = rnaseq.SavedFASTA(reference_fn)
        task_index_bowtie2 = project.add_task(bowtie2index,
                                              Assets(Assets.Source(fa_file),
                                                     None))
        torun.append(task_index_bowtie2)

        # STAR
        starindex = env.activities.INDEX.starindex
        staralign = env.activities.ALIGN.staralign
        Assets = starindex.Assets
        fa_file = rnaseq.SavedFASTA(reference_fn)
        task_index_star = project.add_task(starindex, 
                                           Assets(Assets.Source(fa_file),
                                                  None))
        torun.append(task_index_star)

        # TopHat2
        # (index from bowtie2 used)
        #tophat2 = env.activities.ALIGN.tophat2

        # HTSeqCount
        htseqcount = env.activities.QUANTIFY.featurecount

        # Merge columns (obtained from counting)
        merge = env.activities.UTILITY.columnmerger

        # EdgeR, DESeq, DESeq2, and LIMMA voom
        edger = env.activities.DIFFEXP.edger
        deseq = env.activities.DIFFEXP.deseq
        deseq2 = env.activities.DIFFEXP.deseq2
        voom = env.activities.DIFFEXP.limmavoom
        

        # Now explore the different alignment presets in bowtie2, and vanilla star
        from itertools import cycle
        from collections import namedtuple
        Options = namedtuple('Options', 'aligner assets_index parameters')
        # Try various presets for bowtie2
        bowtie2_parameters = (('--very-fast', ), ('--fast', ), 
                              ('--sensitive', ), ('--very-sensitive', ))
        options = [Options(*x) for x in zip(cycle((bowtie2align,)),
                                            cycle((task_index_bowtie2.call.assets.target,)),
                                            bowtie2_parameters)]

        # add bowtie
        options.append(Options(bowtie1align, task_index_bowtie1.call.assets.target, tuple()))
        # add STAR (vanilla, no specific options beside the size of index k-mers)
        options.append(Options(staralign, 
                               task_index_star.call.assets.target, 
                               ('--genomeChrBinNbits', '12')))
        # add TopHat2
        #options.append(Options(tophat2, task_index_bowtie2.call.assets.target, tuple()))

        # loop over the options
        for option in options:
            sample_counts = list()
            # loop over the samples
            for sample_i in range(nsamples):
                read1_fh, read2_fh = samplereads[sample_i]
                # align
                Assets = option.aligner.Assets
                assets = Assets(Assets.Source(option.assets_index.indexfilepattern,
                                              rnaseq.FASTQPossiblyGzipCompressed(read1_fh.name), 
                                              rnaseq.FASTQPossiblyGzipCompressed(read2_fh.name)),
                                Assets.Target.createundefined())
                task_align = project.add_task(option.aligner,
                                              assets,
                                              parameters=option.parameters)
                torun.append(task_align)

                # quantify
                # (non-default parameters to fit our demo GFF)
                Assets = htseqcount.Assets
                assets = Assets(Assets.Source(task_align.call.assets.target.alignment,
                                              rnaseq.SavedGFF(referenceannotation)),
                                Assets.Target.createundefined())
                task_quantify = project.add_task(htseqcount,
                                                 assets)
                torun.append(task_quantify)

                # keep a pointer to the counts, as we will use it in the merge step
                sample_counts.append(task_quantify.call.assets)

            # merge the sample data into a table (so differential expression can be computed)
            Assets = merge.Assets
            source = Assets.Source(rnaseq.SavedCSVSequence(tuple(x.target.counts\
                                                                 for x in sample_counts)))
            assets_merge = Assets(source,
                                  Assets.Target.createundefined())
            task_merge = project.add_task(merge,
                                          assets_merge,
                                          parameters=("0","1"))
            torun.append(task_merge)

            # differential expression with edgeR, deseq2, and voom
            # (deseq is too whimsical for tests)
            for diffexp in (edger, deseq, deseq2, voom):
                Assets = diffexp.Assets
                assets = Assets(Assets.Source(task_merge.call.assets.target.counts,
                                              core.File(sampleinfo_fh.name)),
                                Assets.Target.createundefined())
                task_de = project.add_task(diffexp,assets)
                torun.append(task_de)

        # run the tasks
        for task in torun:
            if task.info[1] != hortator._TASK_DONE:
                try:
                    task.execute()
                    status = easy.hortator._TASK_DONE
                except:
                    status = easy.hortator._TASK_FAILED
            project.todo._cache.step_concrete_state(task.task_id,
                                                    easy.hortator._TASK_STATUS_LIST[status])
        # -- recipeloop-test-end


def suite():
    suite = unittest.TestLoader.loadTestsFromTestCase(ModelTestCase)
    suite.addTest(unittest.TestLoader.loadTestsFromTestCase(ModelSimulateTestCase))
    suite.addTest(unittest.TestLoader.loadTestsFromTestCase(AssetsTestCase))
    suite.addTest(unittest.TestLoader.loadTestsFromTestCase(ModelUtilsTestCase))
    suite.addTest(unittest.TestLoader.loadTestsFromTestCase(UnifexTestCase))
    suite.addTest(unittest.TestLoader.loadTestsFromTestCase(ModelAlignTestCase))
    suite.addTest(unittest.TestLoader.loadTestsFromTestCase(ModelAlignQuantificationTestCase))
    suite.addTest(unittest.TestLoader.loadTestsFromTestCase(ModelHandleBAMTestCase))
    suite.addTest(unittest.TestLoader.loadTestsFromTestCase(ModelQuantificationTestCase))
    suite.addTest(unittest.TestLoader.loadTestsFromTestCase(ModelColumnMergerTestCase))
    suite.addTest(unittest.TestLoader.loadTestsFromTestCase(ModelDExpressionTestCase))
    suite.addTest(unittest.TestLoader.loadTestsFromTestCase(ModelCRCHeadTailTestCase))
    suite.addTest(unittest.TestLoader.loadTestsFromTestCase(PersistanceTestCase))
    suite.addTest(unittest.TestLoader.loadTestsFromTestCase(StepTestCase))
    suite.addTest(unittest.TestLoader.loadTestsFromTestCase(StepGraphTestCase))
    suite.addTest(unittest.TestLoader.loadTestsFromTestCase(EasyTestCase))
    suite.addTest(unittest.TestLoader.loadTestsFromTestCase(RecipeTestCase))
    return suite

if __name__ == '__main__':
    print("""
Note: The tests are generating random sequencing data, occasionally causing STAR, DESeq, DESeq, or egdeR
to fail.
""")
    unittest.main()
