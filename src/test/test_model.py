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

import tempfile, shutil, os, subprocess, csv, random, re, itertools, gzip
import unittest
from railroadtracks import environment, unifex, core, rnaseq

from railroadtracks.model.simulate import (PHAGEGFF,
                                           PHAGEGTF,
                                           PHAGEFASTA)
import railroadtracks.model.aligners
import railroadtracks.model.simulate
import railroadtracks.model.diffexp
import railroadtracks.model.quantify

try:
    import ngs_plumbing
    has_ngsp = True
except ImportError:
    has_ngsp = False


def _build_StepIndex(cls, executable, 
                     tempdir, parameters=()):
    """ 
    :param cls: class used to model the step
    :param executable: name (or path) of executable 
    :param tempdir: temporary directory 

    :return (assets, cmd, returncode):"""
    runner = cls(executable)
    AssetsIndexer = railroadtracks.model.aligners.AssetsIndexer
    referencegenome = rnaseq.FASTAFile(PHAGEFASTA)
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
    if cls_alignedreads is rnaseq.SAMFile:
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
                     read1_fh, tempdir,
                     parameters_index = (),
                     align_max_attempts = 1):
    assets_index, cmd, returncode = _build_StepIndex(cls_index, executable_index,
                                                     tempdir, 
                                                     parameters = parameters_index)
    for n_attempt in range(1, align_max_attempts+1):
        try:
            version, assets, cmd, returncode, fh = _build_StepAlign(assets_index,
                                                                    cls_align, executable_align,
                                                                    read1_fh)
            break
        except AssertionError as ae:
            if n_attempt == align_max_attempts:
                raise(ase)

    return (assets, cmd, returncode, fh)


def _build_StepQuantify(assets_align,
                        cls_count, executable_count,
                        parameters_count):
    #FIXME: should be sorted by ID !
    runner = cls_count(executable_count)
    AssetsQuantifier = quantify.AssetsQuantifier
    alignmentfiles = tuple(assets_align.target.alignment)
    #FIXME: len == 1 (tested elsewhere)
    f_align = rnaseq.BAMFile(os.path.join(os.path.dirname(assets_align.target.alignment.name),
                                          alignmentfiles[0]))
    referenceannotation = rnaseq.GFFFile(PHAGEGFF)
    source = AssetsQuantifier.Source(f_align,
                                     referenceannotation)
    fh = tempfile.NamedTemporaryFile(suffix='.csv')
    counts = rnaseq.CSVFile(fh.name)
    target = AssetsQuantifier.Target(counts)
    assets = AssetsQuantifier(source, target)
    cmd, returncode = runner.run(assets, parameters = parameters_count)
    return (assets, cmd, returncode, fh) # fh returned to be protected from premature destruction
    


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
        assets = A(A.Source(rnaseq.GFFFile(annotation_fn),
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
        self.assertEqual(0, returncode)
        # check that the target file is not empty
        self.assertGreater(os.stat(assets.target.parameters_and_pro.name + '.pro').st_size, 0)

    @unittest.skipIf(not environment.Executable.ispresent('flux-simulator'),
                     'flux-simulator is not in the PATH')
    def test_FluxsimulatorDerivedExpression(self):
        fse = rnaseq.FluxsimulatorExpression()
        assets = self._test_FluxsimulatorExpression(fse)
        cmd, returncode = fse.run(assets, parameters=('--NB-MOLECULES', '1000'))
        # now devire a new expression profile
        fse = rnaseq.FluxsimulatorDerivedExpression()
        A = fse.Assets
        paramsandpro_prefix_source = os.path.join(self._dir, 'modified')
        paramsandpro_prefix_target = os.path.join(self._dir, 'modified_again')
        assets = A(A.Source(rnaseq.FluxSimulatorParametersAndPro(paramsandpro_prefix_source)),
                   A.Target(rnaseq.FluxSimulatorParametersAndPro(paramsandpro_prefix_target)))
        cmd, returncode = fse.run(assets)
        
        self.assertEqual(0, returncode)
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
        self.assertEqual(0, returncode)
        # check that the target files are not empty
        self.assertGreater(os.stat(assets.target.parameters_and_pro.name + '.pro').st_size, 0)
        self.assertGreater(os.stat(read1_fn).st_size, 0)
        self.assertGreater(os.stat(read2_fn).st_size, 0)
        self.assertGreater(os.stat(lib_fn).st_size, 0)



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
        self._read1_fh, self._read2_fh = railroadtracks.model.simulate.randomPEreads(read1_fh, 
                                                                                     read2_fh, 
                                                                                     reference,
                                                                                     n = 500)
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
        self.assertEqual(fh_sam.name, fh_sam_again.name)

    @unittest.skipIf(not environment.Executable.ispresent('bowtie2-build'),
                     'bowtie2-build is not in the PATH')
    def test_ensurebam(self):
        assets_align = self._assets_align
        fn = assets_align.target.alignment.name
        self.assertTrue(fn.endswith('.bam'))
        fh_sam = railroadtracks.model.files.ensure_sam(fn)
        fh_bam = railroadtracks.model.files.ensure_bam(fh_sam.name)
        fh_bam_again = railroadtracks.model.files.ensure_bam(fh_bam.name)
        self.assertEqual(fh_bam.name, fh_bam_again.name)
        self.assertTrue(fh_bam.name.endswith('.bam'))
        # this is round trip for the SAM file
        with open(fn, mode='rb') as fh_bam_orig:
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
        target = Assets.Target(rnaseq.SAMFile(sam_fh.name))
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
                                 Assets.Target(rnaseq.SAMFile(sam_fh.name)))
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
        inputs = (tempfile.NamedTemporaryFile(mode='w+'), 
                  tempfile.NamedTemporaryFile(mode='w+'))
        inputs[0].write('abc')
        inputs[0].flush()
        inputs[1].write('def')
        inputs[1].flush()
        outputs = (tempfile.NamedTemporaryFile(),)
        sources = {'input': tuple(x.name for x in inputs, )}
        targets = {'output': tuple(x.name for x in outputs)}
        res = cs.run(sources, targets, '')
        out_res = outputs[0].readlines()
        self.assertEqual([b'abcdef'], out_res)


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
        self.assertTrue(isinstance(version, bytes))
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
        self.assertTrue(isinstance(version, bytes))
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
                                                   parameters=('--genomeChrBinNbits', '12',
                                                               '--genomeSAindexNbases', '8'))
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
                                                   os.path.join(self.tempdir, 'reference'),
                                                   parameters = ('-k', '8'))
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

    @unittest.skipIf(not environment.Executable.ispresent('salmon'),
                     'salmon is not in the PATH')
    def test_StepSalmonIndex(self):
        execname = 'salmon'
        SalmonIndex = railroadtracks.model.aligners.SalmonIndex
        runner = SalmonIndex(execname)
        #FIXME: test version
        version = runner.version
        self.assertTrue(isinstance(version, str))
        self.assertEqual(set(runner.activities), set((rnaseq.ACTIVITY.INDEX, )))
        assets, cmd, returncode = _build_StepIndex(SalmonIndex, execname, 
                                                   os.path.join(self.tempdir, 'reference'),
                                                   parameters = ('--threads', '1'))
        self.assertEqual(0, returncode)


class ModelAlignTestCase(unittest.TestCase):

    def setUp(self):
        NFRAGMENTS_NOMATCH = 10
        read1_fh = tempfile.NamedTemporaryFile(prefix='read1', suffix='.fq')
        read2_fh = tempfile.NamedTemporaryFile(prefix='read2', suffix='.fq')
        with open(PHAGEFASTA) as fasta_fh:
            reference = next(railroadtracks.model.simulate.readfasta_iter(fasta_fh))
        self._read1_fh, self._read2_fh = railroadtracks.model.simulate.randomPEreads(read1_fh, 
                                                                                     read2_fh, 
                                                                                     reference,
                                                                                     n=500)
        # reads from a random reference (should not match)
        randint = railroadtracks.model.simulate.random.randint
        rand_reference = railroadtracks.model.simulate.Entry('> Random DNA',
                                                             bytearray(''.join('ATCG'[randint(0, 3)] for x in range(500)), 'ascii'))
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
                                cls_align, executable_align,
                                parameters_index = ()):
        assets_index, cmd, index_res = _build_StepIndex(cls_index, executable_index,
                                                        os.path.join(self.tempdir, 'reference'),
                                                        parameters = parameters_index)
        runner = cls_align(executable_align)
        self.assertEqual(set(runner.activities), set((rnaseq.ACTIVITY.ALIGN, )))
        AssetsAligner = railroadtracks.model.aligners.AssetsAligner
        # no reads
        self.assertRaises(AssertionError,
                          AssetsAligner.Source,
                          assets_index.target.indexfilepattern)

    def _test_StepAlign_singlereads(self, 
                                    cls_index, executable_index,
                                    cls_align, executable_align,
                                    parameters_index = (),
                                    align_max_attempts = 1):
        runner = cls_align(executable_align)
        self.assertEqual(set(runner.activities), set((rnaseq.ACTIVITY.ALIGN, )))
        # single reads
        assets, cmd, returncode, fh = _build_UpToAlign(cls_index, executable_index,
                                                       cls_align, executable_align, self._read1_fh,
                                                       os.path.join(self.tempdir, 'reference'),
                                                       parameters_index = parameters_index,
                                                       align_max_attempts = align_max_attempts)
        self.assertTrue(isinstance(runner.version, bytes))
        self.assertEqual(0, returncode)
        # FIXME: check that the alignment file contains what it should
        self.assertTrue(os.path.exists(fh.name))
        # check that the target file is not empty
        self.assertGreater(os.stat(fh.name).st_size, 0)

    def _test_StepAlign_pairedreads(self, 
                                    cls_index, executable_index,
                                    cls_align, executable_align,
                                    parameters_index = (),
                                    align_max_attempts = 1):
        assets_index, cmd, returncode = _build_StepIndex(cls_index, executable_index,
                                                         os.path.join(self.tempdir, 'reference'),
                                                         parameters = parameters_index)
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
        if cls_alignedreads is rnaseq.SAMFile:
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
            for n_attempt in range(1, align_max_attempts+1):
                cmd, returncode = runner.run(assets)
                if returncode == 0:
                    break
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
    def test_StepGsnapAlign_singlereads(self):
        self._test_StepAlign_singlereads(rnaseq.GsnapIndex, 'gmap_build',
                                         rnaseq.GsnapAlign, 'gsnap',
                                         parameters_index = ('-k', '8'))
    @unittest.skipIf(not environment.Executable.ispresent('gsnap'),
                     'gsnap is not in the PATH')
    @unittest.skipIf(not environment.Executable.ispresent('gmap_build'),
                     'gmap_build is not in the PATH')
    def test_StepGsnapAlign_pairedreads(self):
        self._test_StepAlign_pairedreads(rnaseq.GsnapIndex, 'gmap_build',
                                         rnaseq.GsnapAlign, 'gsnap',
                                         parameters_index = ('-k', '8'))

    @unittest.skipIf(not environment.Executable.ispresent('STAR'),
                     'STAR is not in the PATH')
    def test_StepStarAlign_singlereads(self):
        # STAR is whimsical. Try running the tests up to 5 times.
        self._test_StepAlign_singlereads(rnaseq.StarIndex, 'STAR',
                                         rnaseq.StarAlign, 'STAR',
                                         parameters_index = ('--genomeChrBinNbits', '12',
                                                             '--genomeSAindexNbases', '8'),
                                         align_max_attempts = 3)

    @unittest.skipIf(not environment.Executable.ispresent('STAR'),
                     'STAR is not in the PATH')
    def test_StepStarAlign_pairedreads(self):
        # STAR is whimsical. Try running the tests up to 5 times.
        self._test_StepAlign_pairedreads(rnaseq.StarIndex, 'STAR',
                                         rnaseq.StarAlign, 'STAR',
                                         parameters_index = ('--genomeChrBinNbits', '12',
                                                             '--genomeSAindexNbases', '8'),
                                         align_max_attempts = 3)


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
                                                             bytearray(''.join('ATCG'[randint(0, 3)] for x in range(500)), 'ascii'))
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
        SamtoolsSorterByID = rnaseq.SamtoolsSorterByID
        AssetsSorter = SamtoolsSorterByID.Assets
        executable_sort = 'samtools'
        samsort_byID = SamtoolsSorterByID(executable_sort)
        self.assertEqual(set((rnaseq.ACTIVITY.SORT,)), set(samsort_byID.activities))
        # mostly to check that version is working
        self.assertTrue(isinstance(samsort_byID.version, bytes))
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
        self.assertTrue(isinstance(extractunaligned.version, bytes))
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
        self.assertTrue(isinstance(readfilter.version, bytes))
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
    def setUp(self):
        self.tempdir = tempfile.mkdtemp()
        read1_fh = tempfile.NamedTemporaryFile(prefix='read1', suffix='.fq', dir=self.tempdir)
        read2_fh = tempfile.NamedTemporaryFile(prefix='read2', suffix='.fq', dir=self.tempdir)
        with open(PHAGEFASTA) as fasta_fh:
            reference = next(railroadtracks.model.simulate.readfasta_iter(fasta_fh))
            # set the random seed
            random.seed(a=123)
        self._read1_fh, self._read2_fh = railroadtracks.model.simulate.randomPEreads(read1_fh, read2_fh, reference)

    def tearDown(self):
        self._read1_fh.close()
        self._read2_fh.close()
        shutil.rmtree(self.tempdir)

    @unittest.skipIf(not environment.Executable.ispresent('sailfish'),
                     'sailfish is not in the PATH')        
    def test_pairedend_Sailfish(self):
        execname = 'sailfish'
        assets, cmd, returncode = _build_StepIndex(railroadtracks.model.aligners.SailfishIndex(execname),
                                                   execname, 
                                                   os.path.join(self.tempdir, 'reference'),
                                                   parameters = runner.PARAMETERS_DEFAULT)
        if returncode != 0:
            raise Exception("Index not built")
        self.assets = assets
        sf = railroadtracks.model.quantify.SailfishQuant(execname)
        A = sf.Assets
        count_csv = os.path.join(self.tempdir, 'counts.csv')
        output_dir = os.path.join(self.tempdir, 'sailfish_output_dir')
        assets = A(A.Source(self.assets.target.indexfilepattern,
                            rnaseq.FASTQPossiblyGzipCompressed(self._read1_fh.name),
                            rnaseq.FASTQPossiblyGzipCompressed(self._read2_fh.name)),
                   A.Target(rnaseq.CSVFile(count_csv),
                            rnaseq.FilePattern(output_dir)))
        cmd, returncode = sf.run(assets, parameters = ('--libtype', sf.LIBRARY_PE))

    @unittest.skipIf(not environment.Executable.ispresent('salmon'),
                     'salmon is not in the PATH')
    def test_pairedend_Salmon(self):
        execname = 'salmon'
        SalmonIndex = railroadtracks.model.aligners.SalmonIndex
        assets, cmd, returncode = _build_StepIndex(SalmonIndex, execname,
                                                   os.path.join(self.tempdir, 'reference'),
                                                   parameters = ('--threads', '1'))
        if returncode != 0:
            raise Exception("Index not built")
        self.assets = assets
        sm = railroadtracks.model.quantify.SalmonAlignQuant(execname)
        A = sm.Assets
        count_csv = os.path.join(self.tempdir, 'counts.csv')
        output_dir = os.path.join(self.tempdir, 'salmon_output_dir')
        assets = A(A.Source(self.assets.target.indexfilepattern,
                            rnaseq.FASTQPossiblyGzipCompressed(self._read1_fh.name),
                            rnaseq.FASTQPossiblyGzipCompressed(self._read2_fh.name)),
                   A.Target(rnaseq.CSVFile(count_csv),
                            rnaseq.FilePattern(output_dir)))
        cmd, returncode = sm.run(assets, parameters = ('--libType', sm.LIBRARY_PE))

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
        AssetsQuantifier = railroadtracks.model.quantify.AssetsQuantifier
        HTSeqCount = rnaseq.HTSeqCount
        referenceannotation = rnaseq.GFFFile(PHAGEGFF)
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
            
            counts = rnaseq.CSVFile(fh.name)
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
        AssetsQuantifier = railroadtracks.model.quantify.AssetsQuantifier
        FeatureCount = rnaseq.FeatureCount
        referenceannotation = rnaseq.GFFFile(PHAGEGFF)
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
        counts = rnaseq.CSVFile(fh.name)
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
        with tempfile.NamedTemporaryFile(mode='w+', suffix=".csv") as fh_col1, \
             tempfile.NamedTemporaryFile(mode='w+', suffix=".csv") as fh_col2, \
             tempfile.NamedTemporaryFile(mode='w+', suffix=".csv") as fh_merged:
            fh_col1_csv = csv.writer(fh_col1)
            fh_col2_csv = csv.writer(fh_col2)
            for x in ('1','2','3'):
                fh_col1_csv.writerow([x])
            fh_col1.seek(0)
            for x in ('4','5','6'):
                fh_col2_csv.writerow([x])
            fh_col2.seek(0)
            ColumnMerger = rnaseq.ColumnMerger
            Assets = ColumnMerger.Assets
            merger = ColumnMerger()
            filestomerge = tuple(rnaseq.CSVFile(x.name) for x in (fh_col1, fh_col2))
            source = Assets.Source(rnaseq.CSVFileSequence(filestomerge))
            assets = Assets(source,
                            Assets.Target(rnaseq.CSVFile(fh_merged.name)))
            cmd, returncode = merger.run(assets,
                                         parameters=('None', '0'))
            self.assertEqual(0, returncode)
            # just in case... rewind to 0
            fh_merged.seek(0)
            self.assertEqual(('1,4',
                              '2,5',
                              '3,6'), tuple(x.rstrip() for x in fh_merged))

    def test_merge_withIDs(self):
        with tempfile.NamedTemporaryFile(suffix=".csv", mode='w+') as fh_col1, \
             tempfile.NamedTemporaryFile(suffix=".csv", mode='w+') as fh_col2, \
             tempfile.NamedTemporaryFile(suffix=".csv", mode='w+') as fh_merged:
            fh_col1_csv = csv.writer(fh_col1)
            fh_col2_csv = csv.writer(fh_col2)
            for row in (('1','a'),('2','b'),('3','c')):
                fh_col1_csv.writerow(row)
            fh_col1.flush()
            fh_col1.seek(0)
            for row in (('1','d'),('2','e'),('3','f')):
                fh_col2_csv.writerow(row)
            fh_col2.flush()
            fh_col2.seek(0)
            ColumnMerger = rnaseq.ColumnMerger
            Assets = ColumnMerger.Assets
            merger = ColumnMerger()
            filestomerge = tuple(rnaseq.CSVFile(x.name) for x in (fh_col1, fh_col2))
            source = Assets.Source(rnaseq.CSVFileSequence(filestomerge))
            assets = Assets(source,
                            Assets.Target(rnaseq.CSVFile(fh_merged.name)))
            cmd, returncode = merger.run(assets, parameters=('0', '1'))
            self.assertEqual(0, returncode)
            # just in case... rewind to 0
            fh_merged.seek(0)
            self.assertEqual(('1,a,d',
                              '2,b,e',
                              '3,c,f'),
                             tuple(x.rstrip() for x in fh_merged))

class ModelDExpressionTestCase(unittest.TestCase):
    """
    """
    
    def setUp(self):
        self.tempdir = tempfile.mkdtemp()
        sampleinfo_fh = tempfile.NamedTemporaryFile(dir=self.tempdir, 
                                                    mode='w+',
                                                    suffix='.csv', 
                                                    delete=False)
        csv_w = csv.writer(sampleinfo_fh)
        csv_w.writerow(['sample_id', 'group'])
        for i in range(6):
            csv_w.writerow([str(i), ('A','B')[i%2]])
        sampleinfo_fh.flush()
        self._sampleinfo_fh = sampleinfo_fh
        with open(PHAGEFASTA) as fasta_fh:
            reference = next(railroadtracks.model.simulate.readfasta_iter(fasta_fh))
        random.seed(123)    
        with tempfile.NamedTemporaryFile(dir=self.tempdir,
                                         mode='w+',
                                         suffix='.csv', 
                                         delete=False) as fh_merged:
            csv_w = csv.writer(fh_merged)
            self.asset_merge = rnaseq.CSVFile(fh_merged.name)
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
        fh = tempfile.NamedTemporaryFile(mode='w+')
        source = AssetsDE.Source(self.asset_merge, 
                                 core.File(self._sampleinfo_fh.name))
        target = AssetsDE.Target(core.File(fh.name))
        assets = AssetsDE(source,
                          target)
        cmd, returncode = deseq.run(assets, parameters=('--dispersion-fittype=local', ))
        self.assertEqual(0, returncode)

    @unittest.skipIf(not (environment.Executable.ispresent('R') and \
                          environment.R('R').packageversion_or_none('DESeq2') is not None),
                     'R (with package "DESeq2") must be in the PATH')
    def test_StepDESeq2(self):
        deseq2 = rnaseq.DESeq2(self._r_exec)
        self.assertEqual(set((rnaseq.ACTIVITY.DIFFEXP,)), set(deseq2.activities))
        AssetsDE = railroadtracks.model.diffexp.AssetsDifferentialExpression
        fh = tempfile.NamedTemporaryFile(mode='w+')
        source = AssetsDE.Source(self.asset_merge,
                                 core.File(self._sampleinfo_fh.name))
        target = AssetsDE.Target(core.File(fh.name))
        assets = AssetsDE(source,
                          target)
        cmd, returncode = deseq2.run(assets)
        self.assertEqual(0, returncode)


    @unittest.skipIf(not (environment.Executable.ispresent('R') and \
                          environment.R('R').packageversion_or_none('DESeq2') is not None),
                     'R (with package "edgeR") must be in the PATH')
    def test_StepEdgeR(self):
        edger = rnaseq.EdgeR(self._r_exec)
        self.assertEqual(set((rnaseq.ACTIVITY.DIFFEXP,)), set(edger.activities))
        AssetsDE = railroadtracks.model.diffexp.AssetsDifferentialExpression
        fh = tempfile.NamedTemporaryFile(mode='w+')
        source = AssetsDE.Source(self.asset_merge,
                                 core.File(self._sampleinfo_fh.name))
        target = AssetsDE.Target(core.File(fh.name))
        assets = AssetsDE(source,
                          target)
        cmd, returncode = edger.run(assets)
        self.assertEqual(0, returncode)


    @unittest.skipIf(not (environment.Executable.ispresent('R') and \
                          environment.R('R').packageversion_or_none('limma') is not None),
                     'R (with package "limma") must be in the PATH')
    def test_StepLimmaVoom(self):
        voom = rnaseq.LimmaVoom(self._r_exec)
        self.assertEqual(set((rnaseq.ACTIVITY.DIFFEXP,)), set(voom.activities))
        AssetsDE = railroadtracks.model.diffexp.AssetsDifferentialExpression
        fh = tempfile.NamedTemporaryFile(mode='w+')
        source = AssetsDE.Source(self.asset_merge,
                                 core.File(self._sampleinfo_fh.name))
        target = AssetsDE.Target(core.File(fh.name))
        assets = AssetsDE(source,
                          target)
        cmd, returncode = voom.run(assets)
        self.assertEqual(0, returncode)


    @unittest.skipIf(not (environment.Executable.ispresent('R') and \
                          environment.R('R').packageversion_or_none('EBSeq') is not None),
                     'R (with package "EBSeq") must be in the PATH')
    def test_StepEBSeq(self):
        ebseq = railroadtracks.model.diffexp.EBSeq(self._r_exec)
        self.assertEqual(set((railroadtracks.model.diffexp.ACTIVITY.DIFFEXP,)), set(ebseq.activities))
        AssetsDE = railroadtracks.model.diffexp.AssetsDifferentialExpression
        fh = tempfile.NamedTemporaryFile(mode='w+')
        source = AssetsDE.Source(self.asset_merge,
                                 core.File(self._sampleinfo_fh.name))
        target = AssetsDE.Target(core.File(fh.name))
        assets = AssetsDE(source,
                          target)
        cmd, returncode = ebseq.run(assets)
        self.assertEqual(0, returncode)


class ModelCRCHeadTailTestCase(unittest.TestCase):

    def _test(self, data, crc):
        fh = tempfile.NamedTemporaryFile(mode='w+')
        fh.write(data)
        fh.flush()
        out_fh = tempfile.NamedTemporaryFile(mode='w+', suffix='.csv')
        assets = crc.Assets(crc.Assets.Source(core.File(fh.name)),
                            crc.Assets.Target(rnaseq.CSVFile(out_fh.name)))
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


class GzipFastqFilePairTestCase(unittest.TestCase):

    def setUp(self):
        self.tempdir = tempfile.mkdtemp()
        read1_fh = tempfile.NamedTemporaryFile(prefix='read1', suffix='.fq.gz', dir=self.tempdir, delete=False)
        read1_fh.close()
        self.read1_fn = read1_fh.name

        read2_fh = tempfile.NamedTemporaryFile(prefix='read2', suffix='.fq.gz', dir=self.tempdir, delete=False)
        read2_fh.close()
        self.read2_fn = read2_fh.name

    def tearDown(self):
        shutil.rmtree(self.tempdir)

    @unittest.skipIf(not has_ngsp,
                     'The Python package ngs-plumbing is missing.')
    def test_GzipFastqFilePair(self):
        NFRAGMENTS_MATCH = 300
        with open(PHAGEFASTA) as fasta_fh:
            reference = next(railroadtracks.model.simulate.readfasta_iter(fasta_fh))
            read1_io = gzip.GzipFile(self.read1_fn, mode='w')
            read2_io = gzip.GzipFile(self.read2_fn, mode='w')
            # reads from the target genome
            read1_fh, read2_fh=railroadtracks.model.simulate.randomPEreads(read1_io,
                                                                           read2_io,
                                                                           reference,
                                                                           n = NFRAGMENTS_MATCH)
        read1_fh.close()
        read2_fh.close()
        fqp = rnaseq.GzipFastqFilePair(self.read1_fn,
                                       self.read2_fn)
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
