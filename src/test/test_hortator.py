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

import sys, os, tempfile, subprocess
import unittest

from railroadtracks import core, rnaseq, hortator
import railroadtracks.model.aligners

if sys.version_info[0] < 3:
    linesep = os.linesep
else:
    linesep = bytes(os.linesep, encoding='ascii')

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


class StepTestCase(unittest.TestCase):

    def test_Step(self):
        #'id label classname entityname'
        fh = tempfile.NamedTemporaryFile(mode='w+')
        fh.write('foobarbaz')
        fh.flush()
        out_fh = tempfile.NamedTemporaryFile(suffix='.csv')
        src = (hortator.StoredEntity(1,'file', 'File', fh.name), )
        targets = (hortator.StoredEntity(1,'crc', 'CSVFile', out_fh.name), )
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
        self.cache = hortator.PersistentTaskGraph(self.cache_file.name, model, 
                                                 force_create=True)

    @unittest.skip('test not implemented')
    def test_StepGraph(self):
        raise NotImplementedError()
        #FIXME: where are the tests ?

    @unittest.skip('test not implemented')
    def test_addstep(self):
        raise NotImplementedError()
        #FIXME: where are the tests ?

class PersistanceTestCase(unittest.TestCase):

    cls_to_test = hortator.PersistentTaskGraph
    
    def setUp(self):
        self.cache_file = tempfile.NamedTemporaryFile(mode='w+')
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
                version = res.split(linesep)[0]
                return version
            def run(self, assets, parameters=('%D', )):
                # assets not used
                res = subprocess.check_output([self._execpath, '-c',
                                               "import %s; print(time.strftime('%D', time.localtime()))" % self._default_execpath])
                print(res)
        self.PythonTime = PythonTime

    def tearDown(self):
        self.cache_file.close()

    def test_PersistentTaskGraph(self):
        model = core.Model(tuple()) # getting away with an empty model for this step
        # test the initialization
        cache = self.cls_to_test(self.cache_file.name, model, force_create=True)
        # check that it was created
        self.assertTrue(cache.created)
        # empty set of steps
        s = tuple(cache.iter_steps())
        self.assertEqual(0, len(s))

    def test_statuslist(self):
        model = core.Model(tuple()) # getting away with an empty model for this step
        cache = self.cls_to_test(self.cache_file.name, model, force_create=True)
        statuslist = cache.statuslist
        self.assertEqual(set((y,x) for x,y in hortator._TASK_STATUS_LIST.items()), set(statuslist))

    def test_reopen(self):
        model = core.Model(tuple()) # getting away with an empty model for this step
        #create
        cache = self.cls_to_test(self.cache_file.name, model, force_create=True)
        #reopen
        cache_2 = self.cls_to_test(self.cache_file.name, model)
        # check that it was not created
        self.assertFalse(cache_2.created)
        # check that the statuslist is matching (as it also means that the inserts
        # to set up the DB were committed.
        self.assertEqual(cache.statuslist, cache_2.statuslist)

    def test_id_stored_entity(self):
        model = core.Model(tuple()) # getting away with an empty model for this step
        cache = self.cls_to_test(self.cache_file.name, rnaseq, force_create=True)
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
        cache = self.cls_to_test(self.cache_file.name, rnaseq, force_create=True)
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
        cache = self.cls_to_test(self.cache_file.name, rnaseq, force_create=True)
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
        cache = self.cls_to_test(self.cache_file.name, rnaseq, force_create=True)
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
        cache = self.cls_to_test(self.cache_file.name, rnaseq, force_create=True)
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
        cache = self.cls_to_test(self.cache_file.name, rnaseq, force_create=True)
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
        cache = self.cls_to_test(self.cache_file.name, model, force_create=True)
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
                version = res.split(linesep)[0]
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
        cache = self.cls_to_test(self.cache_file.name, model, force_create=True)
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
        cache = self.cls_to_test(self.cache_file.name, model, force_create=True)
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
        sources = railroadtracks.model.aligners.AssetsIndexer.Source(rnaseq.FASTAFile('foo.fasta'))

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

        # 1-element sources, several parameters
        db_id_2params = cache.id_stepconcrete(stepvariant_db_id.id,
                                              sources, targets, ("%Y", "Z"))
        db_id_same2params = cache.id_stepconcrete(stepvariant_db_id.id,
                                                  sources, targets, ("%Y", "Z"))
        self.assertEqual(db_id_2params.id, db_id_same2params.id)

        db_id_2otherparams = cache.id_stepconcrete(stepvariant_db_id.id,
                                                   sources, targets, ("%Y", "W"))
        self.assertNotEqual(db_id_2params.id, db_id_2otherparams.id)

        # 2-elements sources
        SrcCls = core.assetfactory('Source', 
                                   [core.AssetAttr('reference', rnaseq.FASTAFile, ''),
                                    core.AssetAttr('otherreference', rnaseq.FASTAFile, '')])
        sources = SrcCls(rnaseq.FASTAFile('foo.fasta'),
                         rnaseq.FASTAFile('bar.fasta'))
        db_id_notthesame = cache.id_stepconcrete(stepvariant_db_id.id,
                                                 sources, targets, parameters)
        self.assertNotEqual(db_id.id, db_id_notthesame.id)
        db_id_sameagain = cache.id_stepconcrete(stepvariant_db_id.id,
                                                sources, targets, parameters)
        self.assertEqual(db_id_sameagain.id, db_id_notthesame.id)


        # 1-element source / 1-element target
        sources = railroadtracks.model.aligners.AssetsIndexer.Source(rnaseq.FASTAFile('foo.fasta'))
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

        # fail if target assets are suddenly different
        targets_bar = railroadtracks.model.aligners.AssetsIndexer.Target(rnaseq.FilePattern('bar_idx'))
        self.assertRaises(ValueError, cache.id_stepconcrete, stepvariant_db_id.id,
                          sources, targets_bar, parameters)


    def test_get_assets(self):
        model = core.Model(tuple()) # getting away with an empty model for this step
        cache = self.cls_to_test(self.cache_file.name, model, force_create=True)
        class Activities(core.Enum):
            DATETIME = 'Give date/time'
        PythonTime = self.PythonTime
        python = PythonTime('python')
        stepvariant_db_id = cache.id_step_variant(python,
                                                  (Activities.DATETIME,))
        # 2-elements sources
        SrcCls = core.assetfactory('Source', 
                                   [core.AssetAttr('reference', rnaseq.FASTAFile, ''),
                                    core.AssetAttr('otherreference', rnaseq.FASTAFile, ''),
                                    core.AssetAttr('listoffiles', rnaseq.CSVFileSequence, '')])
        sources = SrcCls(rnaseq.FASTAFile('foo.fasta'),
                         rnaseq.FASTAFile('bar.fasta'),
                         rnaseq.CSVFileSequence((rnaseq.CSVFile('baz.csv'),rnaseq.CSVFile('baz2.csv'))))
        targets = core.AssetSet() # targets
        parameters = tuple()
        db_id = cache.id_stepconcrete(stepvariant_db_id.id,
                                      sources, targets, parameters)
        for storedthing in cache.get_srcassets(db_id.id):
            thing = storedthing.resurrect(rnaseq)

    def test_nconcrete_steps(self):
        model = core.Model(tuple()) # getting away with an empty model for this step
        cache = self.cls_to_test(self.cache_file.name, model, force_create=True)
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
        cache = self.cls_to_test(self.cache_file.name, model, force_create=True)
        split = Split(None)
        parameters = tuple()
        input_file = tempfile.NamedTemporaryFile(mode='w+')
        input_file.write('123')
        input_file.flush()
        head_file = tempfile.NamedTemporaryFile(mode='w+')
        tail_file = tempfile.NamedTemporaryFile(mode='w+')
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
        self.assertEqual(0, len(res))
        # create a new task
        task_id = cache.id_stepconcrete(stepvariant_db_id.id,
                                        sources, targets, parameters)
        # there is only one such 
        res = cache.get_targetsofactivity(ActivitiesSplit.HEADTAIL)
        self.assertEqual(2, len(res))
        #
        head_file = tempfile.NamedTemporaryFile(mode='w+')
        tail_file = tempfile.NamedTemporaryFile(mode='w+')
        targets_other = split.Assets.Target(core.File(head_file.name),
                                      core.File(tail_file.name))
        task_id_other = cache.id_stepconcrete(stepvariant_db_id.id,
                                              sources, targets_other, parameters=(1,))
        res = cache.get_targetsofactivity(ActivitiesSplit.HEADTAIL)
        self.assertEqual(4, len(res))
        # --
        cat = Cat(None)
        concat = tempfile.NamedTemporaryFile(mode='w+')
        sources = cat.Assets.Source(core.FileSequence((targets.head, targets_other.head)))
        targets = cat.Assets.Target(core.File(concat.name))
        # ensure that step variant is tracked
        stepvariant_db_id = cache.id_step_variant(cat,
                                                  cat.activities)
        task_id_cat = cache.id_stepconcrete(stepvariant_db_id.id,
                                            sources, targets, parameters)
        res = cache.get_targetsofactivity(ActivitiesSplit.MERGE)
        self.assertEqual(1, len(res))
                                        
    def test_get_parenttask_of_storedentity(self):
        model = core.Model(tuple()) # getting away with an empty model for this step
        cache = self.cls_to_test(self.cache_file.name, model, force_create=True)
        split = Split(None)
        parameters = tuple()
        input_file = tempfile.NamedTemporaryFile(mode='w+')
        input_file.write('123')
        input_file.flush()
        head_file = tempfile.NamedTemporaryFile(mode='w+')
        tail_file = tempfile.NamedTemporaryFile(mode='w+')
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
            self.assertEqual(task_id.id, res.id)
        #
        head_file_other = tempfile.NamedTemporaryFile(mode='w+')
        tail_file_other = tempfile.NamedTemporaryFile(mode='w+')
        targets_other = split.Assets.Target(core.File(head_file_other.name),
                                            core.File(tail_file_other.name))
        task_id_other = cache.id_stepconcrete(stepvariant_db_id.id,
                                              split.Assets.Source(targets.tail), targets_other, parameters=(1,))
        for sa in cache.get_srcassets(task_id_other):
            res = cache.get_parenttask_of_storedentity(sa)
            self.assertEqual(task_id.id, res.id)
        for sa in cache.get_targetassets(task_id_other):
            res = cache.get_parenttask_of_storedentity(sa)
            self.assertEqual(task_id_other.id, res.id)

        # --
        cat = Cat(None)
        sources_cat = cat.Assets.Source(core.FileSequence((targets.head, targets_other.head)))
        concat = tempfile.NamedTemporaryFile(mode='w+')
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
                    self.assertEqual(t.id, res.id)                    
            else:
                res = cache.get_parenttask_of_storedentity(sa)
                self.assertEqual(task_id_cat.id, res.id)

        head_file_other = tempfile.NamedTemporaryFile(mode='w+')
        tail_file_other = tempfile.NamedTemporaryFile(mode='w+')
        targets_other = split.Assets.Target(core.File(head_file_other.name),
                                            core.File(tail_file_other.name))
        task_id = cache.id_stepconcrete(stepvariant_db_id.id,
                                        split.Assets.Source(targets_cat.result), targets_other, tuple())
        for sa in cache.get_srcassets(task_id):
            res = cache.get_parenttask_of_storedentity(sa)
            self.assertEqual(task_id_cat.id, res.id)
        


    def test_gettargetsoftype(self):
        model = core.Model(tuple()) # getting away with an empty model for this step
        cache = self.cls_to_test(self.cache_file.name, model, force_create=True)
        class SplitCSV(Split):
            _name = 'split'
            activities = (core.DEFAULT_ACTIVITY.MISC, )
            version = '0.1'
            _default_execpath = 'None'
            class Assets(core.AssetsStep):
                Source = core.assetfactory('Source', [core.AssetAttr('file', railroadtracks.model.files.CSVFile, '')])
                Target = core.assetfactory('Target', [core.AssetAttr('head', railroadtracks.model.files.CSVFile, ''),
                                                      core.AssetAttr('tail', railroadtracks.model.files.CSVFile, '')])
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
        input_file = tempfile.NamedTemporaryFile(mode='w+')
        input_file.write('123')
        input_file.flush()
        head_file = tempfile.NamedTemporaryFile(mode='w+')
        tail_file = tempfile.NamedTemporaryFile(mode='w+')
        sources = split.Assets.Source(core.File(input_file.name))
        targets = split.Assets.Target(core.File(head_file.name),
                                      core.File(tail_file.name))
        stepvariant_db_id = cache.id_step_variant(split,
                                                  split.activities)
        res = cache.get_targetsoftype(core.File.__name__)
        self.assertEqual(0, len(res))
        task_id = cache.id_stepconcrete(stepvariant_db_id.id,
                                        sources, targets, parameters)
        res = cache.get_targetsoftype(core.File.__name__)
        self.assertEqual(2, len(res))

        splitcsv = SplitCSV(None)
        parameters = tuple()
        input_file = tempfile.NamedTemporaryFile(suffix=".csv", mode='w+')
        input_file.write('123')
        input_file.flush()
        head_file = tempfile.NamedTemporaryFile(suffix=".csv", mode='w+')
        tail_file = tempfile.NamedTemporaryFile(suffix=".csv", mode='w+')
        sources = splitcsv.Assets.Source(railroadtracks.model.files.CSVFile(input_file.name))
        targets = splitcsv.Assets.Target(railroadtracks.model.files.CSVFile(head_file.name),
                                         railroadtracks.model.files.CSVFile(tail_file.name))
        stepvariant_db_id = cache.id_step_variant(splitcsv,
                                                  splitcsv.activities)
        res = cache.get_targetsoftype(railroadtracks.model.files.CSVFile.__name__)
        self.assertEqual(0, len(res))
        task_id = cache.id_stepconcrete(stepvariant_db_id.id,
                                        sources, targets, parameters)
        res = cache.get_targetsoftype(railroadtracks.model.files.CSVFile.__name__)
        self.assertEqual(2, len(res))

        head_file = tempfile.NamedTemporaryFile(mode='w+')
        tail_file = tempfile.NamedTemporaryFile(mode='w+')
        targets = split.Assets.Target(core.File(head_file.name),
                                      core.File(tail_file.name))
        task_id_other = cache.id_stepconcrete(stepvariant_db_id.id,
                                              sources, targets, parameters=(1,))
        res = cache.get_targetsoftype(core.File.__name__)
        self.assertEqual(4, len(res))


class CachedPersistanceTestCase(PersistanceTestCase):

    cls_to_test = hortator.CachedPersistentTaskGraph
