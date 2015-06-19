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

import tempfile, shutil, os, subprocess
from collections import namedtuple
import unittest

import railroadtracks
from railroadtracks import environment, rnaseq, easy

import railroadtracks.model.simulate
PHAGEFASTA = railroadtracks.model.simulate.PHAGEFASTA

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
        assets = bowtie2.Assets(bowtie2.Assets.Source(rnaseq.FASTAFile(PHAGEFASTA)),
                                bowtie2.Assets.Target(rnaseq.SavedBowtie2Index(os.path.join(self.wd, 'foo'))))
        parameters = tuple()
        step_concrete_id = project.cache.add(bowtie2, assets, parameters = parameters)
        cmd = easy.tasks.command_line(project, step_concrete_id, bowtie2, assets, parameters = parameters)
        returncode = subprocess.check_call(cmd)
        self.assertEqual(0, returncode)

    def _test_EasyProject(self, cached):
        packagedir = os.path.dirname(railroadtracks.__file__)
        PHAGEFASTA = os.path.join(packagedir, 'EF204940.FASTA')
        bowtie2index_exec = environment.Executable('bowtie2-build')
        project = easy.Project(rnaseq, self.wd, cached=cached)
        # Initial number of steps is zero
        self.assertEqual(0, project.persistent_graph.nconcrete_steps)
        bowtie2index = rnaseq.Bowtie2Build(bowtie2index_exec.path)
        reference_fn = PHAGEFASTA
        task_index = project.add_task(bowtie2index, 
                                      bowtie2index.Assets(bowtie2index.Assets.Source(rnaseq.FASTAFile(reference_fn))))
        # Number of steps is one
        self.assertEqual(1, project.persistent_graph.nconcrete_steps)
        task_index_same = project.add_task(bowtie2index, 
                                           bowtie2index.Assets(bowtie2index.Assets.Source(rnaseq.FASTAFile(reference_fn))))
        # Number of steps should remain one (same assets and paremeters)
        self.assertEqual(1, project.persistent_graph.nconcrete_steps)

    @unittest.skipIf(not environment.Executable.ispresent('bowtie2-build'),
                     'bowtie2-build is not in the PATH')
    def test_EasyProject_NoCache(self):
        self._test_EasyProject(False)

    @unittest.skipIf(not environment.Executable.ispresent('bowtie2-build'),
                     'bowtie2-build is not in the PATH')
    def test_EasyProject_Cache(self):
        self._test_EasyProject(True)

    @unittest.skipIf(not environment.Executable.ispresent('bowtie2-build'),
                     'bowtie2-build is not in the PATH')
    def test_Project_get_tasksoftype(self):
        packagedir = os.path.dirname(railroadtracks.__file__)
        PHAGEFASTA = os.path.join(packagedir, 'EF204940.FASTA')
        bowtie2index_exec = environment.Executable('bowtie2-build')
        project = easy.Project(rnaseq, self.wd)
        bowtie2index = rnaseq.Bowtie2Build(bowtie2index_exec.path)
        reference_fn = PHAGEFASTA

        for obj in (rnaseq.Bowtie2Build, bowtie2index):
            taskset = project.get_tasksoftype(obj)
            self.assertEqual(0, len(taskset))
        
        task_index = project.add_task(bowtie2index, 
                                      bowtie2index.Assets(bowtie2index.Assets.Source(rnaseq.FASTAFile(reference_fn))))
        
        for obj in (rnaseq.Bowtie2Build, bowtie2index, task_index):
            taskset = project.get_tasksoftype(obj)
            self.assertEqual(1, len(taskset))



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
                                      bowtie2index.Assets(bowtie2index.Assets.Source(rnaseq.FASTAFile(reference_fn))))

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
        s1 = set(x.task_id for x in alignment_tasks)
        s2 = set(x.task_id for x in res)
        self.assertEqual(0, len(s1.difference(s2)))
        self.assertEqual(0, len(s2.difference(s1)))

        # find the parents of an alignment
        align_task = next(iter(alignment_tasks))
        parent_taskset = align_task.parent_tasks()
        # only one parent (the indexing step)
        self.assertEqual(1, len(parent_taskset))
        parent_task = next(iter(parent_taskset))
        self.assertEqual(bowtie2index.activities, parent_task.call.step.activities)
        self.assertEqual(bowtie2index._name, parent_task.call.step._name)

        # split by status
        # (first create a TaskSet)
        alignment_taskset = easy.TaskSet()
        for t in alignment_tasks:
            alignment_taskset.add(t)
        taskset_groups = alignment_taskset.split_by_status()
        self.assertEqual(1, len(taskset_groups)) # only one status here

        # root tasks
        root_tasks = align_task.primordial_tasks()
        # check that the indexing step is indeed one them
        self.assertTrue(task_index.task_id in (x.task_id for x in root_tasks))

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
                                      bowtie2index.Assets(bowtie2index.Assets.Source(rnaseq.FASTAFile(reference_fn))))

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
        s1 = set(x.task_id for x in alignment_tasks)
        s2 = set(x.task_id for x in res)
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



class IterativeExecutionTestCase(unittest.TestCase):
    def setUp(self):
        wd = tempfile.mkdtemp()
        self.wd2 = tempfile.mkdtemp()
        self.project = easy.Project(rnaseq, wd)
        bowtie2index = rnaseq.Bowtie2Build()
        reference_fn = PHAGEFASTA
        Assets = bowtie2index.Assets
        task = self.project.add_task(bowtie2index, 
                                     Assets(Assets.Source(rnaseq.FASTAFile(reference_fn))))
        self.task = task
        self.bowtie2index = bowtie2index

    def tearDown(self):
        shutil.rmtree(self.project.wd)
        shutil.rmtree(self.wd2)

    def _create_fastq(self, n, d):
        reference_fn = PHAGEFASTA
        with open(reference_fn) as fasta_fh:
            reference = next(railroadtracks.model.simulate.readfasta_iter(fasta_fh))

        ReadPair = namedtuple('ReadPair', 'read1 read2')
        pairnames = tuple(ReadPair(os.path.join(d, 'A_%i_1.fq' % p_i), 
                                   os.path.join(d, 'A_%i_2.fq' % p_i)) for p_i in range(1, n+1))
        for read1_fn, read2_fn in pairnames:
            with open(read1_fn, 'wb') as read1_fh, \
                 open(read2_fn, 'wb') as read2_fh:
                railroadtracks.model.simulate.randomPEreads(read1_fh, read2_fh, reference)
        return pairnames

    def testMap(self):

        mapper = easy.IterativeExecution()

        ts = easy.TaskSet()
        ts.add(self.task)

        mapper.map(ts)

        n_fastq = 3
        pairnames = self._create_fastq(n_fastq, self.wd2)

        index_task = self.task
        bowtie2 = rnaseq.Bowtie2()
        Assets = bowtie2.Assets
        FASTQ = rnaseq.FASTQPossiblyGzipCompressed
        ts = easy.TaskSet()

        for read1_fn, read2_fn in pairnames:
            task = self.project.add_task(bowtie2,
                                         Assets(Assets.Source(index_task.call.assets.target.indexfilepattern,
                                                              FASTQ(read1_fn),
                                                              FASTQ(read2_fn))))
            ts.add(task)
        self.assertEqual(3, ts.status()[railroadtracks.easy._TASK_TODO])
        mapper.map(ts)
        self.assertEqual(3, ts.status()[railroadtracks.easy._TASK_DONE])


    def testIter_chunk(self):
        ts = easy.TaskSet()
        ts.add(self.task)
        Assets = self.bowtie2index.Assets
        reference_fn = PHAGEFASTA
        for tag in range(2, 11):
            task = self.project.add_task(self.bowtie2index, 
                                         Assets(Assets.Source(rnaseq.FASTAFile(reference_fn))),
                                         tag=tag)
            ts.add(task)

        self.assertEqual(10, len(ts))
        chunks = tuple(ts.iter_chunks(3))
        for i in range(3):
            self.assertEqual(3, len(chunks[i]))
        self.assertEqual(1, len(chunks[3]))
