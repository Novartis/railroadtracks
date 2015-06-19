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


import tempfile, shutil, os
from collections import namedtuple
import unittest
from railroadtracks import rnaseq, easy
from railroadtracks.easy import tasksetgraph
import railroadtracks.model.simulate
PHAGEFASTA = railroadtracks.model.simulate.PHAGEFASTA

class TaskSetGraphTestCase(unittest.TestCase):
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


    def testAddTaskSet(self):
        ts = easy.TaskSet()
        ts.add(self.task)
        #
        tsg = tasksetgraph.TaskSetGraph()
        tsindex_id = tsg.add_taskset(ts)
        dict_retrieve = tsg.digraph.node[tasksetgraph.get_tasksetid((self.task,))]
        ts_retrieve = dict_retrieve['taskset']
        self.assertEqual(len(ts), len(ts_retrieve))
        self.assertEqual(tuple(x.task_id for x in ts), tuple(x.task_id for x in ts_retrieve))

        # add connected tasks
        taskset_align = easy.TaskSet()
        aligner = rnaseq.Bowtie2()
        for fq_name in ('foo.fq', 'bar.fq'):
            Assets = aligner.Assets
            assets = Assets(Assets.Source(self.task.call.assets.target.indexfilepattern,
                                          rnaseq.FASTQPossiblyGzipCompressed(fq_name),
                                          None), None)
            task_align = self.project.add_task(aligner, assets)
            taskset_align.add(task_align)
        tsalign_id = tsg.add_taskset(taskset_align)
        lst = tsg.execution_list()
        self.assertEqual(tsindex_id, tasksetgraph.get_tasksetid(lst[0].taskset))
        self.assertEqual(tsalign_id, tasksetgraph.get_tasksetid(lst[1].taskset))

    def testAddTaskSetDuplicateTask(self):
        ts = easy.TaskSet()
        ts.add(self.task)
        tsg = tasksetgraph.TaskSetGraph()
        tsg.add_taskset(ts)

        ts2 = easy.TaskSet()
        ts2.add(self.task)
        self.assertRaises(ValueError, tsg.add_taskset, ts2)

    def testAddTask(self):
        tsg = tasksetgraph.TaskSetGraph()
        tsg.add(self.task)
        dict_retrieve = tsg.digraph.node[tasksetgraph.get_tasksetid((self.task,))]
        ts_retrieve = dict_retrieve['taskset']
        self.assertEqual(1, len(ts_retrieve))
        self.assertEqual((self.task.task_id,), tuple(x.task_id for x in ts_retrieve))


    def testAddTaskDifferentProject(self):
        project2 = easy.Project(rnaseq, self.wd2)
        bowtie2index = rnaseq.Bowtie2Build()
        reference_fn = PHAGEFASTA
        Assets = bowtie2index.Assets
        task2 = project2.add_task(bowtie2index, 
                                 Assets(Assets.Source(rnaseq.FASTAFile(reference_fn))))
        #
        tsg = tasksetgraph.TaskSetGraph()
        tsg.add(self.task)
        self.assertRaises(ValueError, tsg.add, task2)

    def testExecutionList(self):
        tsg = tasksetgraph.TaskSetGraph()
        # initial task/taskset is building a bowtie2 index
        tsindex_id = tsg.add(self.task)
        lst = tsg.execution_list()
        self.assertEqual(1, len(lst))

        # add a second task set that performs the alignment of several pairs of
        # FASTQ files
        n_fastq = 3
        pairnames = self._create_fastq(n_fastq, self.wd2)

        index_task = self.task
        bowtie2 = rnaseq.Bowtie2()
        Assets = bowtie2.Assets
        # Class to model FASTQ files that are optionally compressed
        FASTQ = rnaseq.FASTQPossiblyGzipCompressed
        ts = easy.TaskSet()
        # note that we are included the pair of FASTQ already aligned above
        for read1_fn, read2_fn in pairnames:
            task = self.project.add_task(bowtie2,
                                         Assets(Assets.Source(index_task.call.assets.target.indexfilepattern,
                                                              FASTQ(read1_fn),
                                                              FASTQ(read2_fn))))
            ts.add(task)
        tsalign_id = tsg.add(ts)
        lst = tsg.execution_list()

    def testExecution(self):
        tsg = tasksetgraph.TaskSetGraph()
        tsg.add(self.task)
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
        tsg.add(ts)
        # no task mapper, check that an exception is raised
        self.assertRaises(ValueError, tsg.execute)

        # set the default task mapper to a multiprocessing based one
        n_processes = 1
        mpe = easy.MultiprocessingExecution(n_processes)
        tsg.defaultmapper = mpe
        tsg.execute()

        # use a run-only-once filter (that is only try to execute
        # tasks with the status "TO DO"
        flt = lambda taskset: taskset.filter_on_status(easy._TASK_TODO)
        tsg.defaulttasksetfilter = flt
        lst = tsg.execution_list()
        for elt in lst:
            ts_f = tsg._filtered_taskset(elt)
            # check that nothing is left to be exectuted
            # (since everything has been run earlier)
            self.assertEqual(0, len(ts_f))
    
