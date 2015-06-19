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
import shutil, os, tempfile
from railroadtracks import environment
from railroadtracks.easy import qsub
from railroadtracks.easy import Project, rnaseq, Environment, TaskSet
from railroadtracks.model.simulate import PHAGEFASTA

SETUP = 'source %s' % os.path.abspath('setup.sh')
class QsubTestCase(unittest.TestCase):    
    def setUp(self):
        wd = tempfile.mkdtemp(dir=os.path.expanduser("~"))
        self.qsub_dir = tempfile.mkdtemp(dir=os.path.expanduser("~"))
        self.project = Project(rnaseq, wd=wd)

    def tearDown(self):
        # more complicated than it should: NFS and/or qsub seem to be leaving read-only locks behind
        wd = self.project.wd
        del(self.project)
        shutil.rmtree(wd)
        shutil.rmtree(self.qsub_dir)

    @unittest.skipIf(not environment.Executable.ispresent(qsub.QSUB_EXEC),
                     "No executable 'qsub' in the PATH.")
    def test_map(self):
        env = Environment(rnaseq)
        b2b = env.activities.INDEX.bowtie2build
        Assets = b2b.Assets
        assets = Assets(Assets.Source(rnaseq.FASTAFile(PHAGEFASTA)),
                        None)
        ts = TaskSet()
        task = self.project.add_task(b2b, assets)
        ts.add(task)
        task = self.project.add_task(b2b, assets, ('--ntoa',))
        ts.add(task)
        h2 = qsub.QsubExecution(SETUP, self.qsub_dir)
        h2.map(ts)
        # check that all tasks are complete
        for task in ts:
            self.assertEquals(hortator._TASK_DONE, task.status)


if __name__ == '__main__':
    unittest.main()

    #suite = unittest.TestLoader.loadTestsFromTestCase(QsubTestCase)
