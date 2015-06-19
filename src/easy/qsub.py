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

"""
Extension to rnaseqrs to use SGE (no luck with DRMAA, so this is command-line based.
"""
from __future__ import print_function

import os
import re
import csv
import stat
import subprocess
import tempfile
from railroadtracks import hortator, easy, environment, core

import logging
logger = logging.getLogger(__name__)

QSUB_EXEC = 'qsub'

QSUB_JOB_PAT = re.compile('Job ([0-9]+)\\.([0-9]+) exited with exit code ([0-9]+)\\.')

SCRIPT_TEMPLATE = """
#!/bin/bash
#$ -V
%(qsub_preamble)s

# setup
#source %(setup)s

#export PYTHONPATH=%(pythonpath)s:${PYTHONPATH}

# TASKSET_IDS
IFS=', ' read -a TASKSET_IDS_ARRAY <<< ${TASKSET_IDS}
RRT_TASK_ID=${TASKSET_IDS_ARRAY[$((SGE_TASK_ID-1))]}

mkdir -p ${RRT_EXECUTION_DIR}/${RRT_TASK_ID}
cd ${RRT_EXECUTION_DIR}/${RRT_TASK_ID}
python -c 'from railroadtracks.easy import Project; import %(model)s as model; project=Project(model, wd="'${RRT_PROJECT_WD}'"); task=project.get_task('${RRT_TASK_ID}'); task.execute()'
"""


class QsubExecution(object):
    def __init__(self, setup, executiondir, qsub_preamble=''):
        """
        :param setup: BASH code to run before each task is run (setting environment variables common to each task can be done here)
        :type setup: :class:`str`
        :param executiondir: Parent directory where the tasks will run (Hadwrap will create its directories there)
        :param qsub_preamble: Preamble to add to the script called by qsub. This is the place to add qsub configuration statements.
        """
        # Check the presence of "qsub"
        assert environment.Executable.ispresent(QSUB_EXEC)     
        assert os.path.exists(executiondir) and os.path.isdir(executiondir)
        self._executiondir = executiondir
        self.setup = setup
        self.qsub_preamble = qsub_preamble

    def map(self, taskset, tempdir=None):
        assert isinstance(taskset, easy.TaskSet)
        if len(taskset) == 0:
            return
        mapdir = tempfile.mkdtemp(dir=self._executiondir)
        script_path = os.path.join(mapdir, 'script.sh')
        with open(script_path, 'w') as fh_out:
            script = SCRIPT_TEMPLATE % {'qsub_preamble': self.qsub_preamble,
                                        'setup': self.setup,
                                        'pythonpath': os.getenv('PYTHONPATH'),
                                        'model': taskset._project.model.__name__}
            fh_out.write(script)
        
        ts_range = '1-%i' % len(taskset)

        qsub_exec_abs = subprocess.check_output(('which', QSUB_EXEC)).rstrip()
        cmd = (qsub_exec_abs,
               '-t', ts_range,
               '-sync', 'y',
               '-S', '/bin/bash',
               '-o', mapdir,
               '-e', mapdir,
               script_path)
        env = {'RRT_EXECUTION_DIR': mapdir,
               'RRT_PROJECT_WD': taskset._project.wd,
               'TASKSET_IDS': ', '.join(str(x.task_id.id) for x in taskset)}
        for var in ('SGE_CELL',
                    'SGE_EXECD_PORT',
                    'SGE_QMASTER_PORT',
                    'SGE_ROOT',
                    'SGE_CLUSTER_NAME'):
            env[var] = os.getenv(var)

        wd_orig = os.getcwd()
        try:
            os.chdir(mapdir)
            proc = subprocess.Popen(subprocess.list2cmdline(cmd),
                                    env = env,
                                    stdout = subprocess.PIPE,
                                    bufsize = 1,
                                    shell=True)
        finally:
            os.chdir(wd_orig)
        project = taskset._project
        for output in proc.stdout:
            m = QSUB_JOB_PAT.match(output)
            if m is None:
                continue
            jobid, sge_taskid, returncode = m.groups()
            rrt_taskids = tuple(x.task_id for x in taskset)
            if returncode == '0':
                project.get_task(rrt_taskids[int(sge_taskid)]).status = hortator._TASK_DONE
            else:
                project.get_task(rrt_taskids[int(sge_taskid)]).status = hortator._TASK_FAILED
        stdout_value = proc.communicate()[0]



