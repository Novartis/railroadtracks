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

from collections import Counter
import tempfile, os, time, sys
import subprocess
import multiprocessing

import abc

from railroadtracks import hortator
from railroadtracks.easy.tasks import TaskSet, DbID

from six import with_metaclass

class ParallelExecutionAbstract(with_metaclass(abc.ABCMeta, object)):

    def map(self, taskset):
        raise NotImplementedError()
    def map_async(self, taskset):
        raise NotImplementedError()


def _execute_subprocess_in_dir(x):
    # expand 'x'
    taskid, cmd, wd = x
    import subprocess, tempfile, os, time
    od = os.path.abspath(os.curdir)
    t0 = None
    t1 = None
    env = os.environ.copy()
    if 'PYTHONPATH' in env:
        pythonpath = os.getenv('PYTHONPATH')
        env['PYTHONPATH'] = ':'.join((od, pythonpath))
    else:
        env['PYTHONPATH'] = od
    #env = {}
    try:
        os.chdir(wd)
        t0 = time.time()
        returncode = subprocess.check_call(cmd,
                                           env = env)
        t1 = time.time()
        return (taskid, returncode, t0, t1)
    except:
        return (taskid, 1, t0, t1)
    finally:
        os.chdir(od)

def _set_task_status(project, res):
    """
    :param res: iterable of (taskid, returncode) tuples.
    """
    ct = Counter()
    for x in res:
        try:
            taskid, returncode, t0, t1 = x
        except:
            raise ValueError('Unable to unpack values from:\n%s' % str(x))
        if returncode == 0:
            status = hortator._TASK_DONE
        else:
            status = hortator._TASK_FAILED
        if not hasattr(taskid, 'id'):
            taskid = DbID(taskid)
        project.persistent_graph.step_concrete_state(taskid,
                                                     hortator._TASK_STATUS_LIST[status])
        project.persistent_graph.set_time_t0(taskid, t0)
        project.persistent_graph.set_time_t1(taskid, t1)
        ct.update((status,))
    return ct

class DummyExecution(ParallelExecutionAbstract):
    def __init__(self):
        pass

    def map(self, taskset, tempdir=None):
        ct = Counter()
        for task in taskset:
            ct.update(task.status)
        return ct

class IterativeExecution(ParallelExecutionAbstract):
    """ Task mapper executing tasks in sequence and within the current process. """
    def __init__(self):
        pass

    def map(self, taskset, tempdir=None):
        assert isinstance(taskset, TaskSet)
        if len(taskset) == 0:
            return
        project = next(iter(taskset)).project
        od = os.curdir
        ct = Counter()
        for task in taskset:
            t1 = None
            wd = tempfile.mkdtemp(dir=tempdir)
            os.chdir(wd)
            t0 = time.time()
            try:
                task.execute()
                t1 = time.time()
                returncode = 0
            except:
                returncode = 1
            finally:
                pass
                os.chdir(od)
                res = (task.task_id, returncode, t0, t1)
                ct.update(_set_task_status(project, (res,)))

        return ct


class MultiprocessingExecution(ParallelExecutionAbstract):
    """ Task mapper using Python's multiprocessing """
    def __init__(self, processes, python_executable=sys.executable):
        """ 
        :param processes: number of simultaneous Python processing (underlying tools may use threads or processes themselves)
        """
        self.__processes = processes
        self._pyexecutable = python_executable

    def __get_ressources(self):
        self._pool = multiprocessing.Pool(self.__processes)
        
    def map(self, taskset, tempdir=None):
        assert isinstance(taskset, TaskSet)
        if len(taskset) == 0:
            return
        self.__get_ressources()
        pool = self._pool
        wds = (tempfile.mkdtemp(dir=tempdir) for x in taskset)
        ct = None
        res = None
        try:
            gen = ((x.task_id, x.unifex_cmd(python_executable=self._pyexecutable), y)
                   for x,y 
                   in zip(taskset, wds))
            res = pool.map(_execute_subprocess_in_dir, 
                           gen)
        finally:
            pool.close()
            pool.join()
            for d in wds:
                shutil.rmtree(d)
        if res is not None:
            project = next(iter(taskset)).project
            ct = _set_task_status(project, res)
        return ct

    def __str__(self):
        return "%s using %i processes." % (str(type(self)), self.__processes)


def select_mapper(taskset, mapper_smallset, mapper_largeset):
    if len(taskset) < 3:
        res = mapper_smallset
    else:
        res = mapper_largeset
    return res

class DuplexExecution(ParallelExecutionAbstract):
    """ """
    def __init__(self, mapper_a, mapper_b, func):
        self._mapper_a = mapper_a
        self._mapper_b = mapper_b
        self._func = func
    def map(self, taskset, tempdir=None):
        mapper = self._func(taskset, self._mapper_a, self._mapper_b)
        return mapper.map(taskset)

def exec_taskfilter(taskset):
    """ Returns a taskset in which tasks are candidates for execution. """
    # union of tasks that are either "TO-DO" or "FAILED"
    tasks_todo = taskset.filter_on_status(hortator._TASK_TODO)
    tasks_todo = tasks_todo.union(taskset.filter_on_status(hortator._TASK_FAILED))
    # only keep the tasks for which the parent task was successfully run
    tasks_todo = tasks_todo.filter_on_parent_status(hortator._TASK_DONE)
    return tasks_todo
