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


import os, subprocess, sys, shutil
from railroadtracks import hortator, unifex
from collections import namedtuple, OrderedDict, deque, Counter, defaultdict
import itertools

_TASK_DONE = hortator._TASK_DONE
_TASK_TODO = hortator._TASK_TODO
_TASK_INPROGRESS = hortator._TASK_INPROGRESS
_TASK_FAILED = hortator._TASK_FAILED

DbID = namedtuple('DbID', 'id')

def command_line(project, stepconcrete_id, stepobj, assets, 
                 parameters=(),
                 python_executable=sys.executable,
                 logging_file=None,
                 logging_level=None):
    """ 
    :param stepobj: Step object
    :type stepobj: instance of :class:`StepAbstract`
    :param assets: assets
    :type assets: instance of :class:`AssetsStep`
    :rtype: a named :class:`tuple` with the arguments for a command line call
    """
    import sys
    executable = stepobj._execpath
    model = type(stepobj)._name
    # FIXME: duplicated elsewhere in the same module - make a function
    cmd = [python_executable, '-m', 'railroadtracks.unifex', 'run', model]
    if executable is not None:
        cmd.append(executable)
    cmd.extend(['--module', stepobj.__class__.__module__])
    cmd.append('-s')
    field_items = list()
    for field,src in zip(assets.source._fields, assets.source):
        if src is None:
            continue
        for item in src:
            field_items.append((field,item))
    cmd.extend('%s=%s' % x for x in field_items)
    cmd.append('-t')
    cmd.extend('%s=%s' % (field,item) for field,src in zip(assets.target._fields, assets.target) for item in src)
    if logging_file is not None:
        cmd.append('--logging-file=%s' % logging_file)
    if logging_level is not None:
        cmd.append('--logging-level=%s' % logging_level)
    if len(parameters) > 0:
        cmd.append('-p')
        cmd.extend(' %s' % x for x in parameters)
    return cmd


#FIXME: rename this function
def get_steplist_up(project, stored_entities):
    def func_stored_entity(stored_entity):
        pass
    def func_stored_sequence(stored_entity):
        pass
    def func_step_concrete_factory(res):
        def f(step_concrete):
            res.append(step_concrete)
        return f
    def func_storedentity_stepconcrete(stored_entity, step_concrete):
        pass
    def func_stepconcrete_storedentity(step_concrete, stored_entity):
        pass
    all_stepconcrete = list()
    for stored_entity in stored_entities:
        res = list()
        project.cache.provenancewalk_storedentity(stored_entity,
                                                  func_stored_entity,
                                                  func_stored_sequence,
                                                  func_step_concrete_factory(res),
                                                  func_storedentity_stepconcrete,
                                                  func_stepconcrete_storedentity)
        all_stepconcrete.append(res)
    return all_stepconcrete



class Task(object):
    """ A task. """
    def __init__(self, project, call, task_id):
        self._call = call
        self._project = project
        if hasattr(task_id, 'id'):
            self._task_id = task_id.id
            self.is_new = task_id.new
        else:
            self._task_id = task_id
            self.is_new = None

    @property
    def activities(self):
        return self.call.step.activities

    @property
    def assets(self):
        return self.call.assets
    
    @property
    def call(self):
        """ A :class:`railroadtracks.unifex.Call` object. """
        return self._call
    @property
    def task_id(self):
        """ ID of the task in the project. """
        return self._task_id
    @property
    def info(self):
        """ Status information for the task. """
        # wrap into a DbID
        taskid = DbID(self._task_id)
        step_info = self._project.persistent_graph.step_concrete_info(taskid)
        return step_info
    @property
    def project(self):
        """ Project in which the task is defined. """
        return self._project
    @property
    def dirname(self):
        """ Directory in which the targets of the task are saved. """
        dirname = os.path.join(self.project.wd,
                               self.project.cache.stepconcrete_dirname(DbID(self.task_id)))
        return dirname

    def reset(self, ignore_errors=False):
        """ 'Reset' a task by setting the status to "TO DO" and erasing all result files in the task directory."""
        dirname = self.dirname
        self.status = _TASK_TODO
        for entry in os.listdir(dirname):
            if os.path.isdir(entry):
                shutil.rmtree(os.path.join(dirname, entry), ignore_errors=ignore_errors)
            else:
                os.unlink(os.path.join(dirname, entry))
        
    def unifex_cmd(self, python_executable=sys.executable,
                   logging_file=None,
                   logging_level=None):
        """ Parsed command to run the task. 

        To make a string to copy/paste into a shell script,
        use :func:`subprocess.list2cmdline`.

        """
        cmd = command_line(self._project, 
                           self._task_id, #FIXME: would work with _task_id ?
                           self._call.step, 
                           self._call.assets, 
                           parameters=self._call.parameters,
                           python_executable=python_executable,
                           logging_file=logging_file,
                           logging_level=logging_level)

        return cmd

    def unifex_cmdline(self, python_executable=sys.executable):
        """ 
        Command line as could be run in a shell.
        ( This is a wrapper around subprocess.list2cmdline(self.unifex_cmd()) )
        """
        return subprocess.list2cmdline(self.unifex_cmd(python_executable=python_executable))


    def execute(self):
        """ Execute the task. Note that the status of the task known to the project is not updated. """
        cmd,returncode = self._call.execute()
        if returncode != 0:
            raise unifex.UnifexError('The execution not successful with\n%s\nCheck/activate the logs to trace what is happening.' % str(self))

    def __str__(self):
        res = (super(Task, self).__str__(),
               'task ID: %i' % self.task_id,
               str(self.info),
               str(self.call),
               '---')
        return os.linesep.join(res)


    def parent_tasks(self):
        """ Tasks the current task is depending on. """
        res = TaskSet()
        # wrap into a DbID
        taskid = DbID(self._task_id)
        # fetch the source assets
        entities = self.project.persistent_graph.get_srcassets(taskid)
        for et in entities:
            # for each source asset, look for the tasks that also have it as a target
            if hasattr(et, 'iter_storedentities'):
                for sub_et in et.iter_storedentities():
                    sc_id = self.project.persistent_graph.get_parenttask_of_storedentity(sub_et)
                    if sc_id is None:
                        # we reached a root node
                        continue
                    # only one task at the origin of a target
                    res.add(self.project.get_task(sc_id))
            else:
                sc_id = self.project.persistent_graph.get_parenttask_of_storedentity(et)
                if sc_id is None:
                    # we reached a root node
                    continue
                # only one task at the origin of a target
                res.add(self.project.get_task(sc_id))
        return res

    def primordial_tasks(self):
        """ Direct or indirect parent tasks with source assets that do not
        have parent tasks themselves.
        In other words, root nodes in the dependency graph connected to this task.

        This method is useful to identify the raw data results are derived from.
        """
        res = TaskSet()
        # start with the targets because the current task may also have 
        # primordial source assets
        targets = tuple(self.iter_targetassets())
        up = get_steplist_up(self.project, [x.id for x in targets])
        for tsk in (y for x in up for y in x):
            ancestor_task = self.project.get_task(tsk)
            # The source assets, that is the assets used by the task, can be extracted
            for src_asset in ancestor_task.iter_sourceassets():
                parent_task = src_asset.parenttask
                if parent_task is None:
                    # This source asset does not have a parent task in the project
                    # It is therefore considered a root asset in railroadtracks
                    # (the asset might originate from earlier steps, 
                    # but they are not stored in the project)
                    res.add(ancestor_task)
                    # The task is in, no need to iterate further
                    break
        return res

    def child_tasks(self):
        """ Tasks depending on the current task. """
        res = TaskSet()
        # wrap into a DbID
        taskid = DbID(self._task_id)
        # fetch the target assets
        entities = self.project.persistent_graph.get_targetassets(taskid)
        for et in entities:
            # for each source asset, look for the tasks that also have it as a source
            if hasattr(et, 'iter_storedentities'):
                for sub_et in et.iter_storedentities():
                    for sc_id in self.project.persistent_graph.get_targetstepconcrete(sub_et):
                        if sc_id is None:
                            # root node
                            continue
                        res.add(self.project.get_task(sc_id))
            else:
                for sc_id in self.project.persistent_graph.get_targetstepconcrete(et):
                    if sc_id is None:
                        # root node
                        continue
                    res.add(self.project.get_task(sc_id))
        return res

    def all_child_tasks(self):
        """ Get all the child tasks, direct or indirect. """
        res = TaskSet()
        nodes_visited = set()
        task_queue = deque()
        # get the tasks downstream
        for t in self.child_tasks():
            if t.task_id in nodes_visited:
                continue
            task_queue.append(t)
        while len(task_queue) > 0:
            task = task_queue.pop()
            if task in res:
                continue
            nodes_visited.add(task.task_id)
            res.add(task)
            # get the tasks downstream
            for t in task.child_tasks():
                if t.task_id in nodes_visited:
                    continue
                task_queue.append(t)
        return res

    @property
    def status(self):
        """ Get/Set the status of the task in the project.
        The value must be a valid task status.
        """
        # wrap into a DbID
        taskid = DbID(self._task_id)
        res = self._project.persistent_graph.step_concrete_status(taskid)
        #FIXME: this will likely change at the lower level, the assertion will
        #       act as watcher
        assert(isinstance(res, list))
        assert(len(res) == 1)
        assert(len(res[0]) == 2)
        return res[0][1]
        
    @status.setter
    def status(self, value):
        # wrap into a DbID
        taskid = DbID(self._task_id)
        self._project.persistent_graph.step_concrete_state(taskid,
                                                           hortator._TASK_STATUS_LIST[value])

    @property
    def time_points(self):
        # wrap into a DbID
        taskid = DbID(self._task_id)
        return self._project.persistent_graph.task_time_points(taskid)

    def iter_sourceassets(self):
        return self.project.iter_srcassets(self)

    def iter_targetassets(self):
        return self.project.iter_targetassets(self)

    
class TaskSet(object):
    """ Ordered set of tasks.
    This class can be used for independent tasks for which the order of completion does not matter,
    and provide a container for tasks that can run in parallel.
    """
    _d = None
    _s = None
    def __init__(self, project=None, iterable=tuple()):
        self._d = OrderedDict()
        self._s = set()
        self._project = project
        for item in iterable:
            self.add(item)

    def _ensure_task_project(self, task):
        if self._project is None:
            self._project = task._project
        elif task._project != self._project:
            raise ValueError('Mismatching project.')
    def _ensure_taskset_project(self, taskset):
        if taskset._project is None:
            if len(taskset) > 0:
                raise Exception('We should not be here.')
            else:
                return
        if self._project is None:
            self._project = taskset._project
        elif taskset._project != self._project:
            raise ValueError('Mismatching project.')
    def __iter__(self):
        return iter(self._d.values())

    def iter_chunks(self, maxsize):
        """ Return an iterator of TaskSet of maximum size `maxsize`. """
        it = iter(self)
        l = len(self)
        for chunk_i in range(0, l, maxsize):
            if (chunk_i+maxsize) > l:
                chunksize = l - chunk_i
            else:
                chunksize = maxsize
            yield TaskSet(iterable=itertools.islice(it, 0, chunksize))
            
    def __contains__(self, task):
        self._ensure_task_project(task)
        return task.task_id in self._d
    def __len__(self):
        return len(self._d)
    def __str__(self):
        res = [super(TaskSet, self).__str__(),
               '  with %i element(s) of respective status:' % len(self)]
        res.extend('  - %s: %i' % (k, v) for k,v in self.status().items())
        return '\n'.join(res)
    def add(self, task):
        """ Add a task. """
        assert isinstance(task, Task), 'Element should be in an instance of Task.'
        self._ensure_task_project(task)
        self._s.add(task.task_id)
        self._d[task.task_id] = task
    def remove(self, task):
        """ remove a task. """
        assert isinstance(task, Task), 'Element should be in an instance of Task.'
        self._ensure_task_project(task)
        self._s.remove(task.task_id)
        del(self._d[task.task_id])
    def union(self, taskset):
        assert isinstance(taskset, TaskSet), "'taskset' should be a instance of TaskSet."
        self._ensure_taskset_project(taskset)
        taskids = self._s.union(taskset._s)
        u = TaskSet()
        for task_id in taskids:
            try:
                task = self._d[task_id]
            except KeyError:
                task = taskset._d[task_id]
            u.add(task)
        return u
    def update(self, taskset):
        assert isinstance(taskset, TaskSet), "'taskset' should be a instance of TaskSet."
        self._ensure_taskset_project(taskset)
        self._s.update(taskset._s)
        for task_id in taskset._s:
            self._d.update(taskset._d)

    def intersection(self, taskset):
        assert isinstance(taskset, TaskSet), "'taskset' should be a instance of TaskSet."
        self._ensure_taskset_project(taskset)
        taskids = self._s.intersection(taskset._s)
        i = TaskSet()
        for task_id in taskids:
            task = self._d[task_id]
            i.add(task)
        return i
    def difference(self, taskset):
        assert isinstance(taskset, TaskSet), "'taskset' should be a instance of TaskSet."
        self._ensure_taskset_project(taskset)
        taskids = self._s.difference(taskset._s)
        i = TaskSet()
        for task_id in taskids:
            task = self._d[task_id]
            i.add(task)
        return i
    def isdisjoint(self, taskset):
        assert isinstance(taskset, TaskSet), "'taskset' should be a instance of TaskSet."
        self._ensure_taskset_project(taskset)
        res = self._s.isdisjoint(taskset._s)
        return res
    def issubset(self, taskset):
        assert isinstance(taskset, TaskSet), "'taskset' should be a instance of TaskSet."
        self._ensure_taskset_project(taskset)
        res = self._s.issubset(taskset._s)
        return res
    def issuperset(self, taskset):
        assert isinstance(taskset, TaskSet), "'taskset' should be a instance of TaskSet."
        self._ensure_taskset_project(taskset)
        res = self._s.issuperset(taskset._s)
        return res
    def values(self):
        """ Return the tasks as a tuple (Python 2) or an iterator (Python 3). """
        return self._d.values()
    def status(self):
        """ Return a `Counter` of status labels. """
        ct = Counter(x.status for x in self._d.values())
        return ct
    def split_by_status(self):
        res = defaultdict(TaskSet)
        for task in self:
            #FIXME: simplify access to status
            res[task.status].add(task)
        return res
    def filter_on_status(self, status):
        res = TaskSet()
        for task in self:
            #FIXME: simplify access to status
            if task.status == status:
                res.add(task)
        return res
    def filter_on_parent_status(self, status):
        res = TaskSet()
        for task in self:
            ok = True
            for parent_task in task.parent_tasks():
                if parent_task.status != status:
                    ok = False
                    break
            if ok:
                res.add(task)
        return res

