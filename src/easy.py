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

""" 
Module aiming at making common operations "easy".
"""
import logging
logger = logging.getLogger(__name__)

from collections import namedtuple, OrderedDict, Counter, defaultdict
import collections
from importlib import import_module
from railroadtracks import core, unifex, rnaseq, hortator
import os, shutil
import json
import tempfile
import itertools, operator
import subprocess
import multiprocessing
from abc import ABCMeta

_TASK_DONE = hortator._TASK_DONE
_TASK_TODO = hortator._TASK_TODO
_TASK_INPROGRESS = hortator._TASK_INPROGRESS
_TASK_FAILED = hortator._TASK_FAILED

DbID = namedtuple('DbID', 'id')


class Asset(object):
    """ An asset is either used by a task (then it is a source asset) or is produced
    by a task (then it is a target asset). """
    __slots__ = ('id', '_asset_id', 'entity', '_entity', 'project', '_project', 'parenttask')
    def __init__(self, project, savedentity, asset_id):
        self._asset_id = asset_id
        self._entity = savedentity
        self._project = project
    @property
    def id(self):
        return self._asset_id
    @property
    def entity(self):
        """ Instance of :class:`SavedEntityAbstract` 
        (or rather a child class thereof). """
        return self._entity
    @property
    def project(self):
        return self._project

    @property
    def parenttask(self):
        return self.project.todo._cache.get_parenttask_of_storedentity(self.id)
        
class Task(object):
    """ A task. """
    def __init__(self, project, call, task_id):
        self._call = call
        self._project = project
        if not hasattr(task_id, 'id'):
            task_id = DbID(task_id)
        assert isinstance(task_id.id, int), 'Task ID must be an integer.'
        self._task_id = task_id

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
        step_info = self._project.todo._cache.step_concrete_info(self._task_id)
        return step_info
    @property
    def project(self):
        """ Project in which the task is defined. """
        return self._project
    @property
    def dirname(self):
        """ Directory in which the targets of the task are saved. """
        dirname = os.path.join(self.project.wd,
                               self.project.todo.stepconcrete_dirname(self.task_id))
        return dirname

    def unifex_cmd(self):
        """ Parsed command to run the task. 

        To make a string to copy/paste into a shell script,
        use :func:`subprocess.list2cmdline`.

        """
        cmd = command_line(self._project, 
                           self._task_id.id, #FIXME: would work with _task_id ?
                           self._call.step, 
                           self._call.assets, 
                           self._call.parameters)
        return cmd

    def unifex_cmdline(self):
        """ 
        Command line as could be run in a shell.
        ( This is a wrapper around subprocess.list2cmdline(self.unifex_cmd()) )
        """
        return subprocess.list2cmdline(self.unifex_cmd())


    def execute(self):
        """ Execute the task. Note that the status of the task known to the project is not updated. """
        cmd,returncode = self._call.execute()
        if returncode != 0:
            raise unifex.UnifexError('The execution not successful with\n%s\nCheck/activate the logs to trace what is happening.' % str(self))

    def __str__(self):
        res = (super(Task, self).__str__(),
               'task ID: %i' % self.task_id.id,
               str(self.info),
               str(self.call),
               '---')
        return os.linesep.join(res)

    def parent_tasks(self):
        """ Tasks the current task is depending on. """
        res = TaskSet()
        # fetch the source assets
        entities = self.project.todo._cache.get_srcassets(self._task_id)
        for et in entities:
            # for each source asset, look for the tasks that also have it as a target
            if hasattr(et, 'iter_storedentities'):
                for sub_et in et.iter_storedentities():
                    sc_id = self.project.todo._cache.get_parenttask_of_storedentity(sub_et)
                    if sc_id is None:
                        # we reached a root node
                        continue
                    # only one task at the origin of a target
                    res.add(self.project.get_task(sc_id))
            else:
                sc_id = self.project.todo._cache.get_parenttask_of_storedentity(et)
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
        # fetch the target assets
        entities = self.project.todo._cache.get_targetassets(self._task_id.id)
        for et in entities:
            # for each source asset, look for the tasks that also have it as a source
            if hasattr(et, 'iter_storedentities'):
                for sub_et in et.iter_storedentities():
                    for sc_id in self.project.todo._cache.get_targetstepconcrete(sub_et):
                        if sc_id is None:
                            # root node
                            continue
                        res.add(self.project.get_task(sc_id))
            else:
                for sc_id in self.project.todo._cache.get_targetstepconcrete(et):
                    if sc_id is None:
                        # root node
                        continue
                    res.add(self.project.get_task(sc_id))
        return res

    def all_child_tasks(self):
        """ Get all the child tasks, direct or indirect. """
        res = TaskSet()
        nodes_visited = set()
        task_queue = collections.deque()
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
        return self._project.todo._cache.step_concrete_status(self._task_id)
        
    @status.setter
    def status(self, value):
        self._project.todo._cache.step_concrete_state(self._task_id,
                                                      hortator._TASK_STATUS_LIST[value])

    @property
    def time_points(self):
        return self._project.todo._cache.task_time_points(self._task_id)

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
        ct = Counter(x.status[0][1] for x in self._d.values())
        return ct
    def split_by_status(self):
        res = defaultdict(TaskSet)
        for task in self:
            #FIXME: simplify access to status
            res[task.status[0][1]].add(task)
        return res
    def filter_on_status(self, status):
        res = TaskSet()
        for task in self:
            #FIXME: simplify access to status
            if task.status[0][1] == status:
                res.add(task)
        return res
    def filter_on_parent_status(self, status):
        res = TaskSet()
        for task in self:
            ok = True
            for parent_task in task.parent_tasks():
                if parent_task.status[0][1] != status:
                    ok = False
                    break
            if ok:
                res.add(parent_task)
        return res


class ParallelExecutionAbstract(object):
    __metaclass__ = ABCMeta

    def map(self, taskset):
        raise NotImplementedError()
    def map_async(self, taskset):
        raise NotImplementedError()


def _execute_subprocess_in_dir(x):
    # 'x' is a 2-tuple:
    taskid, cmd, wd = x
    import subprocess, tempfile, os, time
    od = os.curdir
    t0 = None
    t1 = None
    try:
        os.chdir(wd)
        t0 = time.time()
        returncode = subprocess.check_call(cmd)
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
        print('%s: %s' % (str(taskid), str(status)))
        project.todo._cache.step_concrete_state(taskid,
                                                hortator._TASK_STATUS_LIST[status])
        project.todo._cache.set_time_t0(taskid, t0)
        project.todo._cache.set_time_t1(taskid, t1)
        ct.update((status,))
    return ct

class DummyExecution(ParallelExecutionAbstract):
    def __init__(self):
        pass

    def map(self, taskset, tempdir=None):
        ct = Counter()
        for task in taskset:
            ct.update(task.status[0][1])
        return ct

class MultiprocessingExecution(ParallelExecutionAbstract):
    def __init__(self, processes):
        """ 
        :param processes: number of simultaneous Python processing (underlying tools may use threads or processes themselves)
        """
        self.__processes = processes

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
        try:
            res = pool.map(_execute_subprocess_in_dir, ((x.task_id, x.unifex_cmd(), y) for x,y in zip(taskset, wds)))
            pool.close()
            project = next(iter(taskset)).project
            ct = _set_task_status(project, res)
        finally:
            for d in wds:
                shutil.rmtree(d)
        return ct

#FIXME: remove ?
Task2 = namedtuple('Step', 'assets run info step_concrete_id')


ActivityCount = namedtuple('ActivityCount', 'count name status')
class ActivityDoneCount(object):
    __slots__ = ['count_done', 'total', 'name']
    def __init__(self, name, count_done=0, total=0):
        self.count_done = count_done
        self.total = total
        self.name = name


def dict_project_view_storage(project):
    if project.newproject:
        status = 'New "railroadtracks" project'
    else:
        status = 'Reopened existing "railroadtracks" project'
    statvfs= os.statvfs(project.wd)

    avail = statvfs.f_bavail * statvfs.f_frsize / (1024 ** 3)
    total = statvfs.f_blocks * statvfs.f_frsize / (1024 ** 3)
    used = (statvfs.f_blocks - statvfs.f_bfree) * statvfs.f_frsize / (1024 ** 3)

    taskstatuscount = project.todo._cache.nconcrete_steps_status
    totaltasks = sum(x.count for x in taskstatuscount)

    d = {'status': status, 
         'wd': project.wd, 
         'wdspace_avail': avail,
         'wdspace_used': used,
         'wdspace_total': total,
         'db_fn': project.db_fn,
         'db_fn_size': '{:,}'.format(os.stat(project.db_fn).st_size / (1024 ** 2)),
         'taskstatuscount': taskstatuscount,
         'totaltasks': totaltasks}
    return d

def dict_project_view_activities(project, done=hortator._TASK_DONE):
    cursor = project.todo._cache.connection.cursor()
    sql = """
SELECT count(step_activity_id), step_activity.label, ss.label
FROM step_activity
LEFT JOIN step_type2activity AS sta ON sta.step_activity_id=step_activity.id
LEFT JOIN step_type ON sta.step_type_id=step_type.id
LEFT JOIN step_variant AS sv ON sv.step_type_id=step_type.id
LEFT JOIN step_concrete AS sc ON sc.step_variant_id=sv.id
LEFT JOIN step_status AS ss ON ss.id=sc.step_status_id
GROUP BY step_activity_id, step_status_id
ORDER BY step_activity_id
    """
    cursor.execute(sql)
    counts = tuple(ActivityCount(*x) for x in cursor.fetchall())
    activities = list()
    for name, grp in itertools.groupby(counts, operator.attrgetter('name')):
        adc = ActivityDoneCount(name)
        for ac in grp:
            if ac.status == done:
                adc.count_done = ac.count
            adc.total += ac.count
        activities.append(adc)
    d = {'activities': activities }
    return d



def dict_project_view_results(project):
    results = project.todo._cache.iter_finaltargets()
    ct = dict()
    for result in results:
        if result.status_label == _TASK_DONE:
            count_done = 1
        else:
            count_done = 0
        adc = ct.get(result.stored_entity_classname)
        if adc is None:
            adc = ActivityDoneCount(result.stored_entity_classname, 
                                    count_done=count_done,
                                    total = 1)
            ct[result.stored_entity_classname] = adc
        else:
            adc.count_done += count_done
            adc.total += 1
    res = {'results': ct.values()}
    return res

def dict_assetset_view(assetset):
    l = len(assetset)
    AssetSetElement = namedtuple('AssetSetElement', 'isdefined name cls')
    d = {'super_repr': super(core.AssetSet, assetset).__repr__(),
         'elements': []}
    for i in range(l):
        isdefined = getattr(assetset, 
                            assetset._sources[i].name)._defined
        d['elements'].append(AssetSetElement(isdefined,
                                             assetset._sources[i].name,
                                             str(assetset._sources[i].cls)))
    return d

def str_progressbar(val, maxval, width=20, fillchar="#",
                    emptychar=' '):
    assert val <= maxval
    nfill = int(round(1.0 * val / maxval * width))
    res = fillchar * nfill
    if emptychar is not None:
        res += emptychar * (width-nfill)
    return res

def str_project_view_storage(project, width=10):
    d = dict_project_view_storage(project)
    for k in ('wdspace_avail', 'wdspace_total'):
        d[k] = '{:,.2f}'.format(d[k])
    res = """Status: %(status)s
Working directory: %(wd)s
Storage
  Available: %(wdspace_avail)s GB
  Total:     %(wdspace_total)s GB
Tasks
  Total:     %(totaltasks)i
""" % d
    return res


def str_project_view_activities(project, done=hortator._TASK_DONE, width=20):
    d = dict_project_view_activities(project, done=done)
    if len(d['activities']) == 0:
        name_len = 0
    else:
        name_len = max(len(a.name) for a in d['activities'])
    name_len = max(name_len, len('Activity')) # column header
    res = list(('%s    %s' %('Activity'.center(name_len),'Progress'),))
    for a in d['activities']:
        res.append('%s    |%s| %.2f%% (%s/%s)' % (a.name.ljust(name_len), 
                                                  str_progressbar(a.count_done, a.total, width=width),
                                                  100.0*a.count_done/a.total,
                                                  '{:,}'.format(a.count_done),
                                                  '{:,}'.format(a.total))
               )
    res.append('')
    return os.linesep.join(res)

def str_project_view_results(project, width=20):
    d = dict_project_view_results(project)
    if len(d['results']) == 0:
        name_len = 0
    else:
        name_len = max(len(r.name) for r in d['results'])
    name_len = max(name_len, len('Result Type')) # column header
    res = list(('%s    %s' %('Result Type'.center(name_len),'Progress'),))
    for r in d['results']:
        res.append('%s    |%s| %.2f%%' % (r.name.ljust(name_len), 
                                          str_progressbar(r.count_done, r.total, width=width),
                                          100.0*r.count_done/r.total))
    res.append('')
    return os.linesep.join(res)
    return res

def str_project_view(project):
    res = tuple((x(project) for x in (str_project_view_storage,
                                      str_project_view_activities,
                                      str_project_view_results)))
    return ('--%s' % os.linesep).join(res)


class Project(object):
    """ A project, that is a directory containing data as well as a persistant storage of steps and how derived data and final results were obtained. """

    def __init__(self, model, wd='railroadtracks_project', db_fn=None, force_create = False):
        """
        :param wd: Name of a working directory (where all intermediate and final results will be saved.
        :type wd: :class:`str`
        :param db_fn: Name of a file name for the database. If None, use the file "railroadtracks.db" in the directory specified in "wd".
        :type db_fn: :class:`str` or None
        :rtype: a tuple with (:class:`hortator.StepGraph`, working directory as a :class:`str`, file name for the database and a :class:`str`)
        """

        # check the model
        if not hasattr(model, 'ACTIVITY'):
            raise ValueError('The parameter `model` does not appear to be a model (ACTIVITY missing).')

        newproject = force_create
        # derived data files and final results will be in a temporary directory
        if not os.path.isdir(wd):
            raise ValueError("The working directory '%s' should be a directory" % wd)

        # create the database
        if db_fn is None:
            db_fn = os.path.join(wd, 'railroadtracks.db')

        if not os.path.exists(db_fn):
            newproject = True
        cache = hortator.PersistentTaskList(db_fn, model, wd=wd, force_create=force_create)
        # create the dependency graph
        todo = hortator.StepGraph(cache)
        self._todo = todo
        self._wd = wd
        self._db_fn = db_fn
        self._newproject = newproject

    @property
    def todo(self):
        return self._todo

    @property
    def wd(self):
        """Working directory. The directory in which derived data is stored."""
        return self._wd

    @property
    def db_fn(self):
        """Path to the database file"""
        return self._db_fn

    @property
    def newproject(self):
        """Tell whether this is a new project, rather than a existing project (re)opened. """
        return self._newproject

    @property
    def model(self):
        """Model used in the project."""
        return self.todo._cache._model

    def __repr__(self):
        res = """
%s
Working directory: %s
Database file: %s
Number of recorded steps: %i
"""
        if self.newproject:
            status = 'New project'
        else:
            status = 'Reopened existing project'
        return res % (status, self.wd, self.db_fn,
                      self.todo._cache.nconcrete_steps)

    def __str__(self):
        return str_project_view(self)


    def get_targetsofactivity(self, activity):
        """ 
        Retrieve the targets of steps performing a specific activity.
        (calls the method of the same name in the contained :class:`PersistentTaskList`)
        :param activity: an activity
        :type activity: :class:`Enum`
        """
        return self.todo._cache.get_targetsofactivity(activity)

    def get_targetsoftype(self, obj):
        """ 
        Retrieve the targets having a given type.
        (calls the method of the same name in the contained :class:`PersistentTaskList`)
        :param obj: a type, an instance, or a type name
        :type obj: :class:`type`, :class:`Object`, or :class:`str`
        """
        if isinstance(obj, str):
            # name of a class
            cls = getattr(self.todo._cache._model, obj)
        elif isinstance(obj, core.SavedEntityAbstract):
            cls = type(obj)
        elif issubclass(obj, core.SavedEntityAbstract):
            cls = obj
        else:
            raise ValueError("Invalid parameter.")
        return self.todo._cache.get_targetsoftype(cls.__name__)

    def iter_srcassets(self, task):
        """ Return the source files for a given task
        (calls the method of the same name in the contained :class:`PersistentTaskList`)
        :param task: a :class:`Task`.
        """
        for x in self.todo._cache.get_srcassets(task.task_id.id):
            yield Asset(self, x.resurrect(self.model), DbID(x.id))
            

    def iter_targetassets(self, task):
        """ Return the target files for a given task
        (calls the method of the same name in the contained :class:`PersistentTaskList`)
        :param task: a :class:`Task`.
        """
        for x in self.todo._cache.get_targetassets(task.task_id.id):
            yield Asset(self, x.resurrect(self.model), DbID(x.id))

    def add_task(self, step, assets, parameters=(), tag = 1):
        """
        Add a task to the project. If any of the assets' targets is not defined,
        it will be defined automatically.

        :param step:
        :type param: :class:`StepAbstract`
        :param assets:
        :param parameters:
        :param tag: a tag (to differentiate repetitions of the exact same task)
        """
        step_concrete_id = self.todo.add(step, assets, parameters=parameters, tag=tag)
        call = unifex.Call(step, assets, parameters)
        task = Task(self, call, step_concrete_id)
        return task


    def get_task(self, task_id):
        """
        Given a task ID, retrieve the associated task.
        """
        
        if not hasattr(task_id, 'id'):
            task_id = DbID(task_id)
        assert isinstance(task_id.id, int)
        task_dbentry = self.todo._cache._get_stepconcrete(task_id)
        # get the step
        if '.' in task_dbentry.clsname:
            module_name, cls_name = task_dbentry.clsname.rsplit('.', 1)
        else:
            module_name = self.todo.cache._model
            cls_name = task_dbentry.clsname
        module = import_module(module_name)
        cls = getattr(module, cls_name)
        step = cls(task_dbentry.executable)
        # get the assets
        sources = [None, ]*len(cls.Assets.Source._fields)
        model = self.todo._cache._model
        for i, x in enumerate(self.todo._cache.get_srcassets(task_id.id)):
            src_i = cls.Assets.Source._fields.index(x.label)
            sources[src_i] = x.resurrect(model)
        targets = [None, ]*len(cls.Assets.Target._fields)
        for i, x in enumerate(self._todo._cache.get_targetassets(task_id.id)):
            tgt_i = cls.Assets.Target._fields.index(x.label)
            targets[tgt_i] = x.resurrect(model)
        assets = cls.Assets(cls.Assets.Source(*sources),
                            cls.Assets.Target(*targets))
        # get the parameters
        parameters = json.loads(task_dbentry.parameters)
        # build the call
        call = unifex.Call(step, assets, parameters)
        # build the task
        task = Task(self, call, task_dbentry.id)
        return task

    def add_asset(self, asset):
        if hasattr(asset, 'iter_storedentities'):
            raise NotImplementedError('Not yet implemented')
            self.todo._cache.id_stored_sequence()
        else:
            seid = self.todo._cache.id_stored_entity(type(asset), asset.name)
            return hortator.StoredEntity(seid.id, None, asset.__class__.__name__, asset.name)

def call_factory(project, step_concrete_id, stepobj, assets, parameters=()):
    """ 
    :param stepobj: Step object
    :type stepobj: instance of :class:`StepAbstract`
    :param assets: assets
    :type assets: instance of :class:`AssetsStep`
    :rtype: :class:`function`
    """
    def call():
        import hortator
        cmd, returncode = stepobj.run(assets, parameters=parameters)
        if returncode != 0:
            raise RuntimeError('Error (return code %i) while running: %s' % (returncode, ' '.join(cmd)))
        project.todo._cache.step_concrete_state(step_concrete_id,
                                                hortator._TASK_STATUS_LIST[hortator._TASK_DONE])
    return call

def command_line(project, stepconcrete_id, stepobj, assets, parameters=()):
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
    cmd = [sys.executable, '-m', 'railroadtracks.unifex', 'run', model]
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
        project.todo.provenancewalk_storedentity(stored_entity,
                                                 func_stored_entity,
                                                 func_stored_sequence,
                                                 func_step_concrete_factory(res),
                                                 func_storedentity_stepconcrete,
                                                 func_stepconcrete_storedentity)
        all_stepconcrete.append(res)
    return all_stepconcrete


def _delete_tasks_in_set(taskset):
    """ Delete a set of tasks. """
    # get the target assets associated with the task
    #db = task.project.todo._cache
    #storedentitities = list()
    #for task in taskset:
    #    storedentities.extend(db.get_targetassets(task.id))

    sql_create_task_tmp = """
    CREATE TEMP TABLE task_to_delete (
      id INTEGER PRIMARY KEY
    )"""

    sql_add_task_id_to_delete = """
    INSERT INTO task_to_delete
    (id)
    VALUES
    (?)
    """

    sql_create_stored_entity_tmp = """
    CREATE TEMP TABLE se_to_delete
    AS
      SELECT stored_entity.id AS se_id
      FROM stored_entity 
      INNER JOIN step_concrete2targetfile
      ON stored_entity.id=step_concrete2targetfile.stored_entity_id
      INNER JOIN task_to_delete
      ON task_to_delete.id=step_concrete2targetfile.step_concrete_id
      WHERE se_id IN (SELECT task_to_delete.id FROM task_to_delete)"""

    sql_delete_src_associations = """
    DELETE FROM 
    step_concrete2targetfile
    WHERE step_concrete_id IN (SELECT id from task_to_delete)
    """

    sql_delete_target_association = """
    DELETE FROM 
    step_concrete2targetfile
    WHERE step_concrete_id IN (SELECT id from task_to_delete)
    """

    sql_delete_src_association = """
    DELETE FROM 
    step_concrete2srcfile
    WHERE step_concrete_id IN (SELECT id from task_to_delete)
    """

    sql_delete_stored_entry = """
    DELETE FROM 
    stored_entity
    WHERE id IN (SELECT id from se_to_delete)
    """

    sql_delete_tasks = """
    DELETE FROM 
    step_concrete
    WHERE id IN (SELECT id from task_to_delete)
    """

    cursor = taskset._project.todo._cache.connection.cursor()
    cursor.execute(sql_create_task_tmp)
    cursor.executemany(sql_add_task_id_to_delete, ((task.task_id.id,) for task in taskset))
    cursor.execute(sql_create_stored_entity_tmp)
    # delete the associations between stored_entities that are the target of tasks in the taskset and tasks
    cursor.execute(sql_delete_target_association)
    # same for sources
    cursor.execute(sql_delete_src_association)
    # delete the stored_entries 
    cursor.execute(sql_delete_stored_entry)
    # delete the tasks in the taskset
    cursor.execute(sql_delete_tasks)
    # delete the directories
    for task in taskset:
        shutil.rmtree(task.dirname)
    cursor.execute("DROP TABLE task_to_delete")
    cursor.execute("DROP TABLE se_to_delete")
    taskset._project.todo._cache.connection.commit()


class FrozenNamespace(object):
    # Implementation note: a :class:`namedtuple` was not used, because we wanted
    # all class attributes and methods to start with '_', in order to let autocompletion
    # only show instance-defined attributes by default
    def __init__(self, args):
        """ 
        :param args: iterable of (key,value) pairs, the key being the name in the namespace.
        """
        fields = []
        for key, value in args:
            fields.append(key)
            setattr(self, key, value)
        self.__fields = tuple(fields)
    def __setattr__(self, key, value):
        if hasattr(self, key):
            raise ValueError("The attribute '%s' is already defined." % key)
        self.__dict__[key] = value
    @property
    def _fields(self):
        return self.__fields


class Environment(object):
    """
    Represent the current environment in a (presumably) easy way for writing
    recipes.
    """
    def __init__(self, model):
        """ 
        :param model: Python module following the :mod:`railroadtracks` model. 
        """
        classes, activityenum = core.steplist(model), model.ACTIVITY
        stepclasses = list()
        stepinstances = list()
        knownclassnames = set()
        knowninstancenames = set()
        AvailableStepGroup = namedtuple('AvailableStepGroup', 'default_instances classes')
        for cls in unifex._make_stepdict(model).values():
            if not issubclass(cls, core.StepAbstract):
                raise ValueError("Classes used as steps must inherit from core.StepAbstract.")
            stepname = cls._name
            if stepname is None:
                raise ValueError('The step name for class "%s" is not defined.' % cls)
            elif stepname in knownclassnames:
                # check unique class names
                raise ValueError('The step name "%s" is defined twice.' % stepname)
            elif stepname.lower() in knowninstancenames:
                # check unique lower-case for class names (these will be used for instance names)
                raise ValueError('The step name "%s", when converted to lower case, is defined twice.' % stepname)
            stepclasses.append(cls)
            stepinstances.append(cls(cls._default_execpath))
            knownclassnames.add(stepname)
            knowninstancenames.add(stepname.lower())
        
        stepclasses = tuple(stepclasses)
        self._stepclasses = stepclasses
        stepinstances = tuple(stepinstances)
        self._stepinstances = stepinstances

        od = OrderedDict((x, []) for x in activityenum)
        
        for s,si in zip(stepclasses, stepinstances):
            for a in s.activities:
                od[a].append((s, si))

        activitylist = list()
        for activity, classinstancepairs in od.items():
            # iterate over the keys (each being an activity)
            clsnames = list()
            instancenames = list()
            for stepcls, stepinstance in classinstancepairs:
                clsname = stepcls.__name__
                instancename = clsname.lower()
                clsnames.append(clsname)
                instancenames.append(instancename)
            #NT_cls = frozencollection(activity.name+'_cls', instancenames)
            activitylist.append(FrozenNamespace((key, element[1]) for key, element in zip(instancenames, classinstancepairs)))
            #activitylist.append(NT_cls(*(element for element in clsnames)))

        self.__activities = FrozenNamespace((key.name, value) for key, value in zip(od.keys(), activitylist))


    @property
    def activities(self):
        """ Access the activities declared by the model as attributes. """
        return self.__activities

    @property
    def stepclasses(self):
        """ Steps. """
        return self._stepclasses

    @property
    def stepinstances(self):
        """ Default instance for the steps (created from the default executables in the PATH). """
        return self._stepinstances

if __name__ == '__main__':

    def _exec_info(args):
        model = importlib.import_module(args.model)
        project = Project(model, wd=args.wd, db_fn=args.db_fn)
        dbid = DbID(args.task_id)
        res = project.todo._cache.step_concrete_info(dbid)
        print(res)

    def _exec_cmd(args):
        project = Project(wd=args.wd, db_dn=args.db_fn)
        raise NotImplementedError()

    def _exec_run(args):
        project = Project(wd=args.wd, db_dn=args.db_fn)
        task = project.get_task(args.task_id)
        if task.status == hortator._TASK_DONE:
            print('The task %i is already done' % task.id)
            return
        else:
            try:
                res = task.execute()
                task.status = hortator._TASK_DONE
            except:
                task.status = hortator._TASK_FAILED

    import argparse, importlib
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help='Action to perform.')

    parser_info = subparsers.add_parser('status', help="Retrieve information about a step.")
    parser_cmd = subparsers.add_parser('cmd', help="Return the unified execution command line.")
    parser_run = subparsers.add_parser('run', help="Run a step.")

    parser.add_argument('-d', '--working-directory',
                        dest = 'wd',
                        help = 'Working directory for the project.')
    parser.add_argument('-b', '--db-file',
                        dest = 'db_fn',
                        help = 'Database file')
    parser.add_argument('-m', '--model',
                        dest = 'model',
                        help = 'Model (as a Python module)')    
    parser.add_argument(dest = 'task_id',
                        type = int,
                        help = 'Task ID')
    parser_info.set_defaults(func = _exec_info)
    parser_cmd.set_defaults(func = _exec_cmd)
    parser_run.set_defaults(func = _exec_run)


    args = parser.parse_args()
    args.func(args)
