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
Package aiming at making common operations "easy".
"""

import logging
logger = logging.getLogger(__name__)

from collections import namedtuple, OrderedDict, Counter, defaultdict
import collections
from importlib import import_module
from railroadtracks import core, unifex, rnaseq, hortator
import os, shutil, time
import json
import tempfile
import itertools, operator

from railroadtracks.easy.tasks import Task, TaskSet, DbID
from railroadtracks.easy.tasks import (_TASK_DONE,
                                       _TASK_TODO,
                                       _TASK_INPROGRESS,
                                       _TASK_FAILED)

from railroadtracks.easy.execution import (DummyExecution,
                                           IterativeExecution,
                                           MultiprocessingExecution)

TASK_DONE = _TASK_DONE
TASK_TODO = _TASK_TODO
TASK_INPROGRESS = _TASK_INPROGRESS
TASK_FAILED = _TASK_FAILED

class Asset(object):
    """ An asset is either used by a task (then it is a source asset) or is produced
    by a task (then it is a target asset). """
    __slots__ = ('_asset_id', '_entity', '_project')
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
        return self.project.persistent_graph.get_parenttask_of_storedentity(self.id)
        

ActivityCount = namedtuple('ActivityCount', 'count name status')
class ActivityDoneCount(object):
    __slots__ = ['count_done', 'count_failed', 'total', 'name']
    def __init__(self, name, count_done=0, count_failed=0, total=0):
        self.count_done = count_done
        self.count_failed = count_failed
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

    taskstatuscount = project.persistent_graph.nconcrete_steps_status
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

def dict_project_view_activities(project, done=hortator._TASK_DONE, failed=hortator._TASK_FAILED):
    cursor = project.persistent_graph.connection.cursor()
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
            elif ac.status == failed:
                adc.count_failed = ac.count
            adc.total += ac.count
        activities.append(adc)
    d = {'activities': activities }
    return d



def dict_project_view_results(project):
    results = project.persistent_graph.iter_finaltargets()
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


def str_project_view_activities(project, 
                                done=hortator._TASK_DONE, failed=hortator._TASK_FAILED,
                                width=20):
    d = dict_project_view_activities(project, done=done, failed=failed)
    if len(d['activities']) == 0:
        name_len = 0
    else:
        name_len = max(len(a.name) for a in d['activities'])
    name_len = max(name_len, len('Activity')) # column header
    res = list(('%s    %s' %('Activity'.center(name_len),'Progress'),))
    for a in d['activities']:
        res.append('%s    |%s| %.2f%% (%s/%s - failed: %s)' % (a.name.ljust(name_len), 
                                                  str_progressbar(a.count_done, a.total, width=width),
                                                  100.0*a.count_done/a.total,
                                                  '{:,}'.format(a.count_done),
                                                  '{:,}'.format(a.total),
                                                  '{:,}'.format(a.count_failed))
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
    """ A project, that is a directory containing data as well as a persistent storage of steps and how derived data and final results were obtained. """

    def __init__(self, model, wd='railroadtracks_project', db_fn=None, force_create = False, 
                 isolation_level=None, cached=False):
        """
        :param model: module with the definition of the `model` (that is the classes defining the type of tasks computed)
        :param wd: Name of a working directory (where all intermediate and final results will be saved.
        :type wd: :class:`str`
        :param db_fn: Name of a file name for the database. If None, use the file "railroadtracks.db" in the directory specified in "wd".
        :type db_fn: :class:`str` or None
        :param isolation_level: passed to the constructor for :class:`PersistentTaskGraph`
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
        if cached:
            cls_persistent = hortator.CachedPersistentTaskGraph
        else:
            cls_persistent = hortator.PersistentTaskGraph

        persistent_graph = cls_persistent(db_fn, model, wd=wd, force_create=force_create,
                                          isolation_level=isolation_level)

        self._persistent_graph = persistent_graph
        # cache the dependency graph
        cache = hortator.StepGraph(persistent_graph)
        self._cache = cache
        self._wd = wd
        self._db_fn = db_fn
        self._newproject = newproject

    @property
    def persistent_graph(self):
        """ Access to the persistent graph. This is will be slower than access the cached version, but will
        always be accurate."""
        return self._persistent_graph

    @property
    def cache(self):
        """ Cache for the dependency graph. This may not reflect changes made by an independent access
        to the project or a direct access to the persistent graph (see :attr:`persistent_graph`) but
        will be much faster and the preferred way to access the graph when there is only one concurrent access.
        """
        return self._cache

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
        return self.persistent_graph._model

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
                      self.persistent_graph.nconcrete_steps)

    def __str__(self):
        return str_project_view(self)


    def get_targetsofactivity(self, activity):
        """ 
        Retrieve the targets of steps performing a specific activity.
        (calls the method of the same name in the contained :class:`PersistentTaskGraph`)
        :param activity: an activity
        :type activity: :class:`Enum`
        """
        return self.persistent_graph.get_targetsofactivity(activity)

    def get_targetsoftype(self, obj):
        """ 
        Retrieve the targets having a given type.
        (calls the method of the same name in the contained :class:`PersistentTaskGraph`)
        :param obj: a type, an instance, or a type name
        :type obj: :class:`type`, :class:`Object`, or :class:`str`
        """
        if isinstance(obj, str):
            # name of a class
            cls = getattr(self.persistent_graph._model, obj)
        elif isinstance(obj, core.SavedEntityAbstract):
            cls = type(obj)
        elif issubclass(obj, core.SavedEntityAbstract):
            cls = obj
        else:
            raise ValueError("Invalid parameter.")
        return self.persistent_graph.get_targetsoftype(cls.__name__)

    def get_tasksoftype(self, taskorstep):
        """
        Retrieve tasks matching the type of "taskorstep".

        The parameter 'taskorstep' can be either:
        - a class inheriting from StepAbstract
        - an instance of Task
        - an instance of StepAbstract
        In the two last cases, the path to the executable and the version number
        are used to filter the tasks returned.
        """
        cache = self.persistent_graph
        cursor = cache.connection.cursor()

        stepvariant_ids = None
        if issubclass(taskorstep, core.StepAbstract):
            cls = taskorstep
            step_type_id = cache.id_step_type(cls.activities)
            sql = """
            SELECT id
            FROM step_variant
            WHERE step_type_id=?
            AND cls=?
            """
            #FIXME: building the 'cls' string the SQL string
            #       should be factored out as a function
            #       in hortator.py
            cursor.execute(sql, (step_type_id.id, '.'.join((cls.__module__,
                                                         cls.__name__))))
            stepvariant_ids = tuple(x[0] for x in cursor.fetchall())
        else:
            if isinstance(taskorstep, Task):
                step = taskorstep.call.step
            elif isinstance(taskorstep, core.StepAbstract):
                step = taskorstep
            stepvariant_ids = (cache.id_step_variant(step, step.activities).id, )

        if stepvariant_ids is None:
            raise ValueError('"taskorstep" can be either a Task or a core.StepAbstract')

        taskset = TaskSet()

        sql = """
        SELECT step_concrete.id
        FROM step_concrete
        WHERE step_concrete.step_variant_id=?
        """

        for sv_id in stepvariant_ids:
            cursor.execute(sql, (sv_id,))
            for res in cursor.fetchall():
                task_id = res[0]
                taskset.add(self.get_task(task_id))

        return taskset

    def get_taskswithactivity(self, activity):
        """
        Retrieve tasks matching the given activity.
        """
        cache = self.persistent_graph
        cursor = cache.connection.cursor()

        taskset = TaskSet()

        sql = """
        SELECT step_concrete.id
        FROM step_concrete
        INNER JOIN step_variant
        ON step_variant.id=step_concrete.step_variant_id
        INNER JOIN step_type
        ON step_type.id=step_variant.step_type_id
        INNER JOIN step_type2activity
        ON step_type2activity.step_type_id=step_type.id
        INNER JOIN step_activity
        ON step_activity.id=step_type2activity.step_activity_id
        WHERE step_activity.label=?
        """

        cursor.execute(sql, (activity.value,))
        for res in cursor.fetchall():
            task_id = res[0]
            taskset.add(self.get_task(task_id))

        return taskset

    def iter_srcassets(self, task):
        """ Return the source files for a given task
        (calls the method of the same name in the contained :class:`PersistentTaskGraph`)
        :param task: a :class:`Task`.
        """
        for x in self.persistent_graph.get_srcassets(task.task_id):
            yield Asset(self, x.resurrect(self.model), DbID(x.id))
            

    def iter_targetassets(self, task):
        """ Return the target files for a given task
        (calls the method of the same name in the contained :class:`PersistentTaskGraph`)
        :param task: a :class:`Task`.
        """
        for x in self.persistent_graph.get_targetassets(task.task_id):
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
        step_concrete_id = self.cache.add(step, assets, parameters=tuple(parameters), tag=tag)
        call = unifex.Call(step, assets, parameters)
        task = Task(self, call, step_concrete_id)
        return task


    def get_task(self, task_id):
        """
        Given a task ID, retrieve the associated task.
        """
        
        if not hasattr(task_id, 'id'):
            task_id = DbID(task_id)

        task_dbentry = self.persistent_graph._get_stepconcrete(task_id)
        # get the step
        if '.' in task_dbentry.clsname:
            module_name, cls_name = task_dbentry.clsname.rsplit('.', 1)
        else:
            module_name = self.persistent_graph._model
            cls_name = task_dbentry.clsname
        module = import_module(module_name)
        cls = getattr(module, cls_name)
        step = cls(task_dbentry.executable)
        # get the assets
        sources = [None, ]*len(cls.Assets.Source._fields)
        model = self.persistent_graph._model
        for i, x in enumerate(self.persistent_graph.get_srcassets(task_id)):
            src_i = cls.Assets.Source._fields.index(x.label)
            sources[src_i] = x.resurrect(model)
        targets = [None, ]*len(cls.Assets.Target._fields)
        for i, x in enumerate(self._persistent_graph.get_targetassets(task_id)):
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
            self.persistent_graph.id_stored_sequence()
        else:
            seid = self.persistent_graph.id_stored_entity(type(asset), asset.name)
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
        project.persistent_graph.step_concrete_state(step_concrete_id,
                                                hortator._TASK_STATUS_LIST[hortator._TASK_DONE])
    return call




def _delete_tasks_in_set(taskset):
    """ Delete a set of tasks. """
    # get the target assets associated with the task
    #db = task.project.persistent_graph
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

    cursor = taskset._project.persistent_graph.connection.cursor()
    cursor.execute(sql_create_task_tmp)
    cursor.executemany(sql_add_task_id_to_delete, ((task.task_id,) for task in taskset))
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
    taskset._project.persistent_graph.connection.commit()


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
            stepname = cls.__name__.lower()
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
        stepinstances = tuple(stepinstances)

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
        self._stepclasses = stepclasses
        self._stepinstances = stepinstances

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
        res = project.persistent_graph.step_concrete_info(dbid)
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

