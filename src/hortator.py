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
Persistence/memoization for the DAG
"""

from collections import namedtuple, Counter
from importlib import import_module
from railroadtracks.core import StepAbstract
from railroadtracks import __version__
import collections
import json
import os, shutil
import subprocess
import uuid
import warnings
import itertools
import operator
import networkx
import time
from . import core
from . import unifex

import enum
if hasattr(enum, '__version__') and enum.__version__.startswith('0.4'):
    raise ImportError("""
The Python package 'enum34' is required. Unfortunately, whenever the package 'enum' is also
present it can mask 'enum34'. Make sure that both 'enum34' is installed and 'enum'
is uninstalled.
    """)
Enum = enum.Enum

#FIXME: Shouldn't the definition come from separate sources (so the coupling between creation and access
# to the persistance layer is not too much tied to this Python script ?
_TASK_TODO = 'To do'
_TASK_DONE = 'Done'
_TASK_INPROGRESS = 'Work-in-progress'
_TASK_FAILED = 'Failed'
_TASK_STATUS_LIST = {_TASK_TODO: None, 
                     _TASK_INPROGRESS: None,
                     _TASK_DONE: None,
                     _TASK_FAILED: None}

# structure to have database ID, and flag to indicate if new (just created)
DbID = namedtuple('DbID', 'id new')

class Step:
    """
    When used in the context of a StepGraph, a Step is small graph,
    or subgraph, consituted of a vertex, connected downstream to targets and upstream
    to sources. For more information about a StepGraph, see
    the documentation for it.
    """
    __slots__ = ['step', 'sources', 'targets', 'parameters', 'model']
    def __init__(self, step, 
                 sources, targets, parameters,
                 model):
        """
        :param step: A concrete step in the database
        :type step: StepConcrete_DbEntry
        :param sources: a sequence of source assets
        :param targets: a sequence of target assets
        :param model: a model (e.g., :mod:`railroadtracks.rnaseq`)
        """
        assert isinstance(step, StepConcrete_DbEntry)
        self.step = step
        self.sources = sources
        self.targets = targets
        self.parameters = parameters
        self.model = model

    clsname = property(lambda self: self.step.clsname, None, None)
    unifiedname = property(lambda self: getattr(self.model, self.clsname)._name,
                           None, None)
    def iscomplete(self):
        return self.step.status == _TASK_DONE

    def run(self):
        # FIXME: check matching version numbers
        # self.step.version

        sources = ['='.join((x.label, x.entityname)) for x in self.sources]
        targets = ['='.join((x.label, x.entityname)) for x in self.targets]

        uei = core.UnifiedExecInfo(self.step.executable, self.unifiedname,
                                   sources,
                                   targets,
                                   self.parameters,
                                   None, None # logging_file, logging_level
                               )

        res = unifex.unified_exec_run(uei, unifex._make_stepdict(self.model))
        return res



class GRAPHDISPLAY(Enum):
    STEPLEVEL = {'layout': 'dot',
                 'layoutargs': '',
                 'kwargs': {'se_nodeid': '%(clsname)s_%(label)s',
                            'se_nodelabel': '{%(clsname)s | [%(label)s]}',
                            'step_nodeid': '%(clsname)s',
                            'step_nodelabel': '%(clsname)s'}}
    TASKLEVEL = {'layout': 'dot',
                 'layoutargs': '-Goverlap=prism -Gratio=1',
                 'kwargs': {'se_nodeid': '%(id)i | %(clsname)s',
                            'se_nodelabel': '%(id)i | {%(clsname)s|[%(label)s]}',
                            'step_nodeid': '%(id)i | %(clsname)s',
                            'step_nodelabel': '%(id)i | %(clsname)s'}}

# The DAG will be build from 2 types of vertices:
# - steps
# - assets
# Edges will always link of links of 2 different types
# FIXME: if the underlying database is changing, this is not updated
# (have a lock ? monitor changes ?)
class StepGraph(object):
    """ 
    The steps to be performed are stored in a directed acyclic graph (DAG).

    This graph can be thought of as a two-level graph. The higher level represents
    the connectivity between steps (we will call supersteps), and the lower-level expands each step
    into sources, targets, and a step using sources to produce targets.
    
    There is a persistent representation (currently a mysql database),
    and this class is aiming at isolating this implementation detail from a user.
    """

    def __init__(self, persistent_graph):
        # DAG used to resolve order in which steps should be run
        dag = networkx.DiGraph()
        self._dag = dag
        self._persistent_graph = persistent_graph
        # cache of DB IDs
        self._cache_dbid = dict()
        self._cache_stepvariant_dbid = dict()
        # build graph from the persistence layer
        for step in self._persistent_graph.iter_steps():
            # step contrains StepConcrete_DbEntry, sources, targets, parameters)
            self._update_graph(step.step.step_variant_id, step.step.id, step.sources, step.targets)

    @staticmethod
    def stepconcrete_dirname(stepconcrete_id):
        """
        Name of the directory corresponding to an ID.

        :param stepconcrete_id: ID for a directory.
        """
        stepconcrete_dirname = 'step_%i' % stepconcrete_id.id
        return stepconcrete_dirname

    def add(self, step, assets, parameters=tuple(), tag=1, use_cache=True):
        """ Add a step, associated assets, and optional parameters, to the StepGraph.

        The task graph is like a directed (presumably) acyclic multilevel graph.
        Asset vertices are only connected to step vertices (in other words asset vertices
        represent connective layers between steps).

        :param step: The step to be added
        :type step: a :class:`core.StepAbstract` (or of child classes) object
        :param assets: The assets linked to the step added. If :attr:`assets.target`
                       is undefined, the method will define it with unique identifiers
                       and these will assigned in place.
        :type assets: a :class:`core.AssetStep` (or of child classes) object
        :param parameters: Parameters for the step
        :type parameters: A sequence of :class:`str` elements
        :param tag: A tag to differentiate repetitions of the exact same task.
        
        :rtype: :class:`StepConcrete_DbEntry` as the entry added to the database
        """

        # add assets (note: the _cache ensures uniqueness)
        # add step (note: the _cache ensures uniqueness)
        assert isinstance(step, StepAbstract)

        #
        if use_cache:
            task_hashdb = (step.hashdb, assets.source.hashdb, parameters, tag)
            dbid = self._cache_dbid.get(task_hashdb)
            if dbid is not None:
                # We are not done yet as we need to recover the asset.target values
                for known_target in self._persistent_graph.get_targetassets(dbid):
                    candidate_target = getattr(assets.target, known_target.label)
                    if candidate_target._defined:
                        assert candidate_target.name == known_target.name, "Mismatch between target assets for '%s': previously '%s', now '%s'" % (known_target.label, known_target.entityname, candidate_target.name)
                    else:
                        # update the definition of the target
                        candidate_target.name = known_target.entityname
                    
                return DbID(dbid, False)
            

        # obtain the id for variant
        dbid = None
        if use_cache:
            stepvariant_hashdb = (step.hashdb, tuple(x.value for x in step.activities))
            dbid = self._cache_stepvariant_dbid.get(stepvariant_hashdb)
        if dbid is None:
            id_stepvariant = self._persistent_graph.id_step_variant(step,
                                                                    step.activities).id

        # undefined sources is not accepted
        # (exceptions being the source assets that are optionally None)
        if any((not y.allownone and not x._defined) for x,y in zip(assets.source, assets.source._sources)):
            raise ValueError('All sources must be defined.')

        # retrieve or create the task
        stepconcrete_id = self._persistent_graph.id_stepconcrete(id_stepvariant,
                                                                 assets.source,
                                                                 assets.target,
                                                                 parameters,
                                                                 tag = tag)

        stepconcrete_dirname = self.stepconcrete_dirname(stepconcrete_id)
        absdirname = os.path.join(self._persistent_graph._wd, stepconcrete_dirname)
        if stepconcrete_id.new:
            if os.path.isdir(absdirname):
                raise IOError("The directory %s is already existing." % absdirname)
            else:
                os.mkdir(absdirname)
        else:
            if not os.path.isdir(absdirname):
                raise IOError("The directory %s is missing." % absdirname)

        # loop over targets
        #    # targets_db = list()
        #    # for (label, asset) in zip(targets._fields, targets):
        #    #     # each asset can represent several saved objects
        #    #     entity_ids = tuple(self.id_stored_entity(*cn).id for cn in asset.iteritems())
        #    #     targets_db.append(entity_ids)
        #    # is this concrete step already known ?
        #    # Look for the sources

        labelnamepairs = list()
        for field, t in zip(assets.target._fields, assets.target):
            if not t._defined:
                # Generate a presumably unique ID
                uniquename = str(uuid.uuid1())
                if isinstance(t, core.File) and t._extension is not None:
                    uniquename += t._extension[0]
                uniquename = os.path.join(absdirname, uniquename)
                # create an entry in the database
                id_storedentity = self._persistent_graph.id_stored_entity(type(t), 
                                                                          uniquename)
                if not id_storedentity.new:
                    # The newly created entity should be... well, NEW.
                    # if reaching here, we have a (serious) problem
                    raise Exception("Bad, bad, bad... the generated unique name %s is not unique." % uniquename)
                t.name = uniquename
                labelnamepairs.append((field, t))
            else:
                # the asset "t" is defined, and sanity checks are (should be?) in self.id_stepconcrete()
                pass
        self._persistent_graph._insert_stepconcrete2storedentities(labelnamepairs, 
                                                                   'target', 
                                                                   stepconcrete_id.id)
        self._persistent_graph.connection.commit()

        sources = list()
        for x in assets.source:
            if x is not None:
                for y in x:
                    sources.append(y)

        self._update_graph(id_stepvariant, stepconcrete_id.id, 
                           tuple(sources),
                           tuple(y for x in assets.target for y in x))
        # store in the cache
        self._cache_dbid[task_hashdb] = stepconcrete_id.id
        return stepconcrete_id

    def _update_graph(self, id_stepvariant, stepconcrete_id, sources, targets):
        # update the graph
        dag = self._dag

        v1 = '%s-%i' % (id_stepvariant, stepconcrete_id)
        if v1 not in dag:
            dag.add_node(v1)
        # add step details (sources and targets)
        steps_before = set()
        for src in sources:
            assert src is not None
            if src not in dag:
                dag.add_node(src)
            vertex = dag[src]
            # edge between the source file and the step
            dag.add_edge(src, v1)
        steps_after = set()
        for target in targets:
            assert target is not None
            if target not in dag:
                dag.add_node(target)
            vertex = dag[target]
            # edge between the step and the target file
            dag.add_edge(v1, target)

    def stepcrawler(self):
        # This is a very rudimentary way to organise tasks. 
        # Bulk synchronous parallel approachs (Google's PREGEL, Apache's Giraf and Hama, Stanford's GPS...)
        # should be considered if ever becoming a large graph.
        dag = self._dag
        # iterate through connected components
        undag = dag.to_undirected()
        for cg in networkx.connected_components(undag):
            # start from a root
            topo = networkx.topological_sort(cg)
            root = topo[0]
            # depth-first search (to get the initial products ASAP)
            for edge in dag.dfs_edges(root):
                stepnode = edge[1]
                useroot_to_do = list()
                # here the `todo` object contains a list of "immediately next"
                # steps to be performed for that one root.
                # Each step can be computed independently

                #FIXME: shouldn't the logic below be taken out of the crawler,
                #       or in a callback function ?
                for step in root.iter_steps():
                    if step.status == _TASK_DONE:
                        continue
                    elif step.status == _TASK_INPROGRESS:
                        # FIXME: have a mechanism to recover from failure, crash, etc... ?
                        continue
                    elif step.status == _TASK_TODO:
                        useroot_to_do.append(step)
                    elif step.status == _TASK_FAILED:
                        raise Exception("Step previously failed: %s" % step)
                    else:
                        # paranoid check
                        raise ValueError("Unknown step %s" % step)
                yield useroot_to_do
                # step in that list should be performed before the next iteration
        pass


    def provenancewalk_storedentity(self, stored_entity,
                                    func_stored_entity,
                                    func_stored_sequence,
                                    func_step_concrete,
                                    func_storedentity_stepconcrete,
                                    func_stepconcrete_storedentity):
        """ Walk up the path.
        :param stored_entity_id: the stored entity to start from
        :param func_stored_entity: a callback called with each stored entity
        :param func_step_concrete: a callback called with each step concrete
        :param func_storedentity_stepconcrete: a callback called with each link between a stored entity and a step concrete
        :param func_stepconcrete_storedentity: a callback called with each link between a step concrete and a stored entity
        """
        assert hasattr(stored_entity, 'id')
        storedentity_stack = collections.deque()
        storedentity_visited = set()
        step_stack = collections.deque()
        step_visited = set()
        if hasattr(stored_entity, 'iter_storedentities'):
            # sequence:
            for elt in stored_entity.iter_storedentities():
                storedentity_stack.append(elt)
                storedentity_visited.add(elt.id)
                func_stored_entity(elt)
            func_stored_sequence(stored_entity)
        else:
            # not a sequence
            storedentity_stack.append(stored_entity)
            storedentity_visited.add(stored_entity.id)
            func_stored_entity(stored_entity)

        while len(storedentity_stack) > 0:
            stored_entity = storedentity_stack.popleft()
            if hasattr(stored_entity, 'iter_storedentities'):
                # this is a sequence
                for elt in stored_entity.iter_storedentities():
                    storedentity_stack.append(elt)
                    storedentity_visited.add(elt.id)
                    func_stored_entity(elt)
                func_stored_sequence(stored_entity)
                stored_entity = storedentity_stack.popleft()                
            step_concrete = self._persistent_graph.get_parenttask_of_storedentity(stored_entity)
            if step_concrete is not None:
                if step_concrete.id not in step_visited:
                    step_stack.append(step_concrete)
                    func_step_concrete(step_concrete)
                    step_visited.add(step_concrete.id)
            while len(step_stack) > 0:
                step_concrete = step_stack.popleft()
                func_stepconcrete_storedentity(step_concrete, stored_entity)
                storedentities = self._persistent_graph.get_srcassets(step_concrete.id)
                for entity in storedentities:
                    if entity.id not in storedentity_visited:
                        storedentity_stack.append(entity)
                        func_stored_entity(entity)
                        storedentity_visited.add(entity.id)
                        #func_stored_entity(entity)
                    func_storedentity_stepconcrete(entity, step_concrete)


    def _destination_walk(self, storedentity_stack, storedentity_visited,
                          step_stack, step_visited,
                          func_stored_entity,
                          func_step_concrete,
                          func_storedentity_stepconcrete,
                          func_stepconcrete_storedentity):
        while len(storedentity_stack) > 0:
            stored_entity = storedentity_stack.popleft()
            steps_concrete = self._persistent_graph.get_targetstepconcrete(stored_entity)
            for entity in steps_concrete:
                if entity.id not in step_visited:
                    step_stack.append(entity)
                    step_visited.add(entity.id)
            func_stored_entity(stored_entity)
            while len(step_stack) > 0:
                step_concrete = step_stack.popleft()
                func_step_concrete(step_concrete)
                func_stepconcrete_storedentity(step_concrete, stored_entity)
                storedentities = self._persistent_graph.get_targetassets(step_concrete.id)
                for entity in storedentities:
                    if entity.id not in storedentity_visited:
                        storedentity_stack.append(entity)
                        storedentity_visited.add(entity.id)
                    func_storedentity_stepconcrete(entity, step_concrete)
        return

    def destinationwalk_stepconcrete(self, step_concrete_id,
                                     func_stored_entity,
                                     func_stored_sequence,
                                     func_step_concrete,
                                     func_storedentity_stepconcrete,
                                     func_stepconcrete_storedentity):
        """ Walk down the path."""
        storedentity_stack = collections.deque()
        storedentity_visited = set()
        step_stack = collections.deque()
        step_visited = set()

        #FIXME: clean the code by taking the nested 'while' loop inside-out
        step_concrete = self._persistent_graph._get_stepconcrete(DbID(step_concrete_id, 
                                                                      False))
        func_step_concrete(step_concrete)
        step_visited.add(step_concrete_id)
        stored_entities = self._persistent_graph.get_targetassets(step_concrete_id)
        for entity in stored_entities:
            if entity.id not in storedentity_visited:
                storedentity_stack.append(entity)
                storedentity_visited.add(entity.id)
            else:
                raise Exception('Step %s has several identical targets.' % str(step_concrete_id))
            #func_stepconcrete_storedentity(step_concrete, entity)
            func_storedentity_stepconcrete(entity, step_concrete)
        return self._destination_walk(storedentity_stack, storedentity_visited,
                                      step_stack, step_visited,
                                      func_stored_entity,
                                      func_step_concrete,
                                      func_storedentity_stepconcrete,
                                      func_stepconcrete_storedentity)

    def destinationwalk_storedentity(self, 
                                     stored_entity,
                                     func_stored_entity,
                                     func_stored_sequence,
                                     func_step_concrete,
                                     func_storedentity_stepconcrete,
                                     func_stepconcrete_storedentity):
        """ Walk down the path."""
        storedentity_stack = collections.deque()
        storedentity_visited = set()
        step_stack = collections.deque()
        step_visited = set()
        
        storedentity_stack.append(stored_entity)
        storedentity_visited.add(stored_entity.id)
        return self._destination_walk(storedentity_stack, storedentity_visited,
                                      step_stack, step_visited,
                                      func_stored_entity,
                                      func_step_concrete,
                                      func_storedentity_stepconcrete,
                                      func_stepconcrete_storedentity)

    def _graph_storedentity(self, db_id,
                            func, opposite,
                            display = GRAPHDISPLAY.STEPLEVEL):
        """Make a provenance graph."""
        se_nodeid = display.value['kwargs']['se_nodeid']
        se_nodelabel = display.value['kwargs']['se_nodelabel']
        step_nodeid = display.value['kwargs']['step_nodeid']
        step_nodelabel = display.value['kwargs']['step_nodelabel']
        def gst(stored_entity, nodelabel):
            return nodelabel % dict((x, getattr(stored_entity, x)) for x in ('id', 'clsname', 'label'))
        def gsc(task, nodelabel):
            return nodelabel % dict((x, getattr(task, x)) for x in ('id', 'clsname'))
        dag = networkx.DiGraph()
        def func_stored_entity(stored_entity):
            dag.add_node(gst(stored_entity, se_nodeid),
                         label=gst(stored_entity, se_nodelabel),
                         shape='record')
        def func_stored_sequence(stored_sequence):
            #dag.add_node(gst(stored_sequence, se_nodelabel),
            #             label=gst(stored_sequence, se_nodelabel),
            #             shape='record')
            for se in stored_sequence.iter_storedentities():
                dag.add_edge(gst(se, se_nodeid),
                             gst(stored_sequence, se_nodeid))                
        def func_step_concrete(step_concrete):
            if step_concrete.status == _TASK_STATUS_LIST[_TASK_DONE]:
                style = 'filled'
            else:
                style = ''
            dag.add_node(gsc(step_concrete, step_nodeid),
                         label=gsc(step_concrete, step_nodelabel),
                         shape='box', 
                         style = style,
                         xlabel=str(step_concrete.status))
        def func_storedentity_stepconcrete(stored_entity, step_concrete):
            if opposite:
                dag.add_edge(gsc(step_concrete, step_nodeid),
                             gst(stored_entity, se_nodeid))                             
            else:
                dag.add_edge(gst(stored_entity, se_nodeid), 
                             gsc(step_concrete, step_nodeid))
        def func_stepconcrete_storedentity(step_concrete, stored_entity):
            if opposite:
                dag.add_edge(gst(stored_entity, se_nodeid),
                             gsc(step_concrete, step_nodeid))                             
            else:
                dag.add_edge(gsc(step_concrete, step_nodeid), 
                             gst(stored_entity, se_nodeid))

        func(db_id,
             func_stored_entity,
             func_stored_sequence,
             func_step_concrete,
             func_storedentity_stepconcrete,
             func_stepconcrete_storedentity)
        return dag

    def provenancegraph_storedentity(self, stored_entity, 
                                     display = GRAPHDISPLAY.STEPLEVEL):
        return self._graph_storedentity(stored_entity,
                                        self.provenancewalk_storedentity,
                                        False, # opposite
                                        display = display)
    def destinationgraph_stepconcrete(self, step_concrete,
                                      display = GRAPHDISPLAY.STEPLEVEL):
        return self._graph_storedentity(step_concrete,
                                        self.destinationwalk_stepconcrete,
                                        True, # opposite
                                        display = display)
    def destinationgraph_storedentity(self, stored_entity,
                                      display = GRAPHDISPLAY.STEPLEVEL):
        return self._graph_storedentity(stored_entity,
                                        self.destinationwalk_storedentity,
                                        True, # opposite
                                        display = display)




    def cleantargets_stepconcrete(self, step_concrete_id):
        """ Clean the targets downstream of a task (step_concrete),
        which means erasing the target files and (re)setting the tasks' status
        to 'TO DO'.
        
        :param step_concrete_id: A task
        :type step_concrete_if: a :class:`DbID` (or anything with an attribute :attr:`id`).
        """
        def func_stored_entity(stored_entity):
            Cls = getattr(self._persistent_graph._model, 
                          stored_entity.clsname)
            instance = Cls(stored_entity.entityname)
            if hasattr(instance, 'iterlistfiles'):
                nameiter = (os.path.join(os.path.dirname(instance.name), basename) \
                            for elt_cls, basename in instance.iterlistfiles())
            else:
                nameiter = iter(instance)
            for pathname in nameiter:
                if os.path.isfile(pathname):
                    os.remove(pathname)
                elif os.path.isdir(pathname):
                    shutil.rmtree(pathname)
                else:
                    warnings.warn("The pathname '%s' associated with the store entity '%s' cannot be removed." % (pathname, str(stored_entity)))
        def func_step_concrete(step_concrete):
            self._persistent_graph.step_concrete_state(step_concrete,
                                                       _TASK_STATUS_LIST[_TASK_TODO])
        def func_storedentity_stepconcrete(stored_entity, step_concrete):
            pass
        def func_stepconcrete_storedentity(step_concrete, stored_entity):
            pass
        def func_stored_sequence(step_concrete, stored_entity):
            pass

        self.destinationwalk_stepconcrete(
            step_concrete_id.id,
            func_stored_entity,
            func_stored_sequence,
            func_step_concrete,
            func_storedentity_stepconcrete,
            func_stepconcrete_storedentity)
    

import sqlite3
import os

sql_fn = os.path.join(os.path.dirname(__file__), 'cache.sql')

class StepConcrete_DbEntry(object):
    __slots__ = ['id', 'status', 'step_variant_id', 'steptype_id', 'executable', 'clsname', 'version', 'parameters']
    def __init__(self, id, status, step_variant_id, 
                 steptype_id, executable, clsname, version, 
                 parameters=()):
        self.id = id
        self.status = status
        self.step_variant_id = step_variant_id
        self.steptype_id = steptype_id
        self.executable = executable
        self.clsname = clsname
        self.version = version
        self.parameters = parameters

TaskStatusCount = namedtuple('TaskStatus', 'id label count')

class StoredEntity(object):
    __slots__ = ('id', 'label', 'clsname', 'entityname')
    def __init__(self, id, label, clsname, entityname):
        self.id = id
        self.label = label
        self.clsname = clsname
        self.entityname = entityname
    def resurrect(self, model):
        if self.clsname == 'NoneType':
            return None
        elif '.' in self.clsname:
            module_name, cls_name = self.clsname.rsplit('.', 1)
            module = import_module(module_name)
        else:
            module = model
            cls_name = self.clsname
        Cls = getattr(module, cls_name)
        return Cls(self.entityname)
    def __str__(self):
        return '\n'.join((super(StoredEntity, self).__str__(),
                          'id: %i' % self.id))

class StoredSequence(object):
    __slots__ = ('id', 'label', 'clsname', 'storedentities')
    def __init__(self, id, label, clsname, storedentities):
        self.id = id
        self.label = label
        self.clsname = clsname
        self.storedentities = storedentities
    def resurrect(self, model):
        if '.' in self.clsname:
            module_name, cls_name = self.clsname.rsplit('.', 1)
            module = import_module(module_name)
        else:
            module = model
            cls_name = self.clsname
        Cls = getattr(module, cls_name)
        assert issubclass(Cls, core.FileSequence), 'Expected a FileSequence but got %s' % str(Cls)
        return Cls(x.resurrect(model) for x in self.storedentities)
    def iter_storedentities(self):
        return iter(self.storedentities)
    def __iter__(self):
        return self.iter_storedentities()
    def __len__(self):
        return len(self.storedentities)
    def __str__(self):
        return '\n'.join((super(StoredSequence, self).__str__(),
                          'id: %i' % self.id,
                          'length: %i' % len(self)))


class PersistentTaskGraph(object):
    """
    List of tasks stored on disk.
    
    """

    StoredEntityNoLabel = namedtuple('StoredEntityNoLabel', 'id clsname entityname')

    def __init__(self, db_fn, model, wd='.', force_create=False, isolation_level=None):
        """
        :param db_fn: file name for the database.
        :param wd: working directory (where derived data files are stored)
        :param force_create: recreate database if a file with the same name is already
                             present.
        :param isolation_level: passed to :func:`sqlite3.connection`.
        """

        # check whether the file containing the SQLite database exists.
        # 
        if os.path.exists(db_fn):
            if os.path.isfile(db_fn):
                if force_create:
                    create = True
                else:
                    create = False
            else:
                raise ValueError('db_fn should be a file name, not a directory name.')
        else:
            create = True
        self._db_fn = db_fn
        self._model = model
        self._model_steplist = unifex._make_stepdict(model)
        self._wd = wd
        self.created = create
        connection = sqlite3.connect(db_fn, isolation_level=isolation_level)
        self.connection = connection
        if create:
            self._create(db_fn, wd)
        self._fill_tasklist()
        self._statuslist = None # lazy evaluation

    def _create(self, db_fn, wd):
        with open(sql_fn, 'r') as fh:
            sql = fh.readlines()
        sql = ''.join(sql)
        connection = self.connection
        with connection:
            for block in sql.split(';'):
                connection.execute(block)
        # misc
        #    version
        sql = """
        INSERT INTO misc (tag, description, val) VALUES (?, ?, ?)
        """
        with connection:
            res = connection.executemany(sql, 
                                         (('version','Version for the framework',
                                           __version__),
                                          ('workingdirectory', 'Directory in which files are saved',
                                           wd)))

        # populate possible status for tasks
        sql = """
        INSERT INTO step_status (label) VALUES (?)
        """
        with connection:
            res = connection.executemany(sql, ((x,) for x in _TASK_STATUS_LIST.keys()))

    def _fill_tasklist(self):
        cursor = self.connection.cursor()
        sql = """
        SELECT id, label
        FROM step_status
        """
        res = cursor.execute(sql)
        for x in res:
            _TASK_STATUS_LIST[x[1]] = x[0]


    def getversion(self):
        sql = """
        SELECT val 
        FROM misc
        WHERE tag=='version'
        """
        cursor = self.connection.cursor()
        res = cursor.execute(sql)
        return res.fetchone()[0]
    version = property(getversion, None, None,
                       "Version for the database and package (mixing versions comes at one's own risks)")

    def iter_steps(self):
        """ Iterate through the concrete steps """
        sql = """
        SELECT step_concrete.id AS id,
               step_status.label AS status,
               sv.id,
               sv.step_type_id,
               sv.executable,
               sv.cls,
               sv.version
        FROM step_concrete, step_status, step_variant AS sv
        WHERE step_concrete.step_status_id=step_status.id
        AND step_concrete.step_variant_id=sv.id
        ORDER BY step_concrete.id -- ID is autoincremented, so this is like a chronological order
        """
        # If the SQL above changes, the definition of
        # StepConcrete_DbEntry might have to be updated
        cursor = self.connection.cursor()
        steps = cursor.execute(sql)
        # loop through the steps
        for row in steps:
            # for each step, get the sources
            sc = StepConcrete_DbEntry(*row)
            sources = tuple(self.get_srcassets(sc.id))
            targets = tuple(self.get_targetassets(sc.id))
            parameters = tuple(self.get_parameters(sc.id))
            yield Step(sc, sources, targets, parameters, self._model)

    def _get_assets(self, concrete_step, kind):
        if hasattr(concrete_step, 'id'):
            concrete_step_id = concrete_step.id
        else:
            concrete_step_id = int(concrete_step)
        cursor = self.connection.cursor()
        # First yield stored entities
        sql_template = """
        SELECT stored_entity.id as id,
               step_concrete2%(kind)sfile.label as label,
               stored_entity.classname as classname,
               stored_entity.entityname as entityname
        FROM step_concrete2%(kind)sfile, stored_entity
        WHERE step_concrete_id=?
        AND stored_entity_id=stored_entity.id
        """
        sql = sql_template % {'kind': kind}
        res = cursor.execute(sql, (concrete_step_id,))
        for src in res:
            res = StoredEntity(*src)
            yield res
        # Now yield sequences (of stored entities)
        sql_template_ss = """
        SELECT stored_sequence.id,
               sc.label as label,
               stored_sequence.classname as classname
        FROM step_concrete2%(kind)sfile AS sc
        JOIN stored_sequence
            ON stored_sequence.id=sc.stored_sequence_id
        WHERE step_concrete_id=?
        """
        sql = sql_template_ss % {'kind': kind}
        sql_se = """
        SELECT stored_entity.id as id,
               stored_entity.classname as classname,
               stored_entity.entityname as entityname,
               pos
        FROM stored_entity2sequence AS se2s
        JOIN stored_entity
            ON se2s.stored_entity_id=stored_entity.id
        WHERE se2s.stored_sequence_id=?
        ORDER BY id, classname, pos
        """
        res = cursor.execute(sql, (concrete_step_id,))
        cursor2 = self.connection.cursor()
        for src in res:
            res2 = cursor2.execute(sql_se, (src[0],))
            ses = list()
            for pos, x in enumerate(res2.fetchall()):
                if x[-1] != pos:
                    raise Exception('The sequence is missing an item in postion %i (and got position %i instead).' % (pos, x[-1]))            
                ses.append(StoredEntity(x[0], None, x[1], x[2]))
            res = StoredSequence(src[0], src[1], src[2], ses)
            yield res


    def task_time_points(self, task_id):
        cursor = self.connection.cursor()
        sql = """
        SELECT time_creation,
               time_t0,
               time_t1
        FROM step_concrete
        WHERE step_concrete.id=?
        """
        cursor.execute(sql, (task_id.id, ))
        res = cursor.fetchone()
        return res

    def set_time_t0(self, task_id, t0):
        cursor = self.connection.cursor()

        sql = """
        UPDATE step_concrete
        SET time_t0=?
        WHERE step_concrete.id=?
        """
        cursor.execute(sql, (t0, task_id.id, ))
        self.connection.commit()

    def set_time_t1(self, task_id, t1):
        cursor = self.connection.cursor()

        sql = """
        UPDATE step_concrete
        SET time_t1=?
        WHERE step_concrete.id=?
        """
        cursor.execute(sql, (t1, task_id.id, ))
        self.connection.commit()
        
    def _get_stepconcrete(self, step_concrete_id):
        cursor = self.connection.cursor()
        sql = """
        SELECT step_concrete.id,
               step_concrete.step_status_id,
               step_concrete.step_variant_id,
               step_type_id,
               executable,
               cls,
               version,
               step_parameters.json
        FROM step_concrete
        INNER JOIN step_status
        ON step_concrete.step_status_id=step_status.id
        INNER JOIN step_variant AS sv
        ON step_concrete.step_variant_id=sv.id
        INNER JOIN step_concrete2parameters AS sc2p
        ON sc2p.step_concrete_id=step_concrete.id
        INNER JOIN step_parameters
        ON step_parameters_id=step_parameters.id
        WHERE step_concrete.id=?
        """
        cursor.execute(sql, (step_concrete_id.id, ))
        row = cursor.fetchone()
        if row is None:
            # Trouble. Trying to provide an informative error message
            sql = """
            SELECT *
            FROM step_concrete
            WHERE id=?
            """
            cursor.execute(sql)
            if len(cursor.fetchall()) == 0:
                raise ValueError("The task ID %i cannot be found." % step_concrete_id.id)
            else:
                raise ValueError("Error with the database. Likely because of interrupted process while entering tasks.")

        res = StepConcrete_DbEntry(*row)
        return res



    _SQL_TASK_FROMASSET_TEMPLATE = """
        SELECT step_concrete.id,
               step_concrete.step_status_id,
               step_concrete.step_variant_id,
               step_type_id,
               executable,
               cls,
               version,
               step_parameters.json
        FROM step_concrete2%(src_or_target)sfile AS sc2f
        INNER JOIN step_concrete
        ON sc2f.step_concrete_id=step_concrete.id
        INNER JOIN step_variant
        ON step_variant_id=step_variant.id
        INNER JOIN step_concrete2parameters AS sc2p
        ON sc2p.step_concrete_id=step_concrete.id
        INNER JOIN step_parameters
        ON step_parameters_id=step_parameters.id
        WHERE sc2f.stored_%(entity_or_sequence)s_id=?
    """
    _SQL_TASK_TOTARGETSEQUENCE = _SQL_TASK_FROMASSET_TEMPLATE  % {'src_or_target': 'target',
                                                                  'entity_or_sequence': 'sequence'}
    _SQL_TASK_FROMSRCSEQUENCE  = _SQL_TASK_FROMASSET_TEMPLATE % {'src_or_target': 'src',
                                                                 'entity_or_sequence': 'sequence'}
    _SQL_TASK_TOTARGETSENTITY  = _SQL_TASK_FROMASSET_TEMPLATE % {'src_or_target': 'target',
                                                                 'entity_or_sequence': 'entity'}
    _SQL_TASK_FROMSRCSENTITY   = _SQL_TASK_FROMASSET_TEMPLATE % {'src_or_target': 'src',
                                                                 'entity_or_sequence': 'entity'}

    def _get_stepconcrete_from_storedentity(self, stored_entity_id, sql):
        """ 
        Retrieve the step_concrete object(s) linked with a stored_entity_id.
        Whenever :param:`kind` is equal to "target" there will be at most
        one stored_entity (if there is more one, the database has a consistentcy
        issue).
        :param stored_entity_id: the ID for a stored_entity
        :param kind: either 'src' or 'target'
        """
        cursor = self.connection.cursor()
        res = cursor.execute(sql, (stored_entity_id,))
        for src in res:
            res = StepConcrete_DbEntry(*src)
            yield res

    def get_parenttask_of_storedentity(self, stored_entity):
        """ Return the task producing a stored entity. There should obviously only be one such task,
        and an Exception is raised if not the case.
        :param stored_entity: the stored entity in the database
        :type stored_entity: :class:`StoredEntity`
        :rtype: a `StepConcrete_DbEntry` :class:`namedtuple`, or None
        """
        assert hasattr(stored_entity, 'id')
        step_concrete = None
        if hasattr(stored_entity, 'iter_storedentities'):
            raise ValueError('Only atomic stored entities are accepted.')
        else:
            sql = self._SQL_TASK_TOTARGETSENTITY
        for i, step_concrete in enumerate(self._get_stepconcrete_from_storedentity(stored_entity.id, sql)):
            if i > 0:
                raise Exception("""
 Consistency issue with the database. More than one step is claiming to be the source of the stored_entity_id "%s" """ % str(stored_entity.id))
        return step_concrete

    def get_targetstepconcrete(self, stored_entity):
        """ Return the tasks using a given stored entity.
        :param stored_entity: the stored entity in the database.
        :type stored_entity: can be :class:`StoredEntity` or :class:`StoredSequence`
        :rtype: a `SepConcrete_DbEntry` :class:`namedtuple`, or None
        """
        if hasattr(stored_entity, 'iter_storedentities'):
            raise NotImplementedError('Handling sequence as source assets is not yet implemented.')
        else:
            sql = self._SQL_TASK_FROMSRCSENTITY

        return tuple(self._get_stepconcrete_from_storedentity(stored_entity.id, 
                                                              sql))
        
    def get_srcassets(self, concrete_step_id):
        """ Return the source files for a given concrete step ID.
        :param concrete_step_id: ID for the concrete step in the database.
        :rtype: generator
        """
        return self._get_assets(concrete_step_id, 'src')

    def get_targetassets(self, concrete_step_id):
        """ Return the target files for a given concrete step ID.
        :param concrete_step_id: ID for the concrete step in the database.
        :rtype: generator
        """
        return self._get_assets(concrete_step_id, 'target')

    def get_parameters(self, concrete_step_id):
        cursor = self.connection.cursor()
        sql = """
        SELECT step_parameters.id as id,
               step_parameters.json as json
        FROM step_concrete2parameters, step_parameters
        WHERE step_concrete_id=?
        AND step_parameters_id=step_parameters.id
        """
        res = cursor.execute(sql, (concrete_step_id,))
        return res.fetchall()

    def get_statuslist(self):
        if self._statuslist is None:
            cursor = self.connection.cursor()
            sql = """
            SELECT * from step_status
            """
            res = cursor.execute(sql)
            self._statuslist = tuple(res)
        return self._statuslist

    statuslist = property(get_statuslist, None, None,
                          """ Status list """)

    def id_stepparameters(self, parameters):
        """
        Conditionally add parameters (add only if not already present)
        :param parameters: sequence of parameters
        :rtype: ID for the pattern as a :class:`DbID`.
        """
        cursor = self.connection.cursor()
        param_json = json.dumps(parameters)
        sql = """
        SELECT id 
        FROM step_parameters
        WHERE json=?
        """
        res = cursor.execute(sql, (param_json,))
        res = cursor.fetchone()
        if res is None:
            sql = """
            INSERT INTO step_parameters (
            json
            ) VALUES (
            ?);
            """
            cursor.execute(sql, (param_json,))
            db_id = cursor.lastrowid
            self.connection.commit()
            res = DbID(db_id, True)
        else:
            res = DbID(res[0], False)
        return res

        
    def id_stored_entity(self, cls, name):
        """
        Conditionally add a stored entity (add only if not already present)
        :param cls: Python class for the stored entity
        :param name: Parameter "name" for the class "cls".
        :rtype: ID for the pattern as a :class:`DbID`.
        """
        cursor = self.connection.cursor()
        sql = """
        SELECT id 
        FROM stored_entity
        WHERE classname=?
        AND entityname=?;
        """
        res = cursor.execute(sql, (cls.__name__, name))
        res = cursor.fetchone()
        if res is None:
            sql = """
            INSERT INTO stored_entity (
            classname,
            entityname
            ) VALUES (
            ?,
            ?);
            """
            try:
                cursor.execute(sql, (cls.__name__, name))
            except sqlite3.IntegrityError as ie:
                # print an informative message before propagating the exception
                sql = """
                SELECT classname
                FROM stored_entity
                WHERE entityname=?
                """
                res = cursor.execute(sql, (name,))
                res = cursor.fetchone()
                msg = 'ERROR: The descriptor "%s" for the type "%s" is already in use with type "%s"' % (name, cls.__name__, res[0])
                print(msg)
                raise(ie)
            db_id = cursor.lastrowid
            self.connection.commit()
            res = DbID(db_id, True)
        else:
            res = DbID(res[0], False)
        return res

    def id_stored_sequence(self, cls, clsname_sequence):
        """
        Conditionally add a sequence of stored entities (add only if not already present)
        :param clsname_sequence: Sequence of pairs (Python class for the stored entity, parameter "name" for the class "cls")
        :rtype: ID for the sequence as a :class:`DbID`.
        """
        cursor = self.connection.cursor()

        sql = """
        SELECT stored_sequence_id
        FROM stored_entity2sequence
        WHERE stored_entity_id=?
        AND pos=?;
        """
        stored_entity_ids = list()
        candidates = None
        # iterate through the entities in the sequence
        for pos, (cls_elt, name_elt) in enumerate(clsname_sequence):
            # ensure that the entity is tracked
            se_id = self.id_stored_entity(cls_elt, name_elt)
            stored_entity_ids.append(se_id)
            res = cursor.execute(sql, (se_id.id, pos))
            if candidates is None:
                candidates = set(x[0] for x in cursor.fetchall())
            else:
                candidates = candidates.intersection(set(x[0] for x in cursor.fetchall()))
        if candidates is None or len(candidates) == 0:
            sql = """
            INSERT INTO stored_sequence (
            classname
            ) VALUES (
            ?);
            """
            cursor.execute(sql, (cls.__name__, ))
            db_id = cursor.lastrowid
            # 
            sql = """
            INSERT INTO stored_entity2sequence (
            stored_sequence_id,
            stored_entity_id,
            pos
            ) VALUES (
            ?, ?, ?);
            """
            for pos, (se_id) in enumerate(stored_entity_ids):
                res = cursor.execute(sql, (db_id, se_id.id, pos))
            
            self.connection.commit()
            res = DbID(db_id, True)
        elif len(candidates) == 1:
            res = DbID(next(iter(candidates)), False)
        else:
            raise Exception("Consistency issue with the database.")
        return res

    def id_step_activity(self, activity):
        """ 
        Conditionally add an activity (add only if not already present)
        :param activity: one actibity name
        :rtype: ID for the activity as an integer
        """
        cursor = self.connection.cursor()
        sql = """
        SELECT id
        FROM step_activity
        WHERE label=?;
        """
        res = cursor.execute(sql, (activity.value,))
        res = cursor.fetchone()
        if res is None:
            sql = """
            INSERT INTO step_activity
            (label) VALUES (?);
            """
            cursor.execute(sql, (activity.value,))
            res = cursor.lastrowid
            self.connection.commit()
            res = DbID(res, True)
        else:
            res = DbID(res[0], False)
        return res

    def id_step_type(self, activities):
        """
        Conditionally add a step type (add only if not already present).
        :param activities: sequence of activity names
        :rtype: ID for the step type as an integer
        """

        assert len(activities) > 0, "The number of activities must be > 0"

        cursor = self.connection.cursor()

        # first get the activity IDs
        activity_ids = tuple(self.id_step_activity(x) for x in activities)
        steptype_id = None
        if not any(x.new for x in activity_ids):
            # query if any step_type already has all the associated activities
            # The if statement saves that trouble if any activity had to be created
            # (since then the step_type cannot possibly be already in the database)
            sql_activity = """SELECT step_type_id FROM step_type2activity WHERE step_activity_id = %i"""

            sql = """
            SELECT step_type_id
            FROM
              (SELECT step_type_id
               FROM step_type2activity
               WHERE step_activity_id IN
            """ + \
                '    (' + '\nINTERSECT\n'.join((sql_activity % dbid.id) for dbid in activity_ids) + ')' + \
                """
                EXCEPT
                SELECT step_type_id
                FROM step_type2activity
                WHERE step_activity_id NOT IN (%s)
                )
                """ % ','.join(str(dbid.id) for dbid in activity_ids)
            res = cursor.execute(sql)

            steptype_id = res.fetchall()
            
        if steptype_id is None or len(steptype_id) == 0:
            # if no activity ID, we should create one
            sql = """
            INSERT INTO step_type VALUES(null);
            """
            res = cursor.execute(sql)
            steptype_id = cursor.lastrowid

            sql = """
            INSERT INTO step_type2activity
            (step_type_id, step_activity_id)
            VALUES (?, ?);
            """
        
            for activity_id in activity_ids: 
                cursor.execute(sql, (steptype_id, activity_id.id))
            self.connection.commit()
            return DbID(steptype_id, True)
        elif len(steptype_id) > 1:
            #FIXME: if more than one steptype_id found, we have a problem
            raise Exception("Houston, we have problem. We have several step types (%s) with all the activities." % repr(activities))

        steptype_id = steptype_id[0][0]
        return DbID(steptype_id, False)

    def _iter_step_type(self):
        sql = """
        SELECT step_type.id, activity
        FROM step_type
        INNER JOIN step_type2activity
        ON step_type2activity.step_type_id=step_type.id
        INNER JOIN step_activity
        ON step_type2activity.step_activity_id=step_activity.id
        GROUP BY step_type.id
        """
        cursor = self.connection.cursor()
        cursor.execute(sql)
        while True:
            rows = cursor.fetchmany(100)
            if not rows:
                break
            for row in rows:
                yield row
            

    def id_step_variant(self,
                        step,
                        activities, # list of activities the step is covering (step_id inferred from that)
                    ):
        """
        Return a database ID for the step variant (creating a new ID only of the variant is not already tracked)
        
        :param step: a step
        :type step: :class:`core.StepAbstract`
        :param activities: a sequence of activity names

        :rtype: ID for a step variant as an :class:`int`.
        """

        assert isinstance(step, StepAbstract)

        executable = step._execpath
        version = step.version

        step_type_id = self.id_step_type(activities)

        cursor = self.connection.cursor()
        step_variant_id = None
        if not step_type_id.new:
            # not a new step_id, so may be the step_variant is known as well
            if executable is None:
                sql = """
                SELECT id
                FROM step_variant
                WHERE step_type_id=?
                AND executable IS NULL
                AND cls=?
                AND version=?
                """
                t = type(step)
                sql_params = (step_type_id.id,
                              '.'.join((t.__module__, t.__name__)),
                              version)
            else:
                sql = """
                SELECT id
                FROM step_variant
                WHERE step_type_id=?
                AND executable=?
                AND cls=?
                AND version=?
                """
                t = type(step)
                sql_params = (step_type_id.id,
                              executable,
                              '.'.join((t.__module__, t.__name__)),
                              version)
            res = cursor.execute(sql, sql_params)
            step_variant_id = res.fetchone()
            # FIXME: test if more than one ? 
            if step_variant_id is None:
                step_variant_id = None
            else:
                step_variant_id = step_variant_id[0] # ID only
        if step_variant_id is None:
            sql = """
            INSERT INTO step_variant
            (step_type_id, executable, cls, version)
            VALUES
            (?, ?, ?, ?)
            """
            t = type(step)
            res= cursor.execute(sql, (step_type_id.id,
                                      executable,
                                      '.'.join((t.__module__, t.__name__)),
                                      version))
            step_variant_id = cursor.lastrowid
            self.connection.commit()
            return DbID(step_variant_id, True)
        else:
            return DbID(step_variant_id, False)


    @property
    def nconcrete_steps(self):
        cursor = self.connection.cursor()
        sql = """
        SELECT COUNT(id)
        FROM step_concrete
        """
        res = cursor.execute(sql)
        return res.fetchone()[0]

    @property
    def nconcrete_steps_status(self):
        cursor = self.connection.cursor()
        sql = """
        SELECT step_status.id, step_status.label, COUNT(*)
        FROM step_concrete,
             step_status
        WHERE step_concrete.step_status_id=step_status.id
        GROUP BY step_status.id

        """
        res = cursor.execute(sql)
        return tuple(TaskStatusCount(*x) for x in res.fetchall())

    def id_stepconcrete(self, step_variant_id,
                        sources, targets, parameters,
                        tag = 1):
        """
        Conditionally add a task ("concrete" step),
        that is a step variant (executable and parameters)
        to which source and target files, as well as parameters, are added.

        :param step_variant_id: ID for the step variant
        :type step_variant_id: integer
        :param sources: sequence of sources
        :type sources: :class:`AssetSet`
        :param targets: sequence of targets
        :type targets: :class:`AssetSet`
        :param parameters: list of parameters
        :type parameters: a sequence of :class:`str`
        :param tag: a tag, used to performed repetitions of the exact same task
        :type tag: a sequence of :class:`int`

        :rtype: :class:`DbID`

        """

        assert isinstance(step_variant_id, int)

        cursor = self.connection.cursor()

        # Initialize at step_concrete_id to None. If still None later,
        # it will mean that there is no matching task already tracked.
        step_concrete_id = None 
        #FIXME: transaction (in case one of the inserts fails)
        # DB IDs for all sources

        # ID for the parameters
        parameters_id = self.id_stepparameters(parameters)

        if len(sources) == 0:
            # special case: no sources
            # Query whether the task is already known
            sql = """
            SELECT id,
                   step_variant_id,
                   step_status_id
            FROM step_concrete
            WHERE step_variant_id=?
            AND tag=?
            AND id NOT IN (SELECT id FROM step_concrete2srcfile)
            GROUP BY step_concrete.id
            """
            cursor.execute(sql, (step_variant_id, tag))
            res = cursor.fetchall()
            if len(res) == 1:
                step_concrete_id = res[0][0]
            elif len(res) == 0:
                # The task is not known (step_concrete_id is already set to None)
                pass
            elif len(res) > 1:
                # Several tasks are matching. This is not good.
                raise Exception("Serious trouble.")
            else:
                raise Exception("We should never be here.")
        else:
            # Query whether the task is already known. Information about the provenance
            # (including parameters) is enough.
            sources_db = list()
            sql_template = """
            SELECT step_concrete.id,
                   step_variant_id,
                   step_status_id,
                   stored_entity_id,
                   stored_sequence_id
            FROM step_concrete,
                 step_concrete2srcfile,
                 step_concrete2parameters
            WHERE step_variant_id=?
            AND step_parameters_id=?
            AND label=?
            AND step_concrete2parameters.step_concrete_id=step_concrete.id
            AND step_concrete2srcfile.step_concrete_id=step_concrete.id
            AND stored_entity_id%(se)s
            AND stored_sequence_id%(ss)s
            AND step_concrete.tag=?
            """
            sql_se = sql_template % {'se': '=?', 'ss': ' IS NULL'}
            sql_ss = sql_template % {'ss': '=?', 'se': ' IS NULL'}
            candidates = None # candidate known tasks
            for (label, asset, assetattr) in zip(sources._fields, sources, sources._sources):
                # loop over the source assets
                if asset is None and assetattr.allownone:
                    continue
                asset_items = tuple(cn for cn in asset.iteritems())
                if isinstance(asset, core.FileSequence):
                    # if a sequence we also need to ensure that the sequence itself is tracked
                    sequence_id = self.id_stored_sequence(type(asset), asset.iteritems()).id
                    se_id = None
                    etid = (step_variant_id, 
                            parameters_id.id,
                            label,
                            #se_id,
                            sequence_id,
                            tag)
                    sql = sql_ss
                else:
                    # if not a sequence we need to ensure that object (file) is tracked
                    sequence_id = None
                    for se_i, cn in enumerate(asset.iteritems()):
                        if se_i > 0:
                            # we should never reach here, because it was tested above
                            raise Exception("Only FileSequence objects should have several saved entities.")
                        # ensure the that the stored ID for this asset is tracked in the DB
                        se_id = self.id_stored_entity(*cn).id
                    # create a tuple of parameters for the query that selects stored entities.
                    etid = (step_variant_id, 
                            parameters_id.id,
                            label,
                            se_id,
                            #sequence_id,
                            tag)
                    sql = sql_se
                # query with parameters
                cursor.execute(sql, etid)
                # retrieve the task_id using this source asset
                tmp = set(x[0] for x in cursor.fetchall())
                if candidates is None:
                    # first iteration, set is whatever is retrieved
                    candidates = tmp
                else:
                    # otherwise, the candidate task_ids are the intersection of candidates so far and what was returned 
                    candidates = candidates.intersection(tmp)
                if len(candidates) == 0:
                    # no candidate left: the task is not yet tracked.
                    break
            if len(candidates) > 1:
                raise Exception("DB consistency error: several candidate tasks were found.")
            elif len(candidates) == 1:
                step_concrete_id = next(iter(candidates))
            else:
                # step_concrete_id remains None
                if step_concrete_id is not None:
                    raise Exception('')

        if step_concrete_id is not None:
            # the concrete step is already in the database
            res = DbID(step_concrete_id, False)
            # The following is a check of the targets. 
            target_entities = tuple(self.get_targetassets(res))
            if len(target_entities) != len(targets):
                raise ValueError(''.join(('The task was already stored but with a different number of targets ',
                                          '(%i, but now %i).' % (len(target_entities), len(targets)))))
                
            for entity in target_entities:
                targetasset = getattr(targets, entity.label)
                if targetasset._defined:
                    # Check that this is matching the same as what is in the database
                    # Being thorough is relatively important here
                    if isinstance(entity, StoredSequence):
                        raise NotImplementedError('Sequence assets in targets is not yet handled.')
                    # check that constructor's attribute is the same:
                    if targetasset.name != entity.entityname or type(targetasset).__name__ != entity.clsname:
                        raise ValueError(''.join(('The task was already stored but with different values for the target "%s"' % entity.label,
                                                  ' (%s(%s) instead of %s(%s)).' % (entity.clsname, entity.entityname,
                                                                                    type(targetasset).__name__, targetasset.name))))
                else:
                    if isinstance(entity, StoredSequence):
                        #
                        raise NotImplementedError('Sequence assets in targets is not yet handled.')
                    else:
                        # set it to what is in the database
                        targetasset.name = entity.entityname
        else:
            # the concrete step was never seen before
            param_json = json.dumps(parameters_id)
            res = self._add_stepconcrete(step_variant_id,
                                         sources,
                                         targets,
                                         parameters_id,
                                         tag = tag)
        return res
        
    def _add_stepconcrete(self, step_variant_id,
                          sources,
                          targets,
                          parameters_id,
                          tag = 1):
        """
        Add a "concrete" step, that is a step variant (executable and parameters)
        to which source and target files are added.

        :param step_variant_id: ID for the step variant
        :type step_variant_id: integer
        :param sources: sequence of sources
        :type sources: :class:`AssetSet`
        :param targets: sequence of targets
        :type targets: :class:`AssetSet`
        :param parameters_id: ID for the parameters
        :type tag: A tag for the task
        :param tag: :class:`int`
        :rtype: :class:`DbID`

        """

        assert isinstance(step_variant_id, int)

        cursor = self.connection.cursor()
        #FIXME: transaction (in case one of the inserts fails) only committed
        #       before exiting this method, but can this be bypassed by SQLite global settings ?
        #       (I am relying on the current default)
                
        sql = """
        INSERT INTO step_concrete (
        step_variant_id,
        step_status_id,
        time_creation,
        tag
        ) VALUES (
        ?,
        ?,
        ?,
        ?);
        """
        # last update time set to NULL

        #FIXME: try/catch. since step_variant_id is a foreign key,
        #       this will fail if there is no such step_variant_id
        cursor.execute(sql, (step_variant_id, _TASK_STATUS_LIST[_TASK_TODO], time.time(), tag))
        step_concrete_id = DbID(cursor.lastrowid, True)

        labelnamepairs = zip(sources._fields, sources)
        self._insert_stepconcrete2storedentities(labelnamepairs, 'src', step_concrete_id.id)

        # ensure that defined target are tracked
        labelnamepairs = list()
        for field, t in zip(targets._fields, targets):
            if t._defined:
                labelnamepairs.append((field, t))
        self._insert_stepconcrete2storedentities(labelnamepairs, 'target', step_concrete_id.id)
        self.connection.commit()

        self._insert_stepconcrete2parameters(parameters_id.id, step_concrete_id.id)
        self.connection.commit()
        return step_concrete_id
        #return Step(step_concrete_id, sources, targets)

    def _insert_stepconcrete2parameters(self, parameters_id, step_concrete_id):
        """ WARNING: add a commit after the call !!!"""
        cursor = self.connection.cursor()

        sql = """
        INSERT INTO step_concrete2parameters (
        step_concrete_id,
        step_parameters_id
        ) VALUES (
        ?,
        ?);
        """

        cursor.execute(sql,
                       (step_concrete_id, parameters_id))

    def _insert_stepconcrete2storedentities(self, labelnamepairs, what, step_concrete_id):
        """ WARNING: add a commit in caller after this returns !!!"""
        assert what in ('src', 'target')
        cursor = self.connection.cursor()

        sql_template = """
        INSERT INTO step_concrete2%(what)sfile (
        label,
        step_concrete_id,
        stored_%(the)s_id
        ) VALUES (
        ?,
        ?,
        ?);
        """
        # 2 SQL query: one for assets that are sequences and one for assets that are scalars.
        sql_ss = sql_template % {'what': what, 'the': 'sequence'}
        sql_se = sql_template % {'what': what, 'the': 'entity'}

        for (label, asset) in labelnamepairs:
            if isinstance(asset, core.FileSequence):
                # if a sequence we also need to create the sequence itself
                sequence_id = self.id_stored_sequence(type(asset), asset.iteritems()).id
                se_id = None
                etid = (label, 
                        step_concrete_id,
                        # stored_entity
                        sequence_id)
                sql = sql_ss
            else:
                sequence_id = None
                if asset is None:
                    # if the asset is not defined here, 
                    # it can / should only be the case
                    # because it is allowed to be None
                    continue
                    #se_id = self.id_stored_entity(type(asset), None).id
                else:
                    for se_i, cn in enumerate(asset.iteritems()):
                        if se_i > 0:
                            raise Exception("Only FileSequence object should have several saved entities.")
                        # ensure the that the stored IDs are tracked in the DB
                        se_id = self.id_stored_entity(*cn).id
                etid = (label, 
                        step_concrete_id, 
                        se_id
                        # stored_sequence
                )
                sql = sql_se
            res = cursor.execute(sql, etid)
            
        #labelstepandnames = tuple( for cn in asset.iteritems())
        #cursor.executemany(sql,
        #labelstepandnames)

    def step_concrete_info(self, step_concrete_id):
        sql = """
        SELECT
          step_variant_id,
          step_status.label,
          time_creation
        FROM  step_concrete, step_status
        WHERE step_concrete.id=?
        AND step_status.id = step_status_id
        """
        cursor = self.connection.cursor()
        cursor.execute(sql, (step_concrete_id.id, ))
        res = cursor.fetchone()
        return res

    def step_concrete_state(self, step_concrete_id,
                            state_id):
        """ Set the state of a task:

        - step_concrete_id: task ID (DbID)
        - state_id:         state ID
        """
        sql = """
        UPDATE step_concrete
        SET step_status_id=?
        WHERE step_concrete.id=?
        """
        cursor = self.connection.cursor()
        res = cursor.execute(sql, (state_id, step_concrete_id.id))
        self.connection.commit()
        #FIXME: return anything ?

    def step_concrete_status(self, step_concrete_id):
        sql = """
        SELECT step_status.id, step_status.label
        FROM step_concrete
        INNER JOIN step_status
        ON step_status_id=step_status.id
        WHERE step_concrete.id=?
        """
        cursor = self.connection.cursor()
        cursor.execute(sql, (step_concrete_id.id,))
        return cursor.fetchall()

    def _get_assetsofactivity(self, activity, what):
        """ 
        Find the assets associated with steps performing a specific activity.
        :param activity: an activity
        :type activity: :class:`Enum`
        """
        cursor = self.connection.cursor()

        # -- find step type IDs performing a specific activity
        sql = """
        SELECT DISTINCT step_type_id
        FROM step_activity, step_type2activity
        WHERE label=?
        AND step_activity.id=step_type2activity.step_activity_id
        """
        res = cursor.execute(sql, (activity.value,))
        #FIXME: check that activity with that label found ?
        step_type_ids = res.fetchall()

        # -- find concrete step IDs with that activity
        sql = """
        SELECT step_concrete.id
        FROM step_variant, 
             step_concrete,
             step_concrete2%(what)sfile
        WHERE step_type_id=?
        AND step_variant.id=step_variant_id
        AND step_concrete.id=step_concrete2%(what)sfile.step_concrete_id
        """ % {'what': what}
        step_concrete_ids = list()
        for step_type_id in step_type_ids:
            res = cursor.execute(sql, step_type_id)
            step_concrete_ids.extend(res.fetchall())

        # -- find the stored entities resulting from the concrete steps above
        results = list()
        for task_dbentry in step_concrete_ids:
            res = self.get_targetassets(task_dbentry[0])
            results.extend(res)

        # -- 
        
        return tuple(results)
        
    def get_sourcesofactivity(self, activity):
        """ 
        Retrieve the sources of steps performing a specific activity.
        :param activity: an activity
        :type activity: :class:`Enum`
        """
        return self._get_assetsofactivity(activity, 'src')

    def find_targetsofactivity(self, activity):
        warnings.warn('find_targetsofactivity() is deprecated, use get_targetsofactivity().')
        return self._get_assetsofactivity(activity, 'target')

    def get_targetsofactivity(self, activity):
        """ 
        Retrieve the targets of steps performing a specific activity.
        :param activity: an activity
        :type activity: :class:`Enum`
        """
        cursor = self.connection.cursor()

        # -- find step type IDs performing a specific activity
        sql = """
        SELECT DISTINCT step_type_id
        FROM step_activity, step_type2activity
        WHERE label=?
        AND step_activity.id=step_type2activity.step_activity_id
        """
        res = cursor.execute(sql, (activity.value,))
        step_type_ids = res.fetchall()
        if len(step_type_ids) == 0:
            sql = """
            SELECT DISTINCT id
            FROM step_activity
            WHERE label=?
            """
            res = cursor.execute(sql, (activity.value,))
            step_activity_ids = res.fetchall()
            if len(step_activity_ids) == 0:
                raise ValueError("The activity '%s' is not in the database." % activity.value)
            else:
                return tuple()
            
        # -- find task IDs (concrete steps) with that activity
        sql = """
        SELECT DISTINCT step_concrete.id
        FROM step_variant, 
             step_concrete,
             step_concrete2targetfile
        WHERE step_type_id=?
        AND step_variant.id=step_variant_id
        AND step_concrete.id=step_concrete2targetfile.step_concrete_id
        """
        step_concrete_ids = list()
        for step_type_id in step_type_ids:
            res = cursor.execute(sql, step_type_id)
            step_concrete_ids.extend(res.fetchall())

        # -- find the stored entities created by the task above
        results = list()
        for task_dbentry in step_concrete_ids:
            res = self.get_targetassets(task_dbentry[0])
            results.extend(res)
        
        return tuple(results)

    def get_targetsoftype(self, clsname):
        """ Return all targets of a given type. """
        cursor = self.connection.cursor()
        sql = """
        SELECT stored_entity.id,
               classname,
               entityname
        FROM stored_entity
        INNER JOIN step_concrete2targetfile
        ON stored_entity.id=step_concrete2targetfile.stored_entity_id
        WHERE stored_entity.classname=?
        """
        res = cursor.execute(sql, (clsname,))
        results = tuple(self.StoredEntityNoLabel(*x) for x in res.fetchall())
        return results

    def iter_finaltargets(self):
        """ Targets not used as source anywhere else. """
        cursor = self.connection.cursor()
        sql = """
        SELECT se.id, se.classname, se.entityname, step_concrete_id, label
        FROM stored_entity AS se
        INNER JOIN (
          -- select targets not used used as source
          SELECT trg.stored_entity_id, trg.step_concrete_id
          FROM step_concrete2targetfile AS trg
          LEFT JOIN step_concrete2srcfile AS src
          ON trg.stored_entity_id=src.stored_entity_id
          WHERE src.stored_entity_id IS NULL
          )
        ON se.id=stored_entity_id
        INNER JOIN step_concrete
        ON step_concrete.id=step_concrete_id
        INNER JOIN step_status
        ON step_status_id=step_status.id
        ORDER BY se.id
        """
        cursor.execute(sql)
        Target = namedtuple('Target',
                            'stored_entity_id stored_entity_classname stored_entity_entityname step_concrete_id status_label')
        for x in cursor:
            yield Target(*x)

    def finalsteps(self):
        """ Concrete steps for which all targets are final """
        # concrete steps for which all targets are final
        raise NotImplementedError()



class CachedPersistentTaskGraph(PersistentTaskGraph):

    def __init__(self, db_fn, model, wd='.', force_create=False, isolation_level=None):
        super(CachedPersistentTaskGraph, self).__init__(db_fn, model, wd=wd, 
                                                        force_create=force_create,
                                                        isolation_level=isolation_level)
        # cache of DB IDs
        self._cache_steptype_dbid = dict()
        self._cache_task_dbid = dict()
        self._cache_dbid_task_targets = dict()
        self._cache_stepvariant_dbid = dict()
        self._cache_dbid_stepvariant = dict()
        self._cache_stepparameters_dbid = dict()
        self._cache_storedentity_dbid = dict()
        self._cache_storedsequence_dbid = dict()

    def id_step_type(self, activities):
        hashdb = activities
        dbid = self._cache_steptype_dbid.get(hashdb)
        if dbid is not None:
            res = DbID(dbid, False)
        else:
            res = super(CachedPersistentTaskGraph, self).id_step_type(activities)
            self._cache_steptype_dbid[hashdb] = res.id
        return res

    def _cache_id_step_type(self):
        for dbid, activities in itertools.groupby(operator.itemgetter(0),
                                                  self._iter_step_type()):
            hashdb = tuple(activities)
            self._cache_steptype_dbid[hashdb] = dbid

    def id_step_variant(self, step, activities):
        cls = type(step)
        #FIXME: redundancy between the full class name and the activities)
        step_hashdb = (cls.__module__ + '.' + cls.__name__, activities)
        variant_hashdb = (step.version, step._execpath)
        hashdb = (step_hashdb, variant_hashdb)
        dbid = self._cache_stepvariant_dbid.get(hashdb)
        if dbid is not None:
            res = DbID(dbid, False)
        else:
            res = super(CachedPersistentTaskGraph, self).id_step_variant(step, activities)
            self._cache_stepvariant_dbid[hashdb] = res.id
            self._cache_dbid_stepvariant[res.id] = (step_hashdb, variant_hashdb)
        return res


    def id_stepconcrete(self, step_variant_id,
                        sources, targets, parameters,
                        tag = 1):
        """
        Conditionally add a task ("concrete" step),
        that is a step variant (executable and parameters)
        to which source and target files, as well as parameters, are added.

        :param step_variant_id: ID for the step variant
        :type step_variant_id: integer
        :param sources: sequence of sources
        :type sources: :class:`AssetSet`
        :param targets: sequence of targets
        :type targets: :class:`AssetSet`
        :param parameters: list of parameters
        :type parameters: a sequence of :class:`str`
        :param tag: a tag, used to performed repetitions of the exact same task
        :type tag: a sequence of :class:`int`

        :rtype: :class:`DbID`

        """

        assert isinstance(step_variant_id, int)

        step_variant_hashdb = self._cache_dbid_stepvariant.get(step_variant_id)
        if step_variant_hashdb is None:
            raise ValueError('There is no step variant with ID %i' % step_variant_id)
        step_hashdb, variant_hashdb = step_variant_hashdb
        task_hashdb = (step_hashdb, sources.hashdb, parameters, tag)
        dbid = self._cache_task_dbid.get(task_hashdb)
        if dbid is not None:

            res = DbID(dbid, False)

            # Some of the targets might be undefined (because the user leaves it to railroadtracks to
            # fill the blanks). We need to populate these with the values in the database in order to
            # be able to compute the hash.
            if any(not getattr(targets, x.label)._defined for x in self.get_targetassets(res)):
                # FIXME: for now just hit the database (no use of the cache)
                res = super(CachedPersistentTaskGraph, self).id_stepconcrete(step_variant_id, 
                                                                             sources, 
                                                                             targets,
                                                                             parameters, 
                                                                             tag=tag)
            elif dbid in self._cache_dbid_task_targets:
                # check that the target assets have not changed
                if targets.hashdb != self._cache_dbid_task_targets[dbid]:
                    raise ValueError("Target assets not matching the assets already associated with the task.")
            else:
                self._cache_dbid_task_targets[res.id] = targets.hashdb
            
        else:
            res = super(CachedPersistentTaskGraph, self).id_stepconcrete(step_variant_id, sources, targets, parameters, tag=tag)
            self._cache_task_dbid[task_hashdb] = res.id
            if all(x._defined for x in targets):
                self._cache_dbid_task_targets[res.id] = targets.hashdb
        return res


    def id_stepparameters(self, parameters):
        
        dbid = self._cache_stepparameters_dbid.get(parameters)
        if dbid is not None:
            res = DbID(dbid, False)
        else:
            res = super(CachedPersistentTaskGraph, self).id_stepparameters(parameters)
            self._cache_stepparameters_dbid[parameters] = res.id
        return res


    def id_stored_entity(self, cls, name):
        hashdb = core.SavedEntityAbstract._hash_components(cls, name)
        dbid = self._cache_storedentity_dbid.get(hashdb)
        if dbid is not None:
            res = DbID(dbid, False)
        else:
            res = super(CachedPersistentTaskGraph, self).id_stored_entity(cls, name)
            self._cache_storedentity_dbid[hashdb] = res.id
        return res
        

    def id_stored_sequence(self, cls, clsname_sequence):
        names_hashdb = list()
        clsname_sequence_tpl = tuple(clsname_sequence)
        for (elt_cls, elt_name) in clsname_sequence_tpl:
            dbid = self.id_stored_entity(elt_cls, elt_name) 
            item_hashdb = core.SavedEntityAbstract._hash_components(elt_cls, elt_name)
            names_hashdb.append(item_hashdb)
        hashdb = (cls, tuple(names_hashdb))
        dbid = self._cache_storedsequence_dbid.get(hashdb)
        if dbid is not None:
            res = DbID(dbid, False)
        else:
            res = super(CachedPersistentTaskGraph, self).id_stored_sequence(cls, clsname_sequence_tpl)
            self._cache_storedsequence_dbid[hashdb] = res.id
        return res
        
