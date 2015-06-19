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

from railroadtracks import hortator
from railroadtracks.easy.tasks import Task
from collections import namedtuple, deque
import networkx

Node = namedtuple('Node', 'clsname id')


def _task(id, project, cursor):
    return project.get_task(id)

def _asset(id, project, cursor, assettype):
    """ Create a StoredEntity from its ID. """
    sql = """
    SELECT id, label, classname, entityname
    FROM stored_entity
    JOIN step_concrete2%sfile
    ON stored_entity.id=stored_entity_id
    WHERE stored_entity_id=?
    """ % assettype
    cursor.execute(sql, (id, ))
    dbid, label, clsname, entityname = cursor.fetchone()        
    return hortator.StoredEntity(dbid, label, clsname, entityname)

def _assetsequence(id, project, cursor, assettype):
    """ Create a StoredSequence from its ID. """
    sql = """
    SELECT stored_sequence.id, label, stored_sequence.classname, se2s.stored_entity_id
    FROM stored_sequence
    JOIN stored_entity2sequence as se2s
    ON se2s.stored_sequence_id=stored_sequence.id
    JOIN step_concrete2%sfile as sc2sf
    ON sc2sf.stored_sequence_id=se2s.stored_sequence_id
    WHERE se2s.stored_sequence_id=?
    ORDER BY pos
    """ % assettype
    cursor.execute(sql, (id, ))
    storedentities = list()
    for dbid, label, clsname, assetid in cursor:
        storedentities.append(assetid)
            
    return hortator.StoredSequence(dbid, label, clsname, storedentities)

nodeclsname2constructor = {'Task': _task,
                           'Asset': _asset,
                           'AssetSequence': _assetsequence}

def get_digraph(project):
    digraph = networkx.digraph.DiGraph()
    #
    cursor = project.persistent_graph.connection.cursor()
    sql_tasks_template = """
    SELECT step_concrete.id,
           stored_entity_id,
           stored_sequence_id
    FROM step_concrete,
         step_concrete2%sfile
    WHERE step_concrete2%sfile.step_concrete_id=step_concrete.id
    """
    for assettype in ('src', 'target'):
        sql_tasks = sql_tasks_template % (assettype, assettype)
        cursor.execute(sql_tasks)
        #
        cursor_nodes = project.persistent_graph.connection.cursor()
        for taskid, assetid, assetsequenceid in cursor:
            nodeclsname = 'Task'
            task_node = Node(nodeclsname, taskid)
            if task_node not in digraph:
                instance = nodeclsname2constructor[nodeclsname](taskid, project, cursor_nodes)
                digraph.add_node(task_node, instance=instance)
            # sanity check
            assert (assetid is None) != (assetsequenceid is None), \
                                        'Issue with the database: it is either an asset or a sequence of assets.'
            if assetsequenceid is None:
                nodeclsname = 'Asset'
                dbid = assetid
                asset_node = Node(nodeclsname, dbid)
            else:
                # no need to test because of the sanity check above
                nodeclsname = 'AssetSequence'
                dbid = assetsequenceid
                asset_node = Node(nodeclsname, dbid)
            if asset_node not in digraph:
                instance = nodeclsname2constructor[nodeclsname](dbid, project, cursor_nodes, assettype)
                digraph.add_node(asset_node, instance=instance)
            if assettype == 'src':
                digraph.add_edge(asset_node, task_node)
            else:
                digraph.add_edge(task_node, asset_node)
    #
    sql_seq = """
    SELECT se2s.stored_sequence_id as assetsequenceid,
           stored_entity.id as assetid,
           pos
    FROM stored_entity2sequence AS se2s
    JOIN stored_entity
    ON se2s.stored_entity_id=stored_entity.id
    ORDER BY assetsequenceid, pos
    """
    cursor.execute(sql_seq)
    for assetsequenceid, assetid, pos in cursor:
        assetsequence_node = Node('AssetSequence', assetsequenceid)
        asset_node = Node('Asset', assetid)
        if asset_node not in digraph:
            instance = nodeclsname2constructor[nodeclsname](assetid, project, cursor_nodes, 'Asset')
            digraph.add_node(asset_node, instance=instance)
        digraph.add_edge(asset_node, assetsequence_node)
    #

    return digraph

class TaskGraph(object):
    
    def __init__(self, project):
        self._project = project
        self._digraph = get_digraph(project)

    @property
    def digraph(self):
        return self._digraph

    def predecessor_tasks(self, task):
        if isinstance(task, Task):
            node = Node('Task', task.task_id)
        else:
            node = task
        digraph = self._digraph
        assets = digraph.predecessors(node)
        ptasks = list()
        while len(assets) > 0:
            a = assets.pop()
            if a.clsname == 'AssetSequence':
                assets.extend(digraph.predecessors(Node('AssetSequence', a.id)))
                continue
            for p in digraph.predecessors(a):
                assert p.clsname == 'Task'
                ptasks.append(p)
        return ptasks
