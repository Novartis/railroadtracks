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

import networkx
from collections import namedtuple, deque
from railroadtracks import easy
from railroadtracks.easy.execution import exec_taskfilter

TaskSetBundle = namedtuple('TaskSetBundle', 'label taskset mapper filter')

def get_tasksetid(taskset):
    """ Utility function that returns a valid (hashable) key for a task set.
    """
    return frozenset(x.task_id for x in taskset)

class TaskSetGraph(easy.TaskSet):
    """
    Graph of TaskSet nodes.
    """
    def __init__(self, project=None, defaultmapper=None, 
                 defaulttasksetfilter=exec_taskfilter):
        """
        :param project: 
        :param defaultmapper: default task mapper
        :param defaulttasksetfilter: default filter to apply to a taskset before passing it
                                     to the task mapper. This is not used for the moment.
        """
        super(TaskSetGraph, self).__init__(project=project)
        self._digraph = networkx.digraph.DiGraph()
        self._taskid2tasksetid = dict()
        self.defaultmapper = defaultmapper
        self.defaulttasksetfilter = defaulttasksetfilter

    @property
    def digraph(self):
        return self._digraph

    def _ensure_task_unique(self, taskset):
        intersect = self.intersection(taskset)
        if len(intersect) > 0:
            if len(intersect) > 1:
                plural = ('s', 'are')
            else:
                plural = ('', 'is')
            raise ValueError('%i of the task%s in the task set %s already in the graph.' % (len(intersect), plural[0], plural[1]))
    

    def add_taskset(self, taskset, label=None, mapper=None, filter=None):
        """ Add a TaskSet to the graph. 
        The dependency relationships with other :class:`TaskSet` objects
        are inferred from the assets in the individual :class:`Task` objects.

        :param taskset: a :class:`easy.tasks.TaskSet`
        :param label: a label for the `taskset`
        :param mapper: a task mapper to use with this taskset
        :param filter: a task filter to use this this taskset

        """
        self._ensure_taskset_project(taskset)
        self._ensure_task_unique(taskset)
        tasksetid = get_tasksetid(taskset)
        # make a copy to prevent problems if further changes in the TaskSet performed later
        # in the enclosing frame
        taskset_copy = easy.TaskSet()
        if mapper is None:
            mapper = self.defaultmapper
        self._digraph.add_node(tasksetid, attr_dict={'label': label,
                                                     'taskset': taskset_copy,
                                                     'mapper': mapper,
                                                     'filter': filter})
        for task in taskset:
            self._taskid2tasksetid[task.task_id] = tasksetid
            taskset_copy.add(task)
            super(TaskSetGraph, self).add(task)
            # parent tasks
            taskqueue = deque((x, tasksetid) for x in task.parent_tasks())
            while len(taskqueue) > 0:
                parent_task, child_tasksetid = taskqueue.pop()
                if parent_task.task_id in self._taskid2tasksetid:
                    # taskset in which the parent is present
                    p_tasksetid = self._taskid2tasksetid[parent_task.task_id]
                    if not self._digraph.has_edge(p_tasksetid, child_tasksetid):
                        self._digraph.add_edge(p_tasksetid, tasksetid)
                        # add the parents to the queue to reconstruct the graph of tasksets
                        for p in parent_task.parent_tasks():
                            taskqueue.appendleft((p, p_tasksetid))
                else:
                    pass
            # child tasks
            # (note: taskqueue is empty again when reaching here)
            taskqueue = deque((tasksetid, x) for x in task.child_tasks())
            while len(taskqueue) > 0:
                parent_tasksetid, child_task = taskqueue.pop()
                if child_task.task_id in self._taskid2tasksetid:
                    # taskset in which the child is present
                    c_tasksetid = self._taskid2tasksetid[child_task.task_id]
                    if not self._digraph.has_edge(tasksetid, c_tasksetid):
                        self._digraph.add_edge(tasksetid, c_tasksetid)
                        # add the children to the queue to reconstruct the graph of tasksets
                        for c in child_task.child_tasks():
                            taskqueue.appendleft((tasksetid, c))
                else:
                    pass
        return tasksetid

    def add(self, obj, label=None, mapper=None, filter=None):
        """ Generic add method, allowing the addition of either a Task or TaskSet object"""
        if isinstance(obj, easy.Task):
            obj = easy.TaskSet(iterable=(obj,))
        assert isinstance(obj, easy.TaskSet)
        self.add_taskset(obj, label=label, mapper=mapper)


    def execution_list(self):
        """ Return a list of :class:`TaskSetBundle` objects in an order in
        which they can be executed. """
        g = self._digraph
        return tuple(TaskSetBundle(g.node[x]['label'],
                                   g.node[x]['taskset'],
                                   g.node[x]['mapper'],
                                   g.node[x]['filter']) for x in networkx.topological_sort(g))

    def _filtered_taskset(self, tsb):
        if tsb.filter is None:
            taskset_filter = self.defaulttasksetfilter
        else:
            taskset_filter = tsb.taskset_filter
        if taskset_filter is None:
            taskset = tsb.taskset
        else:
            taskset = taskset_filter(tsb.taskset)
        return taskset
        
    def execute(self, 
                update_func=None):
        """ As with the method of the same name in the the parent class
        this will execute the tasks contained, with the notable difference
        that the execution that the tasks will be grouped according to the
        nested task sets defined and in the implicit order defined by the
        connections at the task level.

        This method is little beyond a wrapper around the method
        execution_list(), iterating through the task sets returned and
        executing the tasks with either the given task mapper of the default
        task mapper. 
        """
        lst = self.execution_list()
        for i, tsb in enumerate(lst):
            if tsb.mapper is None:
                mapper = self.defaultmapper
            else:
                mapper = tsb.mapper
            if mapper is None:
                raise ValueError('No task mapper was specified when the task was added, ' + \
                                 'and the default task mapper for the %s is not specified' % type(self).__name__)
            taskset = self._filtered_taskset(tsb)
            if update_func is not None:
                update_func(i, lst, taskset)
            mapper.map(taskset)

def update_func(i, lst, filtered_taskset):
    print('%i - %s (%i/%i tasks)' % (i, lst[i].label, len(filtered_taskset), len(lst[i].taskset)))
