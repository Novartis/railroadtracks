..
   Copyright 2015 Novartis Institutes for Biomedical Research

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

.. _cookbook-label:

********
Cookbook
********

Environment
===========

Checking that 3rd-party executables are present
-----------------------------------------------

When writing a recipe, checking that 3rd-party executables are present can be
highly desirable in order to let a script fail early and provide informative message
to a user about what is required but missing.

.. code-block:: python

   from railroadtracks.environment import Executable, MissingSoftware

   if not Executable.ispresent('bwa'):
       raise MissingSoftware('The bwa is not in the PATH.')


Project
=======

Get all tasks of a given type
-----------------------------

The types, or classes, defined in a model can be used to query a project.

.. code-block:: python

   taskset = project.get_tasksoftype(rnaseq.StarAlign)


Get all tasks matching a list of types
--------------------------------------

.. code-block:: python

   listoftypes = (rnaseq.StarAlign, rnaseq.BWA)
   taskset = reduce(lambda x,y: x.union(project.get_tasksoftype(y)),
                    listoftypes, easy.TaskSet())


Get all tasks fulfilling a given activity
-----------------------------------------

.. code-block:: python

   taskset = project.get_taskswithactivity(rnaseq.ACTIVITY.ALIGN)



TaskSet
=======

:class:`TaskSet` objects provide a convenient abstraction to manipulate groups of tasks by thinking about them
as :class:`set` objects. In addition to :class:`set` methods, they also have :mod:`railroadtracks` specific features.

Subset all tasks with a given status
------------------------------------

Retrieving all tasks with a given status from a :class:`TaskSet` can be achieved with the method
:meth:`TaskSet.filter_on_status`.

For example:

.. code-block:: python

   # all tasks in the set that failed
   taskset_failed = taskset.filter_on_status(hortator._TASK_FAILED))

   # all tasks in the set that succeeded
   taskset_success = taskset.filter_on_status(hortator._TASK_TODO))


Subset all tasks ready for execution
------------------------------------

Tasks ready for execution have all their /parent/ tasks with a status `_TASK_DONE` and have themselves
the status `_TASK_TODO` or `_TASK_FAILED`.

This is implemented very simply, using :meth:`TaskSet.filter_on_status`, :meth:`TaskSet.union`,
and :meth:`TaskSet.filter_on_parent_status`.

.. code-block:: python

   def exec_taskfilter(taskset):
       # union of tasks that are either "TO-DO" or "FAILED"
       tasks_todo = taskset.filter_on_status(hortator._TASK_TODO)
       tasks_todo = tasks_todo.union(taskset.filter_on_status(hortator._TASK_FAILED))
       # only keep the tasks for which the parent task was successfully run
       tasks_todo = tasks_todo.filter_on_parent_status(hortator._TASK_DONE)
       return tasks_todo


That function is part of the code base, so in practice one will only have to write:

.. code-block:: python

   import railroadtracks.easy.execution
   ts_torun = railroadtracks.easy.execution.exec_taskfilter(taskset)



Get all child tasks for a :class:`TaskSet`
------------------------------------------

Each :class:`TaskSet` can be manipulated as a :class:`set`, so if we think of the problem
as the union of the child taks for each task in the set it writes simply as:

.. code-block:: python

   ts_allchilds = easy.TaskSet()
   for task in taskset:
       ts_allchilds.union(task.child_tasks())

The one-line version is:

.. code-block:: python

   ts_allchilds = reduce(lambda x,y: x.union(y),
                         map(lambda x: x.child_tasks(), taskset),
                         easy.TaskSet())


TaskSetGraph
============

:class:`TaskSet` objects are rarely only considered in isolation since they are part of a larger
dependency graph. :class:`TaskSetGraph` helps handling all the tasks by allowing one to manipulate
groups of tasks not connected directly as :class:`TaskSet` objects and in the context
of the task-level dependency graph. Such sets can be added directly to a :class:`TaskSetGraph`
while dependencies between the task sets are inferred from task-level relationships.

.. code-block:: python

   from railroadtracks.easy.tasks import TaskSet
   from railroadtracks.easy.tasksetgraph import TaskSetGraph

   tsg = TaskSetGraph()

   # add a taskset
   tsg.add(taskset_a)

   # add an other taskset
   tsg.add(taskset_b)


Execute all tasks
-----------------

In order to execute tasks, one can decide on a default task mapper (task mappers can also be set at the
:class:`TaskSet` level) and on a task filter (a function that decides which tasks should be executed).
The order in which :class:`TaskSet` object should be executed is inferred from the connections between
the :class:`Tasks` in the :class:`TaskSet` and does not require a user's intervention. The task filter
provides an additional level of refinement by letting a function control the tasks in the set that are
executed.

In the example below, we are using :class:`easy.execution.IterativeExecution` to map the tasks, and
the filter introduced earlier in the section `TaskSet`: only consider for execution the tasks that are
not yet done, and for which the /parent/ tasks were successfully completed (no point trying to
execute a task if its parent has not been succesful).


.. code-block:: python

   from railroadtracks.easy.processing import IterativeExecution   
   # task mapper to run the tasks
   p = IterativeExecution()
   tsg.defaultmapper = p

   # task filter
   tsg.defaultfilter = easy.execute.exec_taskfilter

   # execute all tasks that pass the filter
   tsg.execute()

