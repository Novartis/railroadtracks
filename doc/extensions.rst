..
   Copyright 2014 Novartis Institutes for Biomedical Research

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

Extending the framework
=======================


Writing custom steps
--------------------

Simple case
^^^^^^^^^^^

A lot of the boilerplate handling is handled by parent classes. A custom class for a new tool will have
to implement code to:

- create an instance (constructor)
- get the version
- run the step

The base class for steps is the abstract class :class:`railroadtracks.core.AssetStep`:

.. literalinclude:: ../src/core.py
   :language: python
   :start-after: # -- stepabstract-begin
   :end-before: # -- stepabstract-end

Custom steps will just have to implement the methods and properties.
For example, counting the reads after alignment with htseq-count is defined as a
child class of :class:`railroadtracks.model.quantify.QuantifierAbstract`, where :attr:`Assets` are defined:

.. literalinclude:: ../src/model/quantify.py
   :language: python
   :start-after: # -- quantifierabstract-begin
   :end-before: # -- quantifierabstract-end

.. note:: 

   There is a class for each activity defined (see :ref:`sec-activities-label`), 
   and these can be used as parent
   class for new steps performing one activity
   The :attr:`Assets` is inherited from its parent class, and does not need to be specified further.

   .. literalinclude:: ../src/model/quantify.py
      :language: python
      :start-after: # -- assetsquantifier-begin
      :end-before: # -- assetsquantifier-end

   One can note that specific subclasses of :class:`railroadtracks.core.SavedEntityAbstract` can be specified
   (here :class:`CSVFile` to indicate that file produced by the counting step is a CSV file).
   The definition of assets represents a way to add a static typing flavor to Python,
   as a number of checks are performed behind the hood by :mod:`railroadtracks.model.files`.


The definition of :class:`railroadtracks.model.quantify.HTSeqCount` is then implementing
the remaining methods and attributes

.. literalinclude:: ../src/model/quantify.py
   :language: python
   :start-after: # -- htseqcount-begin
   :end-before: # -- htseqcount-end



Steps performing several activities
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Slightly more work than for the simple case will be required, but not too much either.


Unified execution
-----------------

A new script, extending or overriding the default execution layer in
the package can be written very simply:

.. code-block:: python

   from railroadtracks import core, rnaseq

   class MyStep(core.Stepabstract):
      # Definition of a new step
      pass

   NEW_STEPLIST_CLASSES = list()
   for cls in rnaseq._STEPLIST_CLASSES:
       NEW_STEPLIST.append(cls)
   NEW_STEPLIST_CLASS.append(MyStep)

   if __name__ == '__main__':
       # use the command and subcommands parsing
       # with the new list of steps
       core.unified_exec(classlist = NEW_STEPLIST_CLASSES)


