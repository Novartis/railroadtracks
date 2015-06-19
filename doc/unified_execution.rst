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

.. _unified-execution-label:

Unified execution
=================

Almost each tool in a sequence of steps used for RNA-Seq processing is using
idiosyncratic parameters or arguments. This is making the task of a person
wanting to use them both prone to errors (as these tools are also not always
checking thoroughly that the input and parameters make sense) and to
a significant time spent reading scattered sources of information
(the documentation for the tools if often insufficient and much of the
the knowledge about them is spread throughout internet forums and mailing-lists)

.. -- note-unifex-begin

.. code-block:: bash

   $ alias unifex='python -m railroadtracks.unifex'

.. -- note-unifex-end

Three modes exist: `run`, `version`, and `activities`

.. code-block:: bash

   $ unifex
   usage: unifex.py [-h] {run,version,activities} ...
   rnaseq.py: error: too few arguments

A step is then defined by the name of the executable (either as a name to 
be found in the `${PATH}`, or as an absolute path), as well as by the
name of a modeling class.

.. note::

   Specifying both the executable and the modeling class is required because
   a number of executables can perform different tasks/activities.

   For example, the command STAR can be used to build an index reference
   for subsequent alignment, or perform an alignment. Specifying the model
   associated with the executable helps to remove the ambiguity.

   .. code-block:: bash

      $ unifex run star-index STAR
      $ unifex run star-align STAR # astrology anyone ?



Version number
--------------

Obtaining the version number is achieved the same way for all steps


.. -- note-unifex-version-begin

.. code-block:: bash

   $ unifex version bowtie2-build bowtie2
   2.1.0
   $ unifex version star-align STAR
   STAR_2.3.0e_r291
   $ unifex version limma-voom R
   3.18.3
   $ unifex version edger R
   3.4.0

.. -- note-unifex-version-end

Running a step
--------------

Running a step is also achieved essentially the same way for all steps.
All steps have "sources" (that is input files) and "targets"
(that is destination/output files).


.. -- note-unifex-run-begin

.. code-block:: bash

   $ # bowtie2 (create index)
   $ unifex run bowtie2-build bowtie2 \
         -s <name>=<source file(s)> \
         -t <name>=<target file(s)>
   $ # STAR (create index)
   $ unifex run star-index STAR \
         -s <name>=<source file(s)> \
         -t <name>=<target file(s)>

.. -- note-unifex-run-end

.. note::

   Whenever the sources and targets expected by a given tool are not specified,
   the bash command fails and print the list of missing parameters

   .. code-block:: bash

      $ unifex run edger R
      The following sources must be defined (and are missing):
      - counttable_fn
      - sampleinfo_fn
      The following targets must be defined (and are missing):
      - diffexp_fn


Using a scheduler/Queueing system
---------------------------------

The persistent layer can generate the bash commands for running the unified execution.
This is making the use of any existing queuing or scheduling system able to take `bash` script
straightforward.

An implementaion using SGE's qsub is provided with :mod:`railroadtracks.easy.qsub`.

Docstrings
----------

.. automodule:: railroadtracks.unifex
   :members:
   :special-members: __init__
