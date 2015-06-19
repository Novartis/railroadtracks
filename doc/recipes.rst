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

.. _recipes-label:

Recipes
=======


The model, together with the persistence layer, is designed to make the writing
of sequences of steps for RNA-Seq data processing rather simple and
and make changing to an alternative step (e.g., aligner, differential expression method, etc...)
trivial. This is also designed to make coexisting variants something a user does not
have to worry about (unless inclined to - the system is open).

.. note::

   The general principles to remember are limited to:

   - Steps require assets to run (and optionally parameters)

   - Assets are constituted of two groups: *source* and *target* elements

   - Parameters are optionally given

   The bundle of source and target assets, parameters, and a step represents a *task*.
   Explanations about step and assets follow.

   >>> source = SomeStep.Assets.Source( inputentity )
   >>> target = SomeStep.Assets.Target( outputentity )
   >>> assets = SomeStep.Assets(source, target)
   >>> step = SomeStep( pathtoexecutable )
   >>> step.run( assets )

Steps and assets
----------------

Concept
^^^^^^^

All steps in the process are connected through intermediate data files which we call assets.
Bioinformatics tools are almost always designed to operate on files, with the occasional pipes being used. 

A *step* can be represented as the step itself with input files (the *source* assets) and output files (the *target* assets).

.. note-stepassets-begin

.. graphviz::

   digraph ExampleAssets {
     "src1"->"step";
     "src1"[label="Some data", shape=invhouse, style=filled, fillcolor=lightgrey];
     "src2"->"step";
     "src2"[label="Other data", shape=invhouse, style=filled, fillcolor=lightgrey];
     "step"[label="Step"];
     "step"->"tgt1";
     "tgt1"[label="Product", shape=box];
   }

.. note-stepassets-end

This group of nodes (here 4 nodes: 2 source, 1 step, 1 target) can also be collapsed as a *step*, and the representation
of a workflow be the graph connecting this summary representation.

For example, the initial steps for aligner-based RNA-Seq can be represented as:

.. note-sa-example-high-begin

.. graphviz::

   digraph HighLevel {
     subgraph cluster_3 {
       label="High-level";
       style=filled;
       color=lightgrey;
       node [style=filled,color=white];
       "Index"->"Align"->"Count";
     }
   }

.. note-sa-example-high-end

The Step-and-Assets graph will then look like this:

.. note-sa-example-begin

.. graphviz::

   digraph StepAssets {
     subgraph cluster_0 {
       label="Index";
       color=lightgrey;
       "ReferenceGenome"->"Indexer";
       "ReferenceGenome"[label="Reference Genome", shape=invhouse, colorscheme=set36, style=filled, fillcolor=1];
       "Indexer"->"IndexedGenomeS";
       "Indexer"[colorscheme=set36, style=filled, fillcolor=2];
       "IndexedGenomeS"[label="Indexed Genome", shape=box, colorscheme=set36, style=filled, fillcolor=2, penwidth=2];
     }
     subgraph cluster_1 {
       label="Align";
       color=lightgrey;
       "IndexedGenomeT"->"Aligner";
       "IndexedGenomeT"[label="Indexed Genome", shape=box, colorscheme=set36, style=filled, fillcolor=2, penwidth=2];
       "Reads1" -> "Aligner";
       "Reads1"[label="Reads 1", shape=invhouse, colorscheme=set36, style=filled, fillcolor=5];
       "Reads2" -> "Aligner";
       "Reads2"[label="Reads 2", shape=invhouse, colorscheme=set36, style=filled, fillcolor=5];
       "Aligner" -> "AlignedReadsS";
       "Aligner"[colorscheme=set36, style=filled, fillcolor=2];
       "AlignedReadsS"[label="Aligned Reads", shape=box, colorscheme=set36, style=filled, fillcolor=2, penwidth=2];
     }
     "IndexedGenomeS"->"IndexedGenomeT"[style="dotted",arrowhead="none"];
     subgraph cluster_2 {
       label="Count";
       color=lightgrey;
       "AlignedReadsT"->"Counter";
       "AlignedReadsT"[label="Aligned Reads", shape=box, colorscheme=set36, style=filled, fillcolor=2, penwidth=2];
       "ReferenceAnnotation" -> "Counter";
       "ReferenceAnnotation"[label="Reference Annotation", shape=invhouse, colorscheme=set36, style=filled, fillcolor=5];
       "Counter" -> "DigitalGeneExpression";
       "Counter"[colorscheme=set36, style=filled, fillcolor=2];
       "DigitalGeneExpression"[label="Digital Gene Expression", shape=box, colorscheme=set36, style=filled, fillcolor=2, penwidth=2];
     }
     "AlignedReadsS"->"AlignedReadsT"[style="dotted",arrowhead="none"];
   }

.. note-sa-example-end

The nodes `Indexed Genome` and `Aligned Reads` are duplicated
for clarity with the grouping of nodes,
but it is the same entity saved on disk produced and used by the two steps.

When using the framework, the easiest way to think about it is to start from
a step (child class
of :class:`StepAbstract`), and look for the attribute :attr:`Assets`.
That attribute is a class representing the step's assets, and is itself 
containing 2 class attributes :attr:`Source`
and :attr:`Target` (themselves classes as well).

For example with the class modeling the aligner "STAR":

.. code-block:: python

   import railroadtracks.rnaseq import StarIndex, StarAlign

   # assets for the indexing of references (indexed for the alignment)
   StarIndex.Assets
   # assets for aligning reads against indexed references
   StarAlign.Assets

The assets are divided into two sets, the `Source` and the `Target`, for the files the step is using and the files
the step is producing respectively. :attr:`Source` and :attr:`Target` are themselves classes.

.. code-block:: python

   import railroadtracks.rnaseq import StarIndex, StarAlign

   # assets for the indexing of references (indexed for the alignment)
   fastaref = rnaseq.SavedFASTA('/path/to/referencegenome.fasta')
   indexedref = rnaseq.FilePattern('/path/to/indexprefix')
   StarIndex.Assets(StarIndex.Assets.Source(fastaref),
                    StarIndex.Assets.Target(indexedref))

The results are somewhat verbose, but if an IDE or an advanced interactive shell such as 
`ipython` is used, autocompletion will make writing such statements relatively painless, and mostly intuitive.
One can just start from the step considered (:class:`StarIndex` in the example above) and everything can be derived
from it.


Unspecified assets
^^^^^^^^^^^^^^^^^^

In a traditional approach, for example a sequence of commands in `bash`, the name
of files must be specified at each step.

We are proposing an alternative with which a full *recipe* can be written without having to take care of file names for
derived data (which can be a relative burden when considering to process the same data in many alternative ways).
Only the original data files such as the reference genome, or the experimental data
from the samples and sequencing, are specified and the other file names will be generated automatically.

Objects inheriting from :class:`AssetSet` are expected to have a method :meth:`AssetSet.createundefined` that
creates an undefined set of assets, and this can be used in recipes (see example below).

.. literalinclude:: ../src/test/test_core.py
   :language: python
   :start-after: # -- createundefined-begin
   :end-before: # -- createundefined-end



Notes about the design
^^^^^^^^^^^^^^^^^^^^^^

Data analysis is often a REPL activity, and we are keeping this in mind. The writing
of a recipes tries to provide:

- autocompletion-based discovery. Documentation is rarely read back-to-back
  (congratulations for making it this far though), and dynamic scientists often proceed by
  trial-and-error with software. The package is trying to provide benevolent support when doing so.

- fail early whenever possible (that is before long computations have already been performed)

- allow the writing of a full sequence of steps and run them unattended (tasks to be performed are stored, and executed when the user wants to)




Setup
-----


.. literalinclude:: ../src/test/test_recipe.py
   :language: python
   :start-after: # -- recipe-init-begin
   :end-before: # -- recipe-init-end


The package is also able to generate a small dataset based on a phage:

.. literalinclude:: ../src/test/test_recipe.py
   :language: python
   :start-after: # -- recipe-data-begin
   :end-before: # -- recipe-data-end


Simple recipe
-------------

.. literalinclude:: ../src/test/test_recipe.py
   :language: python
   :start-after: # -- recipesimple-test-begin
   :end-before: # -- recipesimple-test-end


The variable `wd` contains the directory with all intermediate data and the final results,
and `db_dn` is the database file.

.. note::

   If you clean up after yourself, but want to run the next recipe in this documentation,
   the setup step will have to be run again.


Loops, nested loops, and many variants
--------------------------------------

.. literalinclude:: ../src/test/test_recipe.py
   :language: python
   :start-after: # -- recipeloop-test-begin
   :end-before: # -- recipeloop-test-end



Docstrings
----------

.. automodule:: railroadtracks.easy
   :members:
   :special-members: __init__


Troubleshooting
---------------

The standard Python module :mod:`logging` is used for logging, and its documentation should be checked.
For example, a very simple way to activate logging at the level *DEBUG* on *stdout*:

.. code-block:: python

   import railroadtracks
   logger = railroadtracks.logger
   logger.setLevel(railroadtracks.logging.DEBUG)
   railroadtracks.basicConfig('/tmp/stdout')

