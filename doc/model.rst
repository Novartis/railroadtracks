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

Model
=====

Overview
--------


The steps in a RNA-Seq sequence of operations are captured in a model.
The purpose is to formalize the steps enough to simplify the use of various tools for each step,
yet permit enough flexibility to allow an easy integration of additional tools and
approaches.


A canonical graph of steps for RNA-Seq is shown below. 

.. note::

   There is no enforcement that a sequence of steps adheres strictly to this.
   The model is also defining the notion of `activities`, and it is possible
   that one step performs several activities (this will be the case for
   an in-house monolithic pipeline, for example). This feature allows
   us to integrate such steps into the comparison.

The model is directly used to provide an unified execution scheme (see :ref:`unified-execution-label`),
and constitutes the basis for writing "recipes" (see :ref:`recipes-label`).


.. note-canonical-model-begin

.. graphviz::

     digraph RNASeq {
     "Genomes"->"References";
     "Genomes"[label="repository of genomes", shape=invhouse, colorscheme=set36, style=filled, fillcolor=1];
     "References" -> "IndexedReference" -> "Alignment";
     "References"[shape=invhouse, colorscheme=set36, style=filled, fillcolor=1];
     "IndexedReference"[label="Indexed References", colorscheme=set36, style=filled, fillcolor=2];
     "SequencingReads" -> "Alignment";
     "SequencingReads"[label="Sequencing Reads", shape=invhouse, colorscheme=set36, style=filled, fillcolor=5];
     "Alignment" -> "Counts";
     "TranscriptAnnotation" -> "Alignment";
     "Alignment"[colorscheme=set36, style=filled, fillcolor=2];
     "Genomes"->"TranscriptAnnotation";
     "TranscriptAnnotation" -> "Counts";
     "TranscriptAnnotation"[label="Transcript Annotations", shape=invhouse, colorscheme=set36, style=filled, fillcolor=1];
     "Counts" -> "NormalizedCounts";
     "Counts"[colorscheme=set36, style=filled, fillcolor=3];
     "NormalizedCounts" -> "AbundanceEstimates";
     "NormalizedCounts"[label="Normalized Counts", colorscheme=set36, style=filled, fillcolor=3];
     "SampleInformation" -> "DifferentialExpression";
     "SampleInformation"[label="Sample Information", shape=invhouse, colorscheme=set36, style=filled, fillcolor=5];
     "Sample"->"SampleInformation";
     "Sample"->"SequencingReads";
     "Sample"[shape=invhouse, colorscheme=set36, style=filled, fillcolor=5];
     "AbundanceEstimates"->"DifferentialExpression";
     "TranscriptAnnotation" -> "AbundanceEstimates";
     "AbundanceEstimates"[label="Abundance Estimates", colorscheme=set36, style=filled, fillcolor=4];
     "DifferentialExpression"[label="Differential Expression", shape=invhouse, colorscheme=set36, style=filled, fillcolor=6];
     }

.. note-canonical-model-end


.. _sec-activities-label:

Activities
----------


Steps can perform one or several activities.
The list of possible activities is in :class:`ACTIVITY`

Inheritance diagram
-------------------

.. inheritance-diagram:: railroadtracks.rnaseq
   :parts: 1


Docstrings
----------

.. automodule:: railroadtracks.rnaseq
   :members:
   :special-members: __init__


