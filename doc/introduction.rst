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

Introduction
============

.. epigraph::

      It is true that the most ancient peoples, the first librarians, employed a language quite different from the one we speak today; it is true that a few miles to the right, our language devolves into dialect and that ninety floors above, it becomes incomprehensible.

      -- Jorge Luis Borges, The Library of Babel


The processing of Next-Generation Sequencing (NGS) data is generally achieved through a sequence of steps called a pipeline.
Changing steps in the pipeline, or changing parameters for steps, or adding input data while
wishing to compare the alternative outcomes requires tedious bookkeeping, tailored code for
the alternatives, and possibly computing several times the tasks common to several alternatives.

:mod:`railroadtracks` is a Python toolkit with the following capabilities:

**Track task dependencies**
    The toolkit is keeping persistently the files used by all tasks in a project, and is providing intuitive ways to navigate 
    complex dependency graphs.
**Decouple declaration from execution**
    `railroadtracks` allows the creation of all interconnected computational tasks required
    for a project before any computation is performed. This, added to the model aspect, allows
    consistency checks early instead of discovering issues after a number of computation tasks have
    already been performed. This is also permitting the use different scheduling systems within a project.
    The toolkit is providing task execution engines using :mod:`multiprocessing` and SGE's qsub,
    and custom execution layers can be added.
**Implicit naming of result files**
   The toolkit is making a clear separation between the meta-data for a result file (how was a result file obtained) from its file name. User can
   write complete pipelines, while only specifying the file names for the input data at the root of the dependency graph. 
**Incremental writing of analysis pipelines**
   The persistent tracking of task dependencies permits the incremental writing of data processing pipelines without explicit
   bookkeeping nor unnecessarily repeated computations. Pipelines can be modified, with new input files added, tools used changed,
   while only new required computations are performed.
**Standardization into models**
   By providing a system to wrap tools into model classes, :mod:`railroadtracks` allows a standardization of these tools that allows both an easy
   connection between steps and the ability to swap between similar tools. This and the implicit naming of result files makes the
   combination of tools or parameters trivial.
**REPL-aware**
   :mod:`railroadtracks` was designed with interactive programming (REPL) in mind, and ipython and the ipython notebook as a showcase environment.
   High-level rendering of objects is provided (text, HTML) as well as visualization (provenance and destination graphs) are provided.
   The autocompletion of Python namespaces is also part of our design. 
**Models DNA and RNA sequencing (batteries included)**
   While it is possible to extend :mod:`railroadtracks` with custom models, a number of modeled steps for
   generating synthetic reads, aligning, quantifying, testing for differential expression, are already included.


