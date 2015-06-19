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

Installation
============

Requirements
------------

This was developed for Python 2.7 and Python 3.4

.. note:: 

   While a keen attention has been paid to compatibility with Python 3,
   the use of metaclasses in this packages and the incompatibility of syntax with them 
   between Python 2 and 3 prevents this package from currently working with Python 3.
   This is the major blocker for Python 3-compatibility and we plan on using :mod:`six`
   to solve this.

A number of python packages are required or recommended (all of which are available on `Pypi <https://pypi.python.org/pypi>`):

- `jinja2 <http://jinja.pocoo.org/>`
- `networkx <https://networkx.github.io/>`
- `pygraphviz <https://pygraphviz.github.io/>`
- `enum34 <https://pypi.python.org/pypi/enum34>`
- `ngs_plumbing <http://pythonhosted.org/ngs_plumbing/>`

The model steps included require a number of third-party tools in order to have all steps working.
The tools themselves are not included and a missing tool does not prevent this package from working, but the related functionality (running that tool) will not be working.

At the time of writing, the tools are:

- `bedtools2 <http://bedtools.readthedocs.org/en/latest/index.html>`
- `bwa <http://bio-bwa.sourceforge.net/>`
- `bowtie <http://bowtie-bio.sourceforge.net/>`
- `bowtie2 <http://bowtie-bio.sourceforge.net/bowtie2>`
- `GSNAP <http://research-pub.gene.com/gmap/>`
- `HTSeq <http://www-huber.embl.de/users/anders/HTSeq>` (for the executable `htseq-count`)
- `sailfish <http://www.cs.cmu.edu/~ckingsf/software/sailfish/>`
- `samtools <http://samtools.sourceforge.net>`
- `STAR <https://github.com/alexdobin/STAR>`
- `tophat <http://ccb.jhu.edu/software/tophat/index.shtml>`
- `tophat2 <http://ccb.jhu.edu/software/tophat/index.shtml>`
- `R <http://www.r-project.org>`, with libraries:
    - `rjson <http://cran.r-project.org/package=rjson>` (required for any interaction with R)
  and bioconductor libraries
    - `Rsubread <http://www.bioconductor.org/packages/release/bioc/html/Rsubread.html>` (for `featureCount`)
    - `DESeq <http://www.bioconductor.org/packages/release/bioc/html/DESeq.html>`
    - `DESeq2 <http://www.bioconductor.org/packages/release/bioc/html/DESeq2.html>`
    - `limma <http://www.bioconductor.org/packages/release/bioc/html/limma.html>` (for the method `voom`)
    - `edgeR <http://www.bioconductor.org/packages/release/bioc/html/edgeR.html>`
    - `EBSeq <http://www.bioconductor.org/packages/release/bioc/html/EBSeq.html>`

.. note::

   Running the tests with `-v` (see below) will tell which tools are missing to have all functionality (see below for instructions about running the tests).


The package is including a small genome
(Enterobacteria phage MS2 isolate ST4 retrieved from the NCBI site: http://www.ncbi.nlm.nih.gov/nuccore/EF204940.1 - 
bibliographic reference `Friedman SD et al., J Virol. 2009 Nov;83(21):11233-43. doi: 10.1128/JVI.01308-09. <http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2772794/>`) as FASTA, as well as
GTF/GFF derived from its genome annotation, for the purpose of running tests and examples.

Installation
------------

.. note-installation-slide-begin

This is a regular Python package. It can be installed with

.. code-block:: bash

   python setup.py install


Running the unit-tests can be done after installation

.. code-block:: bash

   python -m unittest discover railroadtracks

.. note-installation-slide-end


.. note::

   While the unit tests use a very small dataset included in the package (a viral phage genome),
   running the test can take a bit of disk space (for instance, at the time of writing
   STAR indexes take a rather large amount of space for an otherwise rather small genome).
