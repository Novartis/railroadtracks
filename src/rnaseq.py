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
Simple model for the expression analysis using RNA-Seq
"""
import os, sys, shutil
import re
import warnings
import csv
import zlib
import collections
import subprocess
import tempfile
import logging
logger = logging.getLogger(__name__)
import argparse
from . import environment
from railroadtracks.model.diffexp import RSOURCES_DIR
from . import core
import gzip
from railroadtracks.model import files
from railroadtracks.model.files import (File,
                                        FASTAFile,
                                        FASTQFile,
                                        GFFFile,
                                        FASTQPossiblyGzipCompressed, 
                                        BAMFile,
                                        SAMFile,
                                        BAMFileSortedByID,
                                        CSVFile,
                                        TSVFile,
                                        FilePattern,
                                        CSVFileSequence,
                                        SAMFileSequence,
                                        SamtoolsSamToBam,
                                        SamtoolsBamToSam,
                                        SamtoolsSorterByID,
                                        BEDFile,
                                        BedtoolsBamToFastq,
                                        BedtoolsBamToFastqPE)
from railroadtracks.model import simulate
from railroadtracks.model.simulate import (FluxsimulatorExpression,
                                           FluxsimulatorDerivedExpression,
                                           FluxsimulatorSequencing,
                                           FluxSimulatorParameters,
                                           FluxSimulatorParametersAndPro)
from railroadtracks.model import aligners
from railroadtracks.model.aligners import (BWAIndex,
                                           BWAIndexFiles,
                                           BWA,
                                           Bowtie,
                                           SavedBowtieIndex,
                                           Bowtie2,
                                           SavedBowtie2Index,
                                           BowtieBuild,
                                           Bowtie2Build,
                                           StarIndex,
                                           StarAlign,
                                           GsnapIndex,
                                           GsnapAlign,
                                           SailfishIndex,
                                           SalmonIndex,
                                           TopHat,
                                           TopHat2,
                                           SamtoolsExtractUnaligned)
from railroadtracks.model import diffexp
from railroadtracks.model.diffexp import (EdgeR,
                                          DESeq,
                                          DESeq2,
                                          LimmaVoom,
                                          EBSeq)

from railroadtracks.model import quantify
from railroadtracks.model.quantify import (SailfishQuant,
                                           SalmonAlignQuant,
                                           FeatureCount,
                                           HTSeqCount)

from railroadtracks.model import misc
from railroadtracks.model.misc import Anyscript

from railroadtracks.unifex import UnifexError, _cmdfromuei
from abc import ABCMeta, abstractproperty
            

class ACTIVITY(simulate.ACTIVITY, 
               aligners.ACTIVITY, 
               quantify.ACTIVITY, 
               diffexp.ACTIVITY,
               files.ACTIVITY,
               misc.ACTIVITY):
    """ Activities that can be performed by the different steps modeled.
    (note: steps can combine several activities - the most obvious example is
    existing monolithic pipelines) """
    __order__ = 'SIMULATE INDEX ALIGN SORT QUANTIFY NORMALIZE DIFFEXP UTILITY MISC'
    NORMALIZE = 'Normalize'



class AssetsCRCHeadTail(core.AssetsStep):
    """
    Assets for :class:`CRCHeadTail`

    """
    Source = core.assetfactory('Source', [core.AssetAttr('file', File, '')])
    Target = core.assetfactory('Target', [core.AssetAttr('crc', CSVFile, '')])

class CRCHeadTail(core.StepAbstract):
    """
    Compute a CRC32-based checksum on the beginning (head) and end (tail) of a file
    as a cheap way to check whether 2 files contain identical data.
    Might be useful with large files, however its main purpose is testing.
    """
    _name = 'crc-file'
    _default_execpath = None
    Assets = AssetsCRCHeadTail
    activities = (ACTIVITY.UTILITY, )

    def __init__(self, executable=None):
        if executable is None:
            executable = type(self)._default_execpath
        self._execpath = executable
        # no third-party executable involved
        if self._execpath is not None:
            raise ValueError('No third-party executable involved. The executable should be None')


    @property
    def version(self):
        from railroadtracks import __version__
        return __version__

    def run(self, assets, parameters = tuple()):
        # FIXME: should these checks be moved to the Assets class ?
        out_counts = assets.target.crc
        out_fh = open(out_counts.name, 'w')
        csv_w = csv.writer(out_fh)
        csv_w.writerow(['name', 'crc'])
        for name in assets.source.file:
            with open(name, 'rb') as in_fh:
                # begining of the file
                data = in_fh.read(10000)
                crc = zlib.crc32(data)
                # end of the file
                in_fh.seek(0, 2)
                n = in_fh.tell()
                if n < 10000:
                    in_fh.seek(0)
                else:
                    in_fh.seek(-10000, 2)
                data = in_fh.read(10000)
                crc = zlib.crc32(data, crc)                
            csv_w.writerow((name, crc))
        out_fh.flush()
        out_fh.close()
        cmd = None
        returncode = 0
        return (cmd, returncode)


# patch ngs_plumbing with class Gzip-friendly
try:
    import ngs_plumbing.fastq
    hasngsp = True
except:
    hasngsp = False
if hasngsp:
    import io
    #FIXME: This should ideally be fixed upstream in ngs_plumbing
    class GzipFastqFilePair(ngs_plumbing.fastq.FastqFilePair):
        """ Pair of FASTQ files
        The default iterator will go through the paired entries in the file
        (not the rows), assuming that they are in the same order.
        No check that this is the case (using IDs) is performed.
        """

        def __init__(self, filename1, filename2, **kwargs):
            Cls = type(self)
            if ngs_plumbing.fastq._python2:
                fh1 = gzip.open(filename=filename1, **kwargs)
                fh2 = gzip.open(filename=filename2, **kwargs)
            else:
                stream = io.FileIO(filename1)
                fh1 = gzip.open(stream, **kwargs)
                stream = io.FileIO(filename2)
                fh2 = gzip.open(stream, **kwargs)
            self._fh1 = fh1
            self._fh2 = fh2

        def __iter__(self):
            return ngs_plumbing.fastq.iter_entry_pairs(self._fh1, self._fh2)

        def iter_seqqual_pairs(self):
            """ Iterate over the pairs of (sequence+quality) in 
            the file. """
            return ngs_plumbing.fastq.iter_seqqual_pairs(self._fh1, self._fh2)

        def readcount(self):
            pos = self.tell()
            for i, entry in enumerate(self):
                pass
            return i+1
    




        

class AssetsColumnMerger(core.AssetsStep):
    Source = core.assetfactory('Source', [core.AssetAttr('counts', CSVFileSequence, '')])
    Target = core.assetfactory('Target', [core.AssetAttr('counts', CSVFile, '')])

class ColumnMerger(core.StepAbstract):
    Assets = AssetsColumnMerger
    _name = 'column-merge'
    _default_execpath = None
    activities = (ACTIVITY.UTILITY, )
    #activity = property(lambda x: x.activities, None, None)

    def __init__(self,executable=None):
        """ First parameter is the column used as ID. If None, no column used as ID.
        Second parameter is the column to merge."""
        if executable is None:
            executable = type(self)._default_execpath
        self._execpath = executable # no third-party executable involved
        if  self._execpath is not None:
            raise ValueError('The executable cannot be set. This is a fully internal step.')

    @property
    def version(self):
        from railroadtracks import __version__
        return __version__

    def run(self, assets, parameters = ('0','1')):
        # FIXME: should these checks be moved to the Assets class ?
        out_counts = assets.target.counts
        l = len(tuple(out_counts))
        if l > 1:
            raise ValueError('Merging columns can only have one target file.')
        else:
            import resource
            max_fh_soft, max_fh_hard = resource.getrlimit(resource.RLIMIT_NOFILE)
            if l > max_fh_hard:
                raise ValueError('On this system the maximum number of columns is'
                                 ' %i (and you are trying to merge %i columns)'
                                 ' [note: opened file handles should be taken out of the'
                                  ' maximum number]' % (max_fh_hard-1, l))
        in_counts = assets.source.counts
        # open all
        it_in_counts = in_counts.iteritems()
        cls, value = next(it_in_counts)
        head = cls(value)
        head_iter = head.iter_entry()
        tail = tuple(cls(value) for cls, value in it_in_counts)
        tail_iter = tuple(x.iter_entry() for x in tail)
        # merge
        if len(parameters) != 2:
            raise ValueError('There should be 2 parameters: the index of the columns with IDs, the index of the column with values.')
        try:
            idcol_i = int(parameters[0])
        except ValueError:
            if parameters[0] == 'None':
                idcol_i = None
            else:
                ValueError('The first parameter should be string representation of an integer or "None".')
        valcol_i = int(parameters[1])
        with open(out_counts.name, 'w') as out_fh:
            csv_w = csv.writer(out_fh)                    
            if idcol_i is None:
                for row in head_iter:
                    counts = [row[valcol_i]]
                    for col_i, h in enumerate(tail_iter):
                        row = next(h)
                        counts.append(row[valcol_i].rstrip())
                    csv_w.writerow(counts)
            else:
                for row in head_iter:
                    identifier_head = row[idcol_i]
                    counts = [row[valcol_i].rstrip()]
                    for col_i, h in enumerate(tail_iter):
                        row = next(h)
                        identifier = row[idcol_i]
                        if identifier == identifier_head:
                            counts.append(row[valcol_i].rstrip())
                        else:
                            raise ValueError('Mismatching IDs (first columns has "%s"'
                                             ' while column %i has "%s")' % (identifier_head,
                                                                             col_i+1, identifier))
                    csv_w.writerow([identifier_head, ] + counts)
        
        cmd = None
        returncode = 0
        return (cmd, returncode)

class AssetsNormalizer(core.AssetsStep):
    Source = core.assetfactory('Source', [core.AssetAttr('count', CSVFile, '')])
    Target = core.assetfactory('Target', [core.AssetAttr('normalizedcounts', CSVFile, '')])

class Normalizer(core.StepAbstract):
    Assets = AssetsNormalizer


def _picard_getversion(execpath):
    cmd = [execpath, '--version']
    try:
        logger.debug(subprocess.list2cmdline(cmd))
        res = subprocess.check_output(cmd)
    except OSError as ose:
        raise UnifexError("""Command: %s
        %s""" % (' '.join(cmd), ose))

    m = re.match('^([0-9].+)$', res)
    if m is None:
        raise RuntimeError('Could not find the version number.')
    return m.groups()[0]



class PicardCollectAlignmentSummaryMetrics(core.StepAbstract):
    activities = (ACTIVITY.UTILITY, )
    _name = 'picard-collect-alignment-summary-metrics'
    _default_execpath = 'picard-collect-alignment-summary-metrics.sh'

    class Assets(core.AssetsStep):
        Source = core.assetfactory('Source', [core.AssetAttr('alignedreads', BAMFileSortedByID, ''),
                                              core.AssetAttr('reference', FASTAFile, '')])
        Target = core.assetfactory('Target', [core.AssetAttr('metrics', FilePattern, '')])

    parser = argparse.ArgumentParser(_name)
    parser.add_argument('ASSUME_SORTED',
                        default = 'false',
                        choices = ('true', 'false'))

    def __init__(self, execpath=None):
        """        
        :param executable: path to the executable
        :type param: :class:`str`
        """
        if execpath is None:
            execpath = type(self)._default_execpath
        self._execpath = execpath
        self._version = None
    
    @property 
    def version(self):
        if self._version is None:
            self._version = _picard_getversion(self.execpath)
        return self._version

    def run(self, assets, parameters = tuple()):
        options, unknown = self.parser.parse_known_args()
        # FIXME: check that input is sorted BAM. If not, forbid the use of ASSUME_SORTED ?

        # build command line
        cmd = [self._execpath, ]
        cmd.extend(parameters)
        cmd.extend(['INPUT_FILE=%s', assets.source.annotationfile.name,
                    'REFERENCE_SEQUENCE=%s', assets.source.reference.name,
                    'OUTPUT=%s', assets.target.metrics.name])

        logger.debug(subprocess.list2cmdline(cmd))
        with open(os.devnull, "w") as fnull:
            returncode = subprocess.check_call(cmd)

        return (cmd, returncode)

