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

import csv
import os
import re
import tempfile
import logging
import warnings
import shutil
import abc
import argparse
logger = logging.getLogger(__name__)
import subprocess
from railroadtracks import core, environment
from railroadtracks.core import SavedEntityAbstract, File, FileSequence

from six import with_metaclass

import sys
if sys.version_info[0] < 3:
    def astring(name):
        return isinstance(name, str) or isinstance(name, unicode)
else:
    from functools import reduce
    def astring(name):
        return isinstance(name, str) or isinstance(name, bytes)


class ACTIVITY(core.Enum):
    """ Activities that can be performed by the different steps modeled.
    (note: steps can combine several activities - the most obvious example is
    existing monolithic pipelines) """
    __order__ = 'CONVERTFORMAT'
    CONVERTFORMAT = 'Convert file formats'
    FILTERREADS = 'Filter reads'
    SORT = 'Sort'

class FASTAFile(File):
    """ FASTA file """
    _extension = ('.fa', '.fna', '.fasta', '.FASTA')

class GFFFile(File):
    """ GFF file """
    _extension = ('.gff', '.gtf')

class FASTQFile(File):
    """ FASTQ file """
    _extension = ('.fq', '.fastq')

class BAMFile(File):
    """ BAM file """
    _extension = ('.bam', )

class SAMFile(File):
    """ SAM file """
    _extension = ('.sam', )

class BEDFile(File):
    """ BED file """
    _extension = ('.bed', )

class FASTQPossiblyGzipCompressed(FASTQFile):
    _extension = tuple(y for ext in FASTQFile._extension for y in (ext, ext+'.gz'))
    @property
    def iscompressed(self):
        return self.name.endswith('.gz')

class SAMFileSortedByID(SAMFile):
    """ SAM file in which entries are sorted by read ID
    (this is required by tools such as HTSeq.). """
    pass

class BAMFileSortedByID(BAMFile):
    """ BAM file in which entries are sorted by read ID
    (this is required by tools such as HTSeq.). """
    pass

class CSVFile(File):
    """ CSV file """
    _extension = ('.csv', )
    def iter_entry(self):
        fh = open(self.name)
        csv_r = csv.reader(fh)
        for row in csv_r:
            yield row
        fh.close()

class TSVFile(File):
    """ Tab-separated files """
    _extension = ('.tsv', )
    def iter_entry(self):
        fh = open(self.name)
        csv_r = csv.reader(delimiter='\t')
        for row in csv_r:
            yield row
        fh.close()

# FIXME: should be Prefix, not Pattern
class FilePattern(SavedEntityAbstract):
    """
    Basename for a file.

    Some of the bioinformatics tools work with "basenames" and will produce or consume an arbitrary number
    of files.

    This object is able to represent a number of files having a given pattern.
    """

    def __init__(self, name):
        if name is not None:
            self._defined = True
            if not astring(name):
                if len(name) == 1:
                    name = next(iter(name))
                else:
                    raise ValueError('"name" should be a string.')
        self._name = name

    def _getname(self):
        assert self._defined, "The %s object is not defined." % type(self).__name__
        return self._name

    def _setname(self, name):
        assert not self._defined, "The %s object is already defined." % type(self).__name__
        self._name = name
        self._defined = True

    name = property(_getname, _setname, None,
                      '"Name" for a family of files. When not defined, this can be None and accessing it will raise an AssertionError (it can only be set then).')

    def __len__(self):
        return 1

    def __iter__(self):
        """ All files with ".name" as a prefix,
        and an extension matching the one in the list if it is defined, are returned """
        for x in self.iteritems():
            yield x[1]
        
    def iteritems(self):
        yield (type(self), self.name)

    def iterlistfiles(self):
        dn = os.path.dirname(self.name)
        if dn=='':
            dn = '.'
        bn = os.path.basename(self.name)
        for fn in os.listdir(dn):
            # if file extensions are defined, filter on extensions
            if (self._extension is not None) and (not any(fn.endswith(x) for x in self._extension)):
                continue            
            if fn.startswith(bn):
                yield (type(self), fn)
        

    def lastmodified(self):
        assert self._defined, "The %s object is not defined" % type(self).__name__
        it = iter(self)
        try:
            lasttime = next(it)
        except StopIteration:
            return None
        for fn in it:
            t = os.path.getmtime(os.path.join(self._dirname, fn))
            if t > lasttime:
                lasttime = t
        return lasttime


class CSVFileSequence(FileSequence):
    """ Sequence of CSV files """
    _type = CSVFile

class SAMFileSequence(FileSequence):
    """ Sequence of SAM files """
    _type = SAMFile


class SamtoolsFilter(core.StepAbstract):
    """
    Filter reads in a BAM file
    """
    _name = 'bam-filter'
    _default_execpath = 'samtools'

    class Assets(core.AssetsStep):
        """
        Assets for :class:`SamtoolsFilter`
        
        """
        Source = core.assetfactory('Source', [core.AssetAttr('bamfile', BAMFile, '')])
        Target = core.assetfactory('Target', [core.AssetAttr('bamfile', BAMFile, '')])

    activities = (ACTIVITY.FILTERREADS, )

    parser = argparse.ArgumentParser(_name)
    FILTERVALUE_INFO = """
Flag        Chr     Description
0x0001      p       the read is paired in sequencing
0x0002      P       the read is mapped in a proper pair
0x0004      u       the query sequence itself is unmapped
0x0008      U       the mate is unmapped
0x0010      r       strand of the query (1 for reverse)
0x0020      R       strand of the mate
0x0040      1       the read is the first read in a pair
0x0080      2       the read is the second read in a pair
0x0100      s       the alignment is not primary
0x0200      f       the read fails platform/vendor quality checks
0x0400      d       the read is either a PCR or an optical duplicate
"""
    parser.add_argument('-f', '--filter-include',
                        dest = 'filter_include',
                        nargs = '+',
                        choices = ('0x0001','0x0002','0x0004','0x0008',
                                   '0x0010','0x0020','0x0040','0x0080',
                                   '0x0100','0x0200','0x0400','0x0800'),
                        help = """ Any of the following filter values can be specified.
Only return the reads satisfying all the filters.

                        %s""" % FILTERVALUE_INFO)
    parser.add_argument('-F', '--filter-exclude',
                        dest = 'filter_exclude',
                        nargs = '+',
                        choices = ('0x0001','0x0002','0x0004','0x0008',
                                   '0x0010','0x0020','0x0040','0x0080',
                                   '0x0100','0x0200','0x0400','0x0800'),
                        help = """ Any of the following filter values can be specified.
Only return the reads not satisfying any the filters.
                        %s""" % FILTERVALUE_INFO)


    def __init__(self, executable=None):
        if executable is None:
            executable = type(self)._default_execpath
        self._execpath = executable
        self._version = None

    @property
    def version(self):
        if self._version is None:
            self._version = samtools_getversion(self._execpath)
        return self._version
        
    def run(self, assets, parameters = tuple()):
        args = self.parser.parse_args(parameters)
        # FIXME: move checks out so can they can be made before calling 'run'
        if (args.filter_exclude is None) and (args.filter_include is None):
            raise ValueError('No filter specified.')

        cmd = [self._execpath, 'view', '-o', assets.target.bamfile.name]
        if (args.filter_exclude is not None):
            filter_exclude = reduce(lambda x,y: x^y, (x for x in args.filter_exclude))
            cmd.extend(('-F', filter_exclude))
        if (args.filter_include is not None):
            filter_include = reduce(lambda x,y: x^y, (x for x in args.filter_include))
            cmd.extend(('-f', filter_include))
            
        cmd.extend(('-b', assets.source.bamfile.name))

        with open(os.devnull, 'w') as fnull: 
            logging.debug(cmd)
            returncode = subprocess.check_call(cmd,
                                               stdout = fnull,
                                               stderr = fnull)
        if not os.path.exists(assets.source.bamfile.name):
            # target is missing. Suspected infamous 'bam.bam' issue
            bam_bam = assets.source.bamfile.name + '.bam'
            if os.path.exists(bam_bam):
                # '.bam.bam' issue
                warnings.warn("'.bam.bam' issue detected. Moving the product to intended target %s." % assets.source.samfile.name)
                shutil.move(bam_bam, assets.source.bamfile.name)
            else:
                raise Exception('The target %s is mysteriously missing.' % assets.source.samfile.name)
        return (cmd, returncode)


class SamtoolsSamToBam(core.StepAbstract):
    """
    Convert a SAM file into a BAM file
    """
    _name = 'sam-to-bam'
    _default_execpath = 'samtools'

    class Assets(core.AssetsStep):
        """
        Assets for :class:`SamtoolsSamToBam`
        
        """
        Source = core.assetfactory('Source', [core.AssetAttr('samfile', SAMFile, '')])
        Target = core.assetfactory('Target', [core.AssetAttr('bamfile', BAMFile, '')])

    activities = (ACTIVITY.CONVERTFORMAT, )

    def __init__(self, executable=None):
        if executable is None:
            executable = type(self)._default_execpath
        self._execpath = executable
        self._version = None

    @property
    def version(self):
        if self._version is None:
            self._version = samtools_getversion(self._execpath)
        return self._version
        
    def run(self, assets, parameters = tuple()):
        cmd = [self._execpath, 'view', '-o', assets.target.bamfile.name, '-Sb', assets.source.samfile.name, ]
        with open(os.devnull, 'w') as fnull: 
            logging.debug(cmd)
            returncode = subprocess.check_call(cmd,
                                               stdout = fnull,
                                               stderr = fnull)
        if not os.path.exists(assets.source.samfile.name):
            # target is missing. Suspected infamous 'bam.bam' issue
            bam_bam = assets.source.samfile.name + '.bam'
            if os.path.exists(bam_bam):
                # '.bam.bam' issue
                warnings.warn("'.bam.bam' issue detected. Moving the product to intended target %s." % assets.source.samfile.name)
                shutil.move(bam_bam, assets.source.samfile.name)
            else:
                raise Exception('The target %s is mysteriously missing.' % assets.source.samfile.name)
        return (cmd, returncode)


class SamtoolsBamToSam(core.StepAbstract):
    """
    Convert a BAM file into a SAM file
    """
    _name = 'bam-to-sam'
    _default_execpath = 'samtools'

    class Assets(core.AssetsStep):
        """
        Assets for :class:`SamtoolsBamToSam`
        
        """
        Source = core.assetfactory('Source', [core.AssetAttr('bamfile', BAMFile, '')])
        Target = core.assetfactory('Target', [core.AssetAttr('samfile', SAMFile, '')])

    activities = (ACTIVITY.CONVERTFORMAT, )

    def __init__(self, executable=None):
        if executable is None:
            executable = type(self)._default_execpath
        self._execpath = executable
        self._version = None

    @property
    def version(self):
        if self._version is None:
            self._version = samtools_getversion(self._execpath)
        return self._version
        
    def run(self, assets, parameters = tuple()):
        cmd = [self._execpath, 'view', '-h', '-o', assets.target.samfile.name, assets.source.bamfile.name, ]
        with open(os.devnull, 'w') as fnull: 
            logging.debug(cmd)
            returncode = subprocess.check_call(cmd,
                                               stdout = fnull,
                                               stderr = fnull)
        return (cmd, returncode)


def samtools_getversion(execpath):
    """ Return the version of 'samtools'. """

    cmd = [execpath,]
    logging.debug(cmd)
    m = None
    with core.Popen(cmd,
                    stderr=subprocess.PIPE) as proc:

        for row in proc.stderr:
            m = re.match(b'^Version: ([^ \n]+).*$', row)
            if m is not None:
                break
    if m is None:
        raise RuntimeError('Could not find the version number.')
    version = m.groups()[0]
    return version


def ensure_bam(filename):
    """ Ensure that a file is in the BAM format.
    If 'filename' is already a BAM file, return a filehandle, otherwise
    create a tempfile and return a filehandle to it.
    """
    # cheap test based on file extension
    if filename.endswith('.bam'):
        return open(filename, 'rb')
    elif filename.endswith('.sam'):
        samtools = environment.Executable('samtools')
        newfile = tempfile.NamedTemporaryFile(suffix='.bam')
        cmd = [samtools.path, 'view', '-o', newfile.name, '-Sb', filename, ]
        with open(os.devnull, 'w') as fnull: 
            logger.debug(subprocess.list2cmdline(cmd))
            returncode = subprocess.check_call(cmd,
                                               stdout = fnull,
                                               stderr = fnull)
        return newfile
    else:
        raise ValueError('The filename "%s" is not .sam neither .bam' % filename)

def ensure_sam(filename):
    """ Ensure that a file is in the SAM format.
    If 'filename' is already a SAM file, return a filehandle, otherwise
    create a tempfile and return a filehandle to it.
    """
    # cheap test based on file extension
    if filename.endswith('.bam'):
        samtools = environment.Executable('samtools')
        newfile = tempfile.NamedTemporaryFile(suffix='.sam')
        cmd = [samtools.path, 'view', '-o', newfile.name, '-h', filename, ]
        with open(os.devnull, 'w') as fnull: 
            logger.debug(subprocess.list2cmdline(cmd))
            returncode = subprocess.check_call(cmd,
                                               stderr = fnull,
                                               stdout = fnull)
        return newfile
    elif filename.endswith('.sam'):
        return open(filename, 'r')
    else:
        raise ValueError('The filename "%s" is not .sam neither .bam' % filename)


class AssetsSorter(core.AssetsStep):
    Source = core.assetfactory('Source', [core.AssetAttr('alignedreads',BAMFile, '')])
    Target = core.assetfactory('Target', [core.AssetAttr('sortedbam', BAMFile, '')])

class SorterAbstract(with_metaclass(abc.ABCMeta, core.StepAbstract)):
    """

    A sorting step.

    """
    activities = (ACTIVITY.SORT, )
    Assets = AssetsSorter

class SamtoolsSorterByID(SorterAbstract):
    """
    """

    _name = 'samtools-sortbyid'
    _default_execpath = 'samtools'

    def __init__(self, executable=None):
        if executable is None:
            executable = type(self)._default_execpath
        self._execpath = executable
        self._version = None

    @property
    def version(self):
        if self._version is None:
            self._version = samtools_getversion(self._execpath)
        return self._version

    def run(self, assets, parameters = tuple()):        
        # build command line
        cmd = ['%s' % self._execpath]
        cmd.append('sort')
        cmd.append('-n')
        cmd.append('-f') # otherwise 'samtools' is interpreting the output as a prefix rather than a complete file name
        cmd.extend(parameters)
        cmd.append('%s' % assets.source.alignedreads.name)
        cmd.append('%s' % assets.target.sortedbam.name) 
        with open(os.devnull, "w") as fnull:
            logger.debug(subprocess.list2cmdline(cmd))
            returncode = subprocess.check_call(cmd, 
                                               stdout = fnull,
                                               stderr = fnull)
        return (cmd, returncode)



class BedtoolsBamToFastq(core.StepAbstract):
    """
    Convert a BAM file into a FASTQ file
    """
    _name = 'BAM-to-FASTQ'
    _default_execpath = 'bedtools'

    class Assets(core.AssetsStep):
        """
        Assets for :class:`BedtoolsBamToFastq`
        
        """
        Source = core.assetfactory('Source', [core.AssetAttr('bamfile', BAMFile, '')])
        Target = core.assetfactory('Target', [core.AssetAttr('fastqfile', FASTQFile, '')])

    activities = (ACTIVITY.CONVERTFORMAT, )

    def __init__(self, executable=None):
        if executable is None:
            executable = type(self)._default_execpath
        self._execpath = executable
        self._version = None

    @property
    def version(self):
        if self._version is None:
            cmd = ('bedtools', '--version')
            output = subprocess.check_output(cmd)
            m = re.match("^bedtools v(.+)", output)
            self._version = m.groups()[0]
        return self._version
        
    def run(self, assets, parameters = tuple()):
        cmd = [self._execpath, 'bamtofastq',
               '-i', assets.source.bamfile.name, 
               '-fq', assets.target.fastqfile.name, ]
        with open(os.devnull, 'w') as fnull: 
            logging.debug(cmd)
            returncode = subprocess.check_call(cmd,
                                               stdout = fnull,
                                               stderr = fnull)
        return (cmd, returncode)


class BedtoolsBamToFastqPE(BedtoolsBamToFastq):


    _name = 'BAM-to-FASTQ-PE'

    class Assets(BedtoolsBamToFastq.Assets):
        Target = core.assetfactory('Target', [core.AssetAttr('fastqfile', FASTQFile, ''),
                                              core.AssetAttr('fastqfile2', FASTQFile, '')])

    activities = (ACTIVITY.CONVERTFORMAT, )
        
    def run(self, assets, parameters = tuple()):
        cmd = [self._execpath, 'bamtofastq',
               '-i', assets.source.bamfile.name, 
               '-fq', assets.target.fastqfile.name,
               '-fq2', assets.target.fastqfile2.name ]
        with open(os.devnull, 'w') as fnull: 
            logging.debug(cmd)
            returncode = subprocess.check_call(cmd,
                                               stdout = fnull,
                                               stderr = fnull)
        return (cmd, returncode)

