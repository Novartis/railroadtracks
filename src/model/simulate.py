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
Very minimal set of functions to simulate NGS data.

There are much more advanced tools available. This is mostly here to run unit-tests and
related activities.

"""

import railroadtracks
import warnings
import argparse
import logging
logger = logging.getLogger(__name__)
from railroadtracks import environment, core, unifex
from railroadtracks.model.files import (File, 
                                        FilePattern,
                                        GFFFile,
                                        FASTQPossiblyGzipCompressed,
                                        BEDFile)
import sys, os, argparse, tempfile, re
import random
import collections
from collections import deque
import csv
import subprocess
import shutil

if sys.version_info[0] < 3:
    linesep = os.linesep
else:
    linesep = bytes(os.linesep, 'ascii')

PHAGEFASTA = os.path.join(os.path.dirname(railroadtracks.__file__), 'EF204940.fa')
PHAGEGFF = os.path.join(os.path.dirname(railroadtracks.__file__), 'ef204940.gff')
PHAGEGTF = os.path.join(os.path.dirname(railroadtracks.__file__), 'ef204940.gtf')

# Conditional definition for Python 3:
if sys.version_info[0] > 2:
    maketrans = bytearray.maketrans
else:
    from string import maketrans
COMPLEMENT_TABLE = maketrans(b'ATGC',b'TACG')

Entry = collections.namedtuple('Entry', 'header sequence')
def readfasta_iter(fh):
    """ Iterate through a FASTA file. """
    header = None
    sequence = None
    for row in fh:
        if row.startswith('>'):
            if header is not None:
                yield Entry(header, bytearray(''.join(sequence), 'ascii'))
            header = row.rstrip()
            sequence = list()
        else:
            sequence.append(row.strip())
    if len(sequence) > 0:
        yield Entry(header, bytearray(''.join(sequence), 'ascii'))

GFFEntry = collections.namedtuple('GFFEntry', 'seqname source feature start end score strand fram attribute')
GFFSTART_I = GFFEntry._fields.index('start')
GFFEND_I = GFFEntry._fields.index('end')
def readgff_iter(fh):
    """ Iterate through a GFF file 
    :param fh: iterable (for example a file object). """
    for row in fh:
        if row.startswith('#'):
            continue
        fields = row.split('\t')
        fields[GFFSTART_I] = int(fields[GFFSTART_I])
        fields[GFFEND_I] = int(fields[GFFEND_I])
        yield GFFEntry(*fields)

def entryfrom_gff(entry, gffentry):
    sequence = entry.sequence[gffentry.start:gffentry.end]
    if gffentry.strand == '+':
        pass
    else:
        sequence = sequence.translate(COMPLETEMENT_TABLE)
    return Entry(gffentry.seqname, sequence)

def randomreadstart(entry, length):
    """ random starting point in the sequence of an entry. """
    l = len(entry.sequence)
    if length > l:
        raise ValueError("Reads cannot be longer than the original reference.")
    return random.randint(0, l-length-1)

#FIXME: reverse-complement as well ?
READHEADER_TEMPLATE = '%(init)sFOO:%(lane)i:%(tile)i:%(x)i:%(y)i#0/%(pair)i'
def randomfastq(entry, n, length):
    dummyqual = b'~'*length
    for r_i in range(n):
        r_start = randomreadstart(entry, length)
        if sys.version_info[0] < 3:
            buf = buffer(entry.sequence, r_start, length)
        else:
            # Python 3
            buf = memoryview(entry.sequence)[r_start:(r_start+length)]
        res = (READHEADER_TEMPLATE % {'init':'@', 'x': r_i, 'y': r_i,
                                      'lane':1, 'tile':2, 'pair':1},
               buf,
               READHEADER_TEMPLATE % {'init':'+', 'x': r_i, 'y': r_i,
                                      'lane':1, 'tile':2, 'pair':1},
               dummyqual)
        yield res

if sys.version_info[0] < 3:
    def _readheader(init, r_i, lane, tile, pair):
        return bytes(READHEADER_TEMPLATE % {'init':init, 
                                            'x': r_i, 'y': r_i,
                                            'lane':lane, 'tile':tile, 
                                            'pair':pair})
else:
    def _readheader(init, r_i, lane, tile, pair):
        return bytes(READHEADER_TEMPLATE % {'init':init, 
                                            'x': r_i, 'y': r_i,
                                            'lane':lane, 'tile':tile, 
                                            'pair':pair}, 'ascii')

def randomfastq_pe(entry, n, length, insert, lane=1, tile=2):
    sequence_rc = entry.sequence[::-1].translate(COMPLEMENT_TABLE)
    quality_template = b'wxyz{|}~'
    dummyqual = (quality_template * ((length//len(quality_template))+1))[:length]

    for r_i in range(n):
        r_start = randomreadstart(entry, length+insert+length)
        if sys.version_info[0] < 3:
            buf = buffer(entry.sequence, r_start, length)
        else:
            # Python 3
            buf = memoryview(entry.sequence)[r_start:(r_start+length)]
        read1 = (_readheader('@', r_i, lane, tile, 1),
                 buf,
                 _readheader('+', r_i, lane, tile, 1),
                 dummyqual)
        if sys.version_info[0] < 3:
            buf = buffer(sequence_rc, 0, length)
        else:
            # Python 3
            buf = sequence_rc[:length]
        read2 = (_readheader('@', r_i, lane, tile, 2),
                 buf,
                 _readheader('+', r_i, lane, tile, 2),
                 dummyqual)
        yield (read1, read2)

def randomPEreads(read1_fh, read2_fh,
                  fastaentry, 
                  n=300, l=100, insert=50):
    """
    Random paired-end reads.
    :param read1_fh:
    :param read2_fh:
    :param fastaentry: FASTA entry
    :param n=200: number of reads
    :param l=100: length for the reads
    :return (read1_fh, read2_fh): pair of file handles 
    """
    # FIXME: rather minimalist (plan use of third-party generator of random sequencing
    # data)
    for read1, read2 in randomfastq_pe(fastaentry, n, l, insert):
        for row in read1:
            read1_fh.write(row)
            read1_fh.write(linesep)
        for row in read2:
            read2_fh.write(row)
            read2_fh.write(linesep)
    read1_fh.flush()
    read2_fh.flush()
    return (read1_fh, read2_fh)



# 

class ACTIVITY(core.Enum):
    """ Activities that can be performed by the different steps modeled.
    (note: steps can combine several activities - the most obvious example is
    existing monolithic pipelines) """
    __order__ = 'SIMULATE'
    SIMULATE = 'Simulate'

def _split_mergedpairs(merged_fn, read1_fn, read2_fn):
    # split
    # (inspired by https://gist.github.com/nathanhaigh/3521724)

    inconsistent = False
    if read1_fn.endswith('.gz'):
        if read2_fn.endswith('.gz'):
            pipe_zip = '| gzip '
        else:
            inconsistent = True
    elif read2_fn.endswith('.gz'):
        inconsistent = True
    else:
        pipe_zip = ''
    if inconsistent:
        # this is presumably a mistake. do not allow it
        raise ValueError('The output files can be either both with the extention .gz, or both without.')
        
    env = os.environ.copy()
    env['mergedreads'] = subprocess.list2cmdline((merged_fn,))
    env['read1'] = subprocess.list2cmdline((read1_fn, ))
    env['read2'] = subprocess.list2cmdline((read2_fn,))
    cmd_str = 'paste - - - - - - - - < $mergedreads | ' + \
              'tee >(cut -f 1-4 | tr "\\t" "\\n" %(pipe_zip)s > $read1) | cut -f 5-8 | tr "\\t" "\\n" %(pipe_zip)s > $read2' % locals()
    with core.Popen(('/bin/bash', '-c', 
                     cmd_str),
                    env=env) as proc:
        proc.communicate()
        if proc.returncode != 0:
            raise RuntimeError()


FLUXSIMULATOR_TEMPLATE = """
REF_FILE_NAME   %(annotation_gtf)s
GEN_DIR         %(genome_dir)s

NB_MOLECULES    %(nb_molecules)i

# error model
ERR_FILE        %(error)s

FRAG_SUBSTRATE  RNA
FRAG_METHOD	UR

READ_NUMBER	%(read_number)i
READ_LENGTH	%(read_length)i
PAIRED_END	true

# create a fastq file
FASTA           YES
"""
STATS_FILE_TEMPLATE = """
STATS_FILE   %(stats_fn)s
"""
TMP_DIR_TEMPLATE = """
TMP_DIR   %(tmp_dir)s
"""
LIB_FILE_TEMPLATE = """
LIB_FILE_NAME   %(lib_fn)s
"""
SEQ_FILE_TEMPLATE = """
SEQ_FILE_NAME   %(reads_fn)s
"""
SIZE_DISTRIBUTION_TEMPLATE = """
SIZE_DISTRIBUTION    (%(size_distribution)s)
"""


class FluxSimulatorParameters(File):
    _extension = ('.par', )

class FluxSimulatorParametersAndPro(FilePattern):
    _extension = ('.par', '.pro')

def _fluxsimulator_version(execpath):
    cmd = (execpath, '--version')
    m = None
    logfile = tempfile.NamedTemporaryFile()
    try:
        logger.debug(subprocess.list2cmdline(cmd))
        returncode = subprocess.call(cmd,
                                     stderr=logfile,
                                     stdout=logfile )
        # BWA is returning 1
        assert returncode==1, "BWA should have been failing with return code 1 (and it did not)"
        # now dig the information out (note: the subprocess has closed the file - open it again)
        with open(logfile.name) as logfile:
            for row in logfile:
                m = re.match('^Flux-Simulator v([^ ]+)', row)
                if m is not None:
                    break
    except OSError as ose:
        raise unifex.UnifexError("""Command: %s
        %s""" % (' '.join(cmd), ose))
    if m is None:
        raise ValueError('The version number cannot be determined.')
    version = m.groups()[0]
    return version

from railroadtracks import core
class FluxsimulatorExpression(core.StepAbstract):
    _name = 'flux-simulator-expression'
    _default_execpath = 'flux-simulator'
    _version = None
    activities = (ACTIVITY.SIMULATE, )
    class Assets(core.AssetsStep):
        Source = core.assetfactory('Source', [core.AssetAttr('annotation_gtf',
                                                             GFFFile, ''),
                                              core.AssetAttr('genome_dir',
                                                             FilePattern, ''),
                                              core.AssetAttr('parameters',
                                                             FluxSimulatorParameters, '')])
        Target = core.assetfactory('Target', [core.AssetAttr('parameters_and_pro',
                                                             FluxSimulatorParametersAndPro, '')])
    parser = argparse.ArgumentParser(_name)
    parser.add_argument('--NB-MOLECULES', dest='NB_MOLECULES',
                        required=False,
                        type=int)
    parser.add_argument('--force',
                        action = 'store_true')

    def __init__(self, execpath=None):
        if execpath is None:
            self._execpath = self._default_execpath
        else:
            self._execpath = execpath

    @property
    def version(self):
        if self._version is None:
            self._version = _fluxsimulator_version(self._execpath)
        return self._version
        
    def run(self, assets, parameters=()):
        options = self.parser.parse_args(parameters)
        parameters_out = assets.target.parameters_and_pro.name + '.par'
        with open(assets.source.parameters.name) as fh_in, \
             open(parameters_out, 'w') as fh_out:
            logger.debug('Copying parameter file %s to %s.' % (fh_in.name, fh_out.name))
            for row in fh_in:
                if row.startswith('REF_FILE_NAME'):
                    raise ValueError('The parameter file should not define REF_FILE_NAME')
                fh_out.write(row)
            fh_out.write('REF_FILE_NAME\t%s\n' % assets.source.annotation_gtf.name)
            fh_out.write('GEN_DIR\t%s\n' % assets.source.genome_dir.name)
            if options.NB_MOLECULES is not None:
                fh_out.write('NB_MOLECULES\t%i\n' % options.NB_MOLECULES)
        cmd = [self._execpath, '-p', assets.target.parameters_and_pro.name + '.par', '-x', ]
        cmd.extend(parameters)
        with open(os.devnull, "w") as fnull:
            logger.debug(subprocess.list2cmdline(cmd))
            returncode = subprocess.check_call(cmd, 
                                               stdout=fnull,
                                               stderr=fnull)
        return (cmd, returncode)


class FluxsimulatorSequencing(core.StepAbstract):
    _name = 'flux-simulator-sequencing'
    _default_execpath = 'flux-simulator'
    _version = None
    activities = (ACTIVITY.SIMULATE, )
    class Assets(core.AssetsStep):
        Source = core.assetfactory('Target', [core.AssetAttr('parameters_and_pro',
                                                             FluxSimulatorParametersAndPro, '')])
        Target = core.assetfactory('Target', [core.AssetAttr('parameters_and_pro',
                                                             FluxSimulatorParametersAndPro, ''),
                                              core.AssetAttr('read1',
                                                             FASTQPossiblyGzipCompressed, ''),
                                              core.AssetAttr('read2',
                                                             FASTQPossiblyGzipCompressed, ''),
                                              core.AssetAttr('lib',
                                                             File, '')])
                                          
    parser = argparse.ArgumentParser(_name)
    parser.add_argument('--read-length',
                        dest = 'read_length',
                        default = 150,
                        type = int,
                        help='Read length')
    parser.add_argument('--read-number',
                        dest = 'read_number',
                        type = int,
                        required = True,
                        help='Total number of reads')
    parser.add_argument('--size-distribution',
                        default = 'N(600,200)',
                        help='Distribution of fragment sizes (default: "%(default)s")')
    parser.add_argument('--error-file',
                        dest = 'error_file',
                        default = '76',
                        type = str,
                        help='Error file (possible values by default are "35" or "76")')
    parser.add_argument('--tmp-dir',
                        dest='tmp_dir',
                        default=None)
    parser.add_argument('--force',
                        action = 'store_true')
    parser.add_argument('--label',
                        required = False)
    
    def __init__(self, execpath=None):
        if execpath is None:
            self._execpath = self._default_execpath
        else:
            self._execpath = execpath
    @property
    def version(self):
        if self._version is None:
            self._version = _fluxsimulator_version(self._execpath)
        return self._version
        
    def run(self, assets, parameters=()):
        options, unknown = self.parser.parse_known_args(parameters)
        parameters_out = assets.target.parameters_and_pro.name + '.par'
        ok_pe = False
        missing_pe = True
        temp_dir = tempfile.mkdtemp()
        syntheticreads_prefix = os.path.join(temp_dir, 'mergedpairs') # is it really a bed file ?
        with open(assets.source.parameters_and_pro.name + '.par') as fh_in, \
             open(parameters_out, 'w') as fh_out:
            logger.debug('Copying parameters from %s to %s.' % (fh_in.name, fh_out.name))
            for row in fh_in:
                if row.startswith('LIB_FILE_NAME'):
                    raise ValueError('The parameter file should not define LIB_FILE_NAME')
                if row.startswith('SEQ_FILE_NAME'):
                    raise ValueError('The parameter file should not define SEQ_FILE_NAME')
                if row.startswith('TMP_DIR') and options.tmp_dir is not None:
                    # silently pass - TMP_DIR is written later
                    pass

                m = re.match('PAIRED_END\t(.+)', row)
                if m is not None:
                    missing_pe = False
                    if m.groups()[0] == 'true':
                        ok_pe = True
                fh_out.write(row)
            fh_out.write('\n')
            fh_out.write('FASTA\tYES\n') # although we want FASTQ, but FluxSimulator wants it that way.
            fh_out.write('READ_NUMBER\t%i\n' % options.read_number)
            fh_out.write('READ_LENGTH\t%i\n' % options.read_length)
            fh_out.write('ERR_FILE\t%s\n' % options.error_file)
            fh_out.write(linesep.join((SEQ_FILE_TEMPLATE,
                                       LIB_FILE_TEMPLATE,
                                       SIZE_DISTRIBUTION_TEMPLATE)) % {'reads_fn': syntheticreads_prefix + '.fastq', 
                                                                          'lib_fn': assets.target.lib.name,
                                                                          'size_distribution': options.size_distribution})
            #fh_out.write(STATS_FILE_TEMPLATE % ({'stats_fn': syntheticreads_prefix + '.txt'}))
            if options.tmp_dir is not None:
                fh_out.write(TMP_DIR_TEMPLATE % ({'tmp_dir': options.tmp_dir}))

            if missing_pe:
                #FIXME: allow paired-end to be False ?
                fh_out.write('PAIRED_END\ttrue\n')
        if not ok_pe and not missing_pe:
            raise ValueError('The input parameter file does not set PAIRED_END to true')
        logger.debug('Copying profile from %s to %s.' % (assets.source.parameters_and_pro.name + '.pro',
                                                         assets.target.parameters_and_pro.name + '.pro'))
        shutil.copy(assets.source.parameters_and_pro.name + '.pro',
                    assets.target.parameters_and_pro.name + '.pro')
        cmd = [self._execpath, '-p', parameters_out, '-l', '-s']
        cmd.extend(parameters)
        #FIXME: channel the ouput into a log file instead ?
        with open(os.devnull, "w") as fnull:
            logger.debug(subprocess.list2cmdline(cmd))
            returncode = subprocess.check_call(cmd, stdout=fnull, stderr=fnull)
        _split_mergedpairs((syntheticreads_prefix + '.fastq'), assets.target.read1.name, assets.target.read2.name)
        #os.unlink(tmp_reads.name) # not need - Python will remove it when object collected
        return (cmd, returncode)

#


FluxsimulatorProEntry = collections.namedtuple('FluxsimulatorProEntry',
                                               'locus transcript_id coding length expressed_fraction expressed_number')


class FluxsimulatorPro(object):
    def __init__(self, filename, warn=True):
        if warn and not filename.endswith('.pro'):
            warnings.warn('The file "%s" might not be a PRO file.')
        self.name = filename

    def __iter__(self):
        l = len(FluxsimulatorProEntry._fields)
        with open(self.name) as fh:
            csv_r = csv.reader(fh, delimiter="\t")
            for row in csv_r:
                yield FluxsimulatorProEntry(*(row[:l]))


class FluxsimulatorDerivedExpression(core.StepAbstract):
    """
    Derive new expression values from an existing Flux-Simulator .PRO file by shuffling a proportion
    of expression values.

    The proportion to be shuffled can be modified as a parameter.

    """
    _name = 'flux-simulator-derivedexpression'
    _default_execpath = None
    _version = railroadtracks.__version__
    activities = (ACTIVITY.SIMULATE, )
    COLNAMES = ('Locus', 'Transcript_ID', 'Coding', 'Length', 'Expressed Fraction', 'Expressed Number')

    class Assets(core.AssetsStep):
        Source = core.assetfactory('Source', [core.AssetAttr('parameters_and_pro',
                                                             FluxSimulatorParametersAndPro, '')])
        Target = core.assetfactory('Target', [core.AssetAttr('parameters_and_pro',
                                                             FluxSimulatorParametersAndPro, '')])
                                          
    parser = argparse.ArgumentParser(_name)
    parser.add_argument('--percentage',
                        default = 20,
                        type = float,
                        help='Percentage of expression values to change')
    parser.add_argument('--tmp-dir',
                        dest='tmp_dir',
                        default=None)

    def __init__(self, execpath=None):
        if execpath is None:
            self._execpath = self._default_execpath
        else:
            self._execpath = execpath
    @property
    def version(self):
        return self._version
        
    def run(self, assets, parameters=()):
        options, unknown = self.parser.parse_known_args(parameters)
        parameters_out = assets.target.parameters_and_pro.name + '.par'
        ok_pe = False
        missing_pe = True
        temp_dir = tempfile.mkdtemp()
        syntheticreads_prefix = os.path.join(temp_dir, 'mergedpairs') # is it really a bed file ?
        with open(assets.source.parameters_and_pro.name + '.par') as fh_in, \
             open(parameters_out, 'w') as fh_out:
            logger.debug('Copying parameters from %s to %s.' % (fh_in.name, fh_out.name))
            for row in fh_in:
                if row.startswith('LIB_FILE_NAME'):
                    raise ValueError('The parameter file should not define LIB_FILE_NAME')
                if row.startswith('SEQ_FILE_NAME'):
                    raise ValueError('The parameter file should not define SEQ_FILE_NAME')
                if row.startswith('TMP_DIR') and options.tmp_dir is not None:
                    # silently pass - TMP_DIR is written later
                    pass

                m = re.match('PAIRED_END\t(.+)', row)
                if m is not None:
                    missing_pe = False
                    if m.groups()[0] == 'true':
                        ok_pe = True
                fh_out.write(row)
            if options.tmp_dir is not None:
                fh_out.write(TMP_DIR_TEMPLATE % ({'tmp_dir': options.tmp_dir}))

        if not ok_pe and not missing_pe:
            raise ValueError('The input parameter file does not set PAIRED_END to true')
        logger.debug('Deriving a new profile from %s into %s.' % (assets.source.parameters_and_pro.name + '.pro',
                                                                  assets.target.parameters_and_pro.name + '.pro'))
        # First read it all (easier to swap)
        ef_i = self.COLNAMES.index('Expressed Fraction')
        en_i = self.COLNAMES.index('Expressed Number')
        allrows = list()
        with open(assets.source.parameters_and_pro.name + '.pro') as fh_in:
            csv_r = csv.reader(fh_in, delimiter='\t')
            for n_entries, row in enumerate(csv_r):
                allrows.append((n_entries, row))
        # determine which rows will be modified:
        # 1- get the indices for all rows
        sample = list(range(n_entries))
        # 2- shuffle the indices (we will take the first indices in the shuffle)
        random.shuffle(sample)
        # number of entries modified
        n_modified = int(round(n_entries * options.percentage * 1.0 / 100))
        # require an even number (because we are going to swap the expression values)
        if n_modified % 2 != 0:
            n_modified += 1
        # the first indices will be our sample
        sample = sample[:n_modified]
        
        # iterate through pairs of indices to be swapped
        for i in range(0, n_modified, 2):            
            spl_i = sample[i]
            spl_j = sample[i+1]
            # swap the "Expressed Fraction" (ef_i)
            allrows[spl_i][1][ef_i], allrows[spl_j][1][ef_i] = allrows[spl_j][1][ef_i], allrows[spl_i][1][ef_i]
            # swap the "Expressed Number" (en_i)
            allrows[spl_i][1][en_i], allrows[spl_j][1][en_i] = allrows[spl_j][1][en_i], allrows[spl_i][1][en_i]

        with open(assets.target.parameters_and_pro.name + '.pro', 'w') as fh_out:
            csv_w = csv.writer(fh_out, delimiter='\t')
            csv_w.writerows(row for row_i, row in allrows)
        return (None, 0)
