# Copyright 2014 Novartis Institutes for Biomedical Research

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
                                        SavedFASTA,
                                        SavedGFF,
                                        FASTQPossiblyGzipCompressed, 
                                        BAMFile,
                                        SavedSAM,
                                        BAMFileSortedByID,
                                        SavedCSV,
                                        SavedTSV,
                                        FilePattern,
                                        SavedCSVSequence,
                                        SavedSAMSequence,
                                        SamtoolsSamToBam,
                                        SamtoolsBamToSam,
                                        BEDFile)
from railroadtracks.model import simulate
from railroadtracks.model.simulate import (FluxsimulatorExpression,
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
                                           TopHat,
                                           TopHat2,
                                           SamtoolsExtractUnaligned)
from railroadtracks.model import diffexp
from railroadtracks.model.diffexp import (EdgeR,
                                          DESeq,
                                          DESeq2,
                                          LimmaVoom)

from railroadtracks.unifex import UnifexError, _cmdfromuei
from abc import ABCMeta, abstractproperty
            

class ACTIVITY(simulate.ACTIVITY, aligners.ACTIVITY, diffexp.ACTIVITY, files.ACTIVITY):
    """ Activities that can be performed by the different steps modeled.
    (note: steps can combine several activities - the most obvious example is
    existing monolithic pipelines) """
    __order__ = 'SIMULATE INDEX ALIGN SORT QUANTIFY NORMALIZE DIFFEXP UTILITY'
    SORT = 'Sort'
    QUANTIFY = 'Quantify'
    NORMALIZE = 'Normalize'



class AssetsCRCHeadTail(core.AssetsStep):
    """
    Assets for :class:`CRCHeadTail`

    """
    Source = core.assetfactory('Source', [core.AssetAttr('file', File, '')])
    Target = core.assetfactory('Target', [core.AssetAttr('crc', SavedCSV, '')])

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
    


class AssetsSorter(core.AssetsStep):
    Source = core.assetfactory('Source', [core.AssetAttr('alignedreads',BAMFile, '')])
    Target = core.assetfactory('Target', [core.AssetAttr('sortedbam', BAMFile, '')])

class SorterAbstract(core.StepAbstract):
    """

    A sorting step.

    """
    __metaclass__ = ABCMeta
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
            self._version = files.samtools_getversion(self._execpath)
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

# -- assetsquantifier-begin
class AssetsQuantifier(core.AssetsStep):
    Source = core.assetfactory('Source', [core.AssetAttr('alignedreads', BAMFile, ''),
                                          core.AssetAttr('annotationfile', SavedGFF, '')])
    Target = core.assetfactory('Target', [core.AssetAttr('counts', SavedCSV, '')])
# -- assetsquantifier-end

# -- quantifierabstract-begin
class QuantifierAbstract(core.StepAbstract):
    __metaclass__ = ABCMeta
    Assets = AssetsQuantifier
    activities = (ACTIVITY.QUANTIFY, )
# -- quantifierabstract-end

# -- htseqcount-begin
class HTSeqCount(QuantifierAbstract):
        
    _name = 'htseqcount'
    _default_execpath = 'htseq-count'
    _separator = '\t'
    # set of parameters to get htseq-count to work with references such
    # as bacteria or viruses (stored as a class attribute for convenience)
    _noexons_parameters = ('--type=CDS', '--idattr=db_xref', '--stranded=no')

    def __init__(self, executable=None):
        if executable is None:
            executable = type(self)._default_execpath
        self._execpath = executable
        self._version = None

    @property
    def version(self):
        if self._version is None:
            cmd = [self._execpath, '-h']
            try:
                logger.debug(subprocess.list2cmdline(cmd))
                res = subprocess.check_output(cmd)
            except OSError as ose:
                raise UnifexError("""Command: %s
                %s""" % (' '.join(cmd), ose))

            m = re.match('^.*( version )?([^ \n]+)\.$', res.split('\n')[-2])
            if m is None:
                raise RuntimeError('Could not find the version number.')
            self._version = m.groups()[1]
        return self._version

    def run(self, assets, parameters = tuple()):
        # FIXME: shouldn't strandedness be a better part of the model ?
        source = assets.source
        sortedbam = source.alignedreads
        if not isinstance(source.alignedreads, BAMFileSortedByID):
            # htseq-count requires sorted entries
            warnings.warn(("The source asset '%s' should ideally be sorted by read IDs. " +\
                           "We are sorting the file; use explicitly a '%s' rather than a '%s' "+\
                           "for better performances, as well as for reproducibility issues "+\
                           "(the sorting will use whatever 'samtools` is first found in the PATH)") \
                          % ("alignedreads", BAMFileSortedByID.__name__, BAMFile.__name__))
            output_dir = os.path.dirname(assets.target.counts.name)
            # temp file name for the sorted output
            sortedbam_fh = tempfile.NamedTemporaryFile(dir=output_dir, suffix=".bam", delete=False)
            # (cleaning temp files handled by Python, except sortedsam)
            # -- sort
            sorter = SamtoolsSorterByID()
            sorter_assets = sorter.Assets(sorter.Assets.Source(source.alignedreads),
                                          sorter.Assets.Target(BAMFile(sortedbam_fh.name)))
            sorter.run(sorter_assets)
            # sanity check:
            if os.stat(sorter_assets.target.sortedbam.name).st_size == 0:
                warnings.warn('The sorted BAM file is empty.')
            sortedbam = sorter_assets.target.sortedbam
        else:
            sortedbam_fh = None

        # BAM to SAM
        cmd_bam2sam = ['samtools', 'view', sortedbam.name]

        # build command line
        cmd_count = [self._execpath, ]
        cmd_count.extend(parameters)
        cmd_count.extend(['-', source.annotationfile.name])
        cmd = cmd_bam2sam + ['|', ] + cmd_count

        logger.debug(subprocess.list2cmdline(cmd))
        with open(os.devnull, "w") as fnull, \
             open(assets.target.counts.name, 'w') as output_fh:
            csv_w = csv.writer(output_fh)
            # HTSeq-count does not use column names in its output, unfortunately,
            # so we correct that
            csv_w.writerow(['ID','count'])
            p_bam2sam = subprocess.Popen(cmd_bam2sam, stdout=subprocess.PIPE, stderr=fnull)
            p_htseq = subprocess.Popen(cmd_count, 
                                       stdin = p_bam2sam.stdout,
                                       stdout = subprocess.PIPE,
                                       stderr = fnull)
            csv_r = csv.reader(p_htseq.stdout, delimiter='\t')
            for row in csv_r:
                csv_w.writerow(row)
                p_htseq.stdout.flush()
            p_htseq.communicate()[0]
        if p_htseq.returncode != 0:
            if sortedbam_fh is not None:
                os.unlink(sortedbam_fh.name)
            raise subprocess.CalledProcessError(p_htseq.returncode, cmd, None)
        if sortedbam_fh is not None:
            os.unlink(sortedbam_fh.name)
        return (cmd, p_htseq.returncode)
# -- htseqcount-end

class SailfishQuant(QuantifierAbstract):
        
    _name = 'sailfish-quant'
    _default_execpath = 'sailfish'

    LIBRARY_SE_DS = "'T=SE:S=U'"
    LIBRARY_PE = "'T=PE:O=><:S=SA'"

    parser = argparse.ArgumentParser(_name)
    parser.add_argument('--libtype',
                        required=True,
                        help='Library type for Sailfish (see the documentation for Sailfish for details).'
                        'Canned values for that parameters are available as class attribute (LIBRARY_*).')

    class Assets(AssetsQuantifier):
        Source = core.assetfactory('Source', 
                                   [SailfishIndex.Assets.Target.getassetattr('indexfilepattern'),
                                    core.AssetAttr('read1', FASTQPossiblyGzipCompressed, ''),
                                    core.AssetAttr('read2', FASTQPossiblyGzipCompressed, '', allownone=True)])
        # 'counts' is coming from QuantifierAbstract
        Target = core.assetfactory('Target', [core.AssetAttr('counts', SavedCSV, ''),
                                              core.AssetAttr('output_dir', FilePattern, '')])
                                              


    def __init__(self, executable=None):
        if executable is None:
            executable = type(self)._default_execpath
        self._execpath = executable
        self._version = None

    @property
    def version(self):
        if self._version is None:
            cmd = [self._execpath, '--version']
            try:
                logger.debug(subprocess.list2cmdline(cmd))
                res = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
            except OSError as ose:
                raise UnifexError("""Command: %s
                %s""" % (' '.join(cmd), ose))

            m = re.match('^version : (.+)$', res)
            if m is None:
                raise RuntimeError('Could not find the version number.')
            self._version = m.groups()[0]
        return self._version

    def run(self, assets, parameters = ('--libtype', LIBRARY_SE_DS, )):
        options, unknown = self.parser.parse_known_args(parameters)
        source = assets.source
        indexfiles = [name for cls, name in source.indexfilepattern.iterlistfiles()]
        if len(indexfiles) == 0:
            raise ValueError("No index files in %s" % source.indexfilepattern)

        # different strings for the command line
        if assets.source.read1 is None or not isinstance(source.read1, core.SavedEntityAbstract):
            raise ValueError("Incorrect value %s for read1" % source.read1)

        m = re.match("^'T=(SE|PE).+'", options.libtype)
        if m is None:
            warnings.warn('Could not extract the library type from "%s".' % options.libtype)
        elif assets.source.read2 is None and m.groups()[0] == 'PE':
            raise ValueError('The library type "%s" is for paired-ends but the "read2" is not defined in the assets.')

        if assets.source.read2 is None:
            # single reads
            cmd_sub = ('-r', source.read1.name)
        else:
            cmd_sub = ('-1', source.read1.name, 
                       '-2', source.read2.name)

        cmd_count = [self._execpath, 'quant']
        #cmd_count.extend(options)
        #cmd_count.extend(unknown)
        cmd_count.extend(parameters)
        cmd_count.extend(('--index', source.indexfilepattern.name))
        if os.path.exists(assets.target.output_dir.name):
            raise ValueError('Sailfish will fail if the output directory "%s" is already present (and it is)' % assets.target.counts.name)
        cmd_count.extend(('--out', assets.target.output_dir.name))
        cmd_count.extend(cmd_sub)
        #FIXME: why can't I make it work without list2cmdline() and shell=True
        with open(os.devnull, 'w') as fnull:
            returncode = subprocess.check_call(subprocess.list2cmdline(cmd_count), shell=True,
                                               stdout = fnull, stderr = fnull)

        # The model is trying to make quantifiers compatible with one an other.
        # We build the target asset `counts`
        with open(assets.target.counts.name, 'w') as fh_out,\
             open(os.path.join(assets.target.output_dir.name, 'quant_bias_corrected.sf')) as fh_in:
            csv_r = csv.reader(fh_in, delimiter='\t')
            csv_w = csv.writer(fh_out)
            for row in csv_r:
                if (len(row) == 0) or row[0].startswith('#'):
                    header = row
                    continue
                break
            csv_w.writerow(['ID', 'count'])
            col_id_i = header.index('# Transcript')
            col_counts_i = header.index('EstimatedNumReads')
            if header != row:
                csv_w.writerow((row[col_id_i], row[col_counts_i]))
            for row in csv_r:
                csv_w.writerow((row[col_id_i], row[col_counts_i]))
        return (cmd_count, returncode)


class FeatureCount(QuantifierAbstract):
    _name = 'featurecount'
    _default_execpath = 'R'
    _rscript_name = 'featurecount.R'
    _rpackagename = 'Rsubread'

    # parser for additional parameters
    #parser = argparse.ArgumentParser('unifex')
    parser = argparse.ArgumentParser(_name)
    # parser.add_argument('action', nargs='?',
    #                     default = 'run')
    # parser.add_argument('model', nargs='?',
    #                     default = _name)
    # parser.add_argument('executable', nargs='?',
    #                     default = _default_execpath)
    parser.add_argument('--paired-ends',
                        dest = 'ispairedend',
                        default = False,
                        action='store_true')
    parser.add_argument('--strand-specific',
                        dest = 'strandspecific',
                        type = int,
                        default = 0,
                        choices = (0,1,2))
    parser.add_argument('--gtf-featuretype',
                        dest = 'gtf_featuretype',
                        default = 'exon',
                        help = "Passed to R's 'Rsubread::featureCounts()' as named parameter 'GTF.featureType' (default: %(default)s).")
    parser.add_argument('--gtf-attrtype',
                        dest = 'gtf_attrtype',
                        default = 'gene_id',
                        help = "Passed to R's 'Rsubread::featureCounts()' as named parameter 'GTF.attrType' (default: %(default)s).")

    # set of parameters to get htseq-count to work with references such
    # as bacteria or viruses (stored as a class attribute for convenience)
    _noexons_parameters = ('--gtf-featuretype=CDS', '--gtf-attrtype=db_xref')
    _pe_parameters = ('--paired-ends', )

    def __init__(self, executable=None):
        """
        :param executable: the executable is R. If None, the class-level
        attribute :attr:`_default_execpath` will be used.
        :type executable: a :class:`str` or a :class:`environment.R`
        """
        if executable is None:
            executable = type(self)._default_execpath
        if not isinstance(executable, environment.R):
            executable = environment.R(executable)
        self.r = executable
        rsource_template = os.path.join(RSOURCES_DIR, 
                                        self._rscript_name)
        assert os.path.isfile(rsource_template), \
            'The needed R script "%s" is not a file' % rsource_template
        self.rsource_template = rsource_template
        self._run_cmd = None
        self._version = None

    _execpath = property(lambda x: x.r.path, None, None)

    @property
    def version(self):
        if self._version is None:
            self._version = self.r.packageversion(self._rpackagename)
        return self._version

    def run(self, assets, parameters = tuple(), magicvariable = 'railroadtracks_import'):
        """ 
        :param assets: 
        :type assests: instance of :class:`core.AssetsStep` (or child class)
        :param targets:
        :type parameters: :class:`tuple`
        """
        assert isinstance(assets, AssetsQuantifier), "The parameter 'assets' must be an %s" % AssetsQuantifier.__name__

        options = self.parser.parse_args(parameters)
        var_in = {'alignedreads_fn': assets.source.alignedreads.name,
                  'annotation_fn': assets.source.annotationfile.name,
                  'results_fn': assets.target.counts.name,
                  'ispairedend': options.ispairedend,
                  'strandspecific': options.strandspecific,
                  'gtf_featuretype': options.gtf_featuretype,
                  'gtf_attrtype': options.gtf_attrtype}

        with open(self.rsource_template) as template:
            code = template.read()

        code = os.linesep.join([code, 'run(%s)' % magicvariable])
        returncode = self.r.run_snippet(code, var_in=var_in)
        # FIXME: should return the unified execution command line
        uei = core.UnifiedExecInfo(self.r.path, self._name, 
                                   assets.source, assets.target, parameters, 
                                   None, None # logging_file and logging_level
        )
        cmd = _cmdfromuei(uei)
        return (cmd, returncode)


        

class AssetsColumnMerger(core.AssetsStep):
    Source = core.assetfactory('Source', [core.AssetAttr('counts', SavedCSVSequence, '')])
    Target = core.assetfactory('Target', [core.AssetAttr('counts', SavedCSV, '')])

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
        out_fh = open(out_counts.name, 'w')
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

        out_fh.close()
        
        cmd = None
        returncode = 0
        return (cmd, returncode)

class AssetsNormalizer(core.AssetsStep):
    Source = core.assetfactory('Source', [core.AssetAttr('count', SavedCSV, '')])
    Target = core.assetfactory('Target', [core.AssetAttr('normalizedcounts', SavedCSV, '')])

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
                                              core.AssetAttr('reference', SavedFASTA, '')])
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

class AssetsAnyscript(core.AssetsStep):
    """
    Assets for :class:`Anyscript`

    """
    Source = core.assetfactory('Source', [core.AssetAttr('script', core.File, ''),
                                          core.AssetAttr('input', core.FileSequence, '')])
    Target = core.assetfactory('Target', [core.AssetAttr('output', FilePattern, '')])

class Anyscript(core.StepAbstract):
    """
    Do anything.
    """
    _name = 'anyscript'
    _default_execpath = None
    Assets = AssetsAnyscript
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
        executable = assets.source.script.name
        # FIXME: check that this is an executable
        assert os.path.exists(executable)
        assert os.path.exists(assets.target.output.name)
        outdir = assets.target.output.name
        assert os.path.isdir(outdir)
        curdir = os.getcwd()
        os.chdir(outdir)
        try:
            cmd = [executable, ]
            cmd.extend(x.name for x in assets.source.input)
            with open(os.devnull, 'w') as fnull:
                subprocess.check_call(cmd, sdtout = fnull, stderr = fnull)
        finally:
            os.chdir(curdir)
        return (cmd, 0)

_STEPLIST_CLASSES = (
    Anyscript,
    FluxsimulatorExpression,
    FluxsimulatorSequencing,
    BWA,
    BWAIndex,
    Bowtie,
    Bowtie2,
    BowtieBuild,
    Bowtie2Build,
    StarIndex,
    StarAlign,
    GsnapIndex,
    GsnapAlign,
    SailfishIndex,
    SailfishQuant,
    HTSeqCount,
    FeatureCount,
    TopHat,
    TopHat2,
    #Cufflinks,
    #Cuffdiff,
    ColumnMerger,
    EdgeR,
    DESeq,
    DESeq2,
    LimmaVoom,
    SamtoolsSorterByID,
    SamtoolsBamToSam,
    SamtoolsSamToBam,
    SamtoolsExtractUnaligned,
    CRCHeadTail,
    PicardCollectAlignmentSummaryMetrics
)

