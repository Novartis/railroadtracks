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

import os, subprocess
import warnings
import abc
import argparse
import tempfile
import logging
import re
import csv
logger = logging.getLogger(__name__)

from railroadtracks import core, environment, unifex
from railroadtracks.model.files import (BAMFile,
                                        BAMFileSortedByID,
                                        SamtoolsSorterByID,
                                        GFFFile,
                                        CSVFile,
                                        FilePattern,
                                        FASTQPossiblyGzipCompressed)

from railroadtracks.model.aligners import SailfishIndex

from railroadtracks.model.diffexp import RSOURCES_DIR

class ACTIVITY(core.Enum):
    """ Activities that can be performed by the different steps modeled.
    (note: steps can combine several activities - the most obvious example is
    existing monolithic pipelines) """
    __order__ = 'QUANTIFY'
    QUANTIFY = 'Quantify'

# -- assetsquantifier-begin
class AssetsQuantifier(core.AssetsStep):
    Source = core.assetfactory('Source', [core.AssetAttr('alignedreads', BAMFile, ''),
                                          core.AssetAttr('annotationfile', GFFFile, '')])
    Target = core.assetfactory('Target', [core.AssetAttr('counts', CSVFile, '')])
# -- assetsquantifier-end

# -- quantifierabstract-begin
class QuantifierAbstract(core.StepAbstract):
    __metaclass__ = abc.ABCMeta
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
                raise unifex.UnifexError("""Command: %s
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
            with core.Popen(cmd_bam2sam, 
                            stdout=subprocess.PIPE,
                            stderr=fnull) as p_bam2sam:
                with p_bam2sam.stdout, core.Popen(cmd_count, 
                                                  stdin = p_bam2sam.stdout,
                                                  stdout = subprocess.PIPE,
                                                  stderr = subprocess.PIPE) as p_htseq:
                    p_bam2sam.stdout.close()
                    # read the output of HTSeq line-per-line
                    csv_r = csv.reader(p_htseq.stdout, delimiter='\t')
                    for row in csv_r:
                        csv_w.writerow(row)
                        p_htseq.stdout.flush()
                    stdout, stderr = p_htseq.communicate()
                if p_htseq.returncode != 0:
                    if sortedbam_fh is not None:
                        os.unlink(sortedbam_fh.name)
                    logger.error(subprocess.list2cmdline(cmd), extra={'stderr': stderr})
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
    TRANSCRIPT_COLNAME = '# Transcript'

    parser = argparse.ArgumentParser(_name)
    parser.add_argument('--libtype',
                        required=True,
                        help='Library type for Sailfish (see the documentation for Sailfish for details).'
                        'Canned values for that parameters are available as class attribute (LIBRARY_*).')
    parser.add_argument('--quantfile',
                        default='quant_bias_corrected.sf',
                        help='Quantification file produced by sailfish')
    parser.add_argument('--count-colname',
                        dest = 'count_colname',
                        default='EstimatedNumReads',
                        help='Column to use to report counts')

    class Assets(AssetsQuantifier):
        Source = core.assetfactory('Source', 
                                   [SailfishIndex.Assets.Target.getassetattr('indexfilepattern'),
                                    core.AssetAttr('read1', FASTQPossiblyGzipCompressed, ''),
                                    core.AssetAttr('read2', FASTQPossiblyGzipCompressed, '', allownone=True)])
        # 'counts' is coming from QuantifierAbstract
        Target = core.assetfactory('Target', [core.AssetAttr('counts', CSVFile, ''),
                                              core.AssetAttr('output_dir', FilePattern, '')])
                                              


    def __init__(self, executable=None):
        if type(self) == SailfishQuant:
            warnings.warn('Model SailfishQuant is deprecated. Please use SalmonAlignQuant.')

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
                raise unifex.UnifexError("""Command: %s
                %s""" % (' '.join(cmd), ose))

            m = re.match('^version : (.+)$', res)
            if m is None:
                raise RuntimeError('Could not find the version number.')
            self._version = m.groups()[0]
        return self._version


    @staticmethod
    def _is_libtype_pe(libtype):
        m = re.match("^'T=(SE|PE).+'", libtype)
        if m is None:
            return None
        elif m.groups()[0] == 'PE':
            return True
        else:
            return False

    def run(self, assets, parameters = ('--libtype', LIBRARY_SE_DS, )):
        options, unknown = self.parser.parse_known_args(parameters)
        source = assets.source
        indexfiles = [name for cls, name in source.indexfilepattern.iterlistfiles()]
        if len(indexfiles) == 0:
            raise ValueError("No index files in %s" % source.indexfilepattern)

        # different strings for the command line
        if assets.source.read1 is None or not isinstance(source.read1, core.SavedEntityAbstract):
            raise ValueError("Incorrect value %s for read1" % source.read1)

        is_pe = self._is_libtype_pe(options.libtype)
        if is_pe is None:
            warnings.warn('Could not extract the library type from "%s".' % options.libtype)
        elif assets.source.read2 is None and is_pe:
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
            raise ValueError('%s will fail if the output directory "%s" is already present (and it is)' % (self._name, assets.target.counts.name))
        cmd_count.extend(('--out', assets.target.output_dir.name))
        cmd_count.extend(cmd_sub)
        #FIXME: why can't I make it work without list2cmdline() and shell=True
        with open(os.devnull, 'w') as fnull:
            returncode = subprocess.check_call(subprocess.list2cmdline(cmd_count), shell=True,
                                               stdout = fnull, stderr = fnull)

        # The model is trying to make quantifiers compatible with one an other.
        # We build the target asset `counts`
        # NOTE: We experienced that the Sailfish family could duplicate entries IDs
        #       (column TRANSCRIPT_COLNAME) while the numerical values differ (and be very close to zero)
        transcript_set = set()
        with open(assets.target.counts.name, 'w') as fh_out,\
             open(os.path.join(assets.target.output_dir.name, options.quantfile)) as fh_in:
            csv_r = csv.reader(fh_in, delimiter='\t')
            csv_w = csv.writer(fh_out)
            for row in csv_r:
                if (len(row) == 0) or row[0].startswith('#'):
                    header = row
                    continue
                break
            csv_w.writerow(['ID', 'count'])
            col_id_i = header.index(self.TRANSCRIPT_COLNAME)
            col_counts_i = header.index(options.count_colname)
            if header != row:
                transcript_id = row[col_id_i]
                transcript_set.add(transcript_id)
                csv_w.writerow((transcript_id, row[col_counts_i]))
            for row in csv_r:
                transcript_id = row[col_id_i]
                if transcript_id in transcript_set:
                    warnings.warn("The ID %s is not unique." % transcript_id)
                    continue
                transcript_set.add(transcript_id)
                csv_w.writerow((transcript_id, row[col_counts_i]))
        return (cmd_count, returncode)

class SalmonAlignQuant(SailfishQuant):
    _name = 'salmon-alignquant'
    _default_execpath = 'salmon'
    LIBRARY_SE_DS = None
    LIBRARY_PE = "IU"
    TRANSCRIPT_COLNAME = "# Name"
    parser = argparse.ArgumentParser(_name)
    parser.add_argument('--libType',
                        dest='libtype',
                        required=True,
                        help='Library type for Salmon (see the documentation for Sal,on for details).'
                        'Canned values for that parameters are available as class attribute (LIBRARY_*).')
    parser.add_argument('--quantfile',
                        default='quant.sf',
                        help='Quantification file produced by salmon')
    parser.add_argument('--count-colname',
                        dest = 'count_colname',
                        default='TPM',
                        help='Column to use to report counts')


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
        cmd = unifex._cmdfromuei(uei)
        return (cmd, returncode)

