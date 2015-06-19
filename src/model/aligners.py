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

import sys
from abc import ABCMeta, abstractproperty
import os
import logging
logger = logging.getLogger(__name__)
import subprocess
import re
import tempfile
import shutil
import argparse
import warnings
from railroadtracks.unifex import UnifexError
from railroadtracks import (core,
                            environment)
from railroadtracks.model.files import (FilePattern,
                                        FASTQPossiblyGzipCompressed,
                                        FASTAFile,
                                        SAMFile,
                                        BAMFile,
                                        samtools_getversion,
                                        SamtoolsSamToBam)


class ACTIVITY(core.Enum):
    """ Activities that can be performed by the different steps modeled.
    (note: steps can combine several activities - the most obvious example is
    existing monolithic pipelines) """
    __order__ = 'INDEX ALIGN UTILITY'
    INDEX = 'Index'
    ALIGN = 'Align'
    UTILITY = 'Utility'


class AssetsIndexer(core.AssetsStep):
    """
    Assets for a general indexing step.
    """
    Source = core.assetfactory('Source', 
                               [core.AssetAttr('reference',
                                               FASTAFile, '')])
    Target = core.assetfactory('Target', 
                               [core.AssetAttr('indexfilepattern',
                                               FilePattern, '')])

class IndexerAbstract(core.StepAbstract):
    """
    Parent class for an indexing step
    """
    __metaclass__ = ABCMeta
    activities = (ACTIVITY.INDEX, )
    Assets = AssetsIndexer


class BWAIndexFiles(FilePattern):
    """ Set of BWA index files. """
    _extension = ('.fai', '.rpac', '.amb', '.ann', '.pac', '.bwt', '.rbwt', '.rsa', '.sa')

class SavedBowtieIndex(FilePattern):
    """ Set of Bowtie index files. """
    _extension = ('.ebwt', )


class BowtieBuild(IndexerAbstract):
    """
    Bowtie as an indexer.
    
    Should the way it is run on the command line change, a child class
    will be implemented (and the relevant methods be overriden).
    """
    #FIXME: Bowtie2Build inherit from BowtieBuild ?
    _name = 'bowtie-build'
    _default_execpath = 'bowtie-build'
    class Assets(AssetsIndexer):
        Target = core.assetfactory('Target',
                                   [core.AssetAttr('indexfilepattern', 
                                                   SavedBowtieIndex, '')])

    def __init__(self, executable=None):
        """        
        :param executable: path to the executable
        :type param: :class:`str`
        """
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
                tmp = subprocess.check_output(cmd)
            except OSError as ose:
                raise UnifexError("""Command: %s
                %s""" % (' '.join(cmd), ose))
            m = re.match(b'^.+? version ([^ \n]+).*', tmp)
            self._version = m.groups()[0]
        return str(self._version)

    def run(self, assets, parameters=tuple()):
        """ The step on the assets given as parameters. 
        :param assets: the assets to run the step with.
        :param parameters: optional parameters as a sequence at the exception of the assets
        (e.g. ('-k 24', '--super-fast') ).
        :type parameters: a sequence. """
        assert environment.Executable.ispresent(self._execpath)
        assert isinstance(assets, core.AssetsStep), "The parameter 'assets' must be an %s" % core.AssetsStep.__name__
        basename_index = assets.target.indexfilepattern.name
        cmd = ['%s' % self._execpath]
        cmd.extend(parameters)
        cmd.extend(['-f',
                    assets.source.reference.name, 
                    assets.target.indexfilepattern.name])
        with open(os.devnull, "w") as fnull:
            logger.debug(subprocess.list2cmdline(cmd))
            returncode = subprocess.call(cmd, stdout = fnull, stderr = fnull)

        #FIXME: test return code ?
        return (cmd, returncode)


class SailfishIndex(IndexerAbstract):
    _name = 'sailfish-index'
    _default_execpath = 'sailfish'

    parser = argparse.ArgumentParser(_name)
    parser.add_argument('--kmerSize',
                        dest = 'kmersize',
                        required = True,
                        type = int)
    parser.add_argument('--threads',
                        type = int,
                        default = 2)
    #FIXME: duplication with parameter declaration above. should be a way to refactor this.
    PARAMETERS_DEFAULT = ('--kmerSize', '16', '--threads', '2')

    def __init__(self, executable=None):
        """
        :param executable: the executable is R. If None, the class-level
        attribute :attr:`_default_execpath` will be used.
        :type executable: a :class:`str` or a :class:`environment.R`
        """
        if type(self) == SailfishIndex:
            warnings.warn('Model SailfishIndex is deprecated. Please use SalmonIndex.')
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
                tmp = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
            except OSError as ose:
                raise UnifexError("""Command: %s
                %s""" % (' '.join(cmd), ose))
            m = re.match(b'^version : (.+)', tmp)
            self._version = m.groups()[0]
        return self._version

    def run(self, assets, parameters = PARAMETERS_DEFAULT):
        """ 
        :param assets: 
        :type assests: instance of :class:`core.AssetsStep` (or child class)
        :param targets:
        :type parameters: :class:`tuple`
        """
        assert environment.Executable.ispresent(self._execpath)
        assert isinstance(assets, self.Assets), \
            "The parameter 'assets' must be a '%s', not a '%s'" % (self.Assets.__name__,
                                                                   type(assets).__name__)
        options, unknown = self.parser.parse_known_args(parameters)
        basename_index = assets.target.indexfilepattern.name
        cmd = [self._execpath, 'index']
        cmd.extend(parameters)
        cmd.extend(['--transcripts',
                    assets.source.reference.name, 
                    '--out',
                    assets.target.indexfilepattern.name])
        
        with open(os.devnull, "w") as fnull:
            logger.debug(subprocess.list2cmdline(cmd))
            returncode = subprocess.call(cmd, stdout = fnull, stderr = fnull)
        return (cmd, returncode)


class SalmonIndex(SailfishIndex):
    _name = 'salmon-index'
    _default_execpath = 'salmon'

    parser = argparse.ArgumentParser(_name)
    parser.add_argument('--threads',
                        type = int,
                        default = 2)

    def run(self, assets, parameters = ()):
        """ 
        :param assets: 
        :type assests: instance of :class:`core.AssetsStep` (or child class)
        :param targets:
        :type parameters: :class:`tuple`
        """
        assert environment.Executable.ispresent(self._execpath)
        assert isinstance(assets, self.Assets), \
            "The parameter 'assets' must be a '%s', not a '%s'" % (self.Assets.__name__,
                                                                   type(assets).__name__)
        options, unknown = self.parser.parse_known_args(parameters)
        basename_index = assets.target.indexfilepattern.name
        cmd = [self._execpath, 'index']
        cmd.extend(parameters)
        cmd.extend(['--transcripts',
                    assets.source.reference.name, 
                    '-i',
                    assets.target.indexfilepattern.name])
        
        with open(os.devnull, "w") as fnull:
            logger.debug(subprocess.list2cmdline(cmd))
            returncode = subprocess.call(cmd, stdout = fnull, stderr = fnull)
        return (cmd, returncode)


def _bwa_version(execpath):
    cmd = [execpath, ]
    m = None
    logfile = tempfile.NamedTemporaryFile()
    with open(os.devnull, 'w') as fnull:
        try:
            logger.debug(subprocess.list2cmdline(cmd))
            returncode = subprocess.call(cmd,
                                         stderr=logfile,
                                         stdout=logfile )
            # BWA is returning 1
            assert returncode==1, "BWA should have been failing with return code 1 (and it did not)"
            # now dig the information out (note: the subprocess has closed the file - open it again)
            with open(logfile.name, mode='rb') as logfile:
                for row in logfile:
                    m = re.match(b'^Version: (.+)$', row)
                    if m is not None:
                        break
        except OSError as ose:
            raise UnifexError("""Command: %s
            %s""" % (' '.join(cmd), ose))
    if m is None:
        raise ValueError('The version number could be determined.')
    return m.groups()[0]


class BWAIndex(IndexerAbstract):
    """
    BWA as an indexer.
    
    """
    _name = 'bwa-index'
    _default_execpath = 'bwa'
    class Assets(AssetsIndexer):
        Target = core.assetfactory('Target',
                                   [core.AssetAttr('indexfilepattern', 
                                                   BWAIndexFiles, '')])

    def __init__(self, executable=None):
        """        
        :param executable: path to the executable
        :type param: :class:`str`
        """
        if executable is None:
            executable = type(self)._default_execpath
        self._execpath = executable
        self._version = None

    @property
    def version(self):
        if self._version is None:
            self._version = _bwa_version(self._execpath)
        return self._version

    def run(self, assets, parameters=tuple()):
        """ The step on the assets given as parameters. 
        :param assets: the assets to run the step with.
        :param parameters: optional parameters as a sequence at the exception of the assets
        :type parameters: a sequence. """
        assert environment.Executable.ispresent(self._execpath)
        assert isinstance(assets, core.AssetsStep), "The parameter 'assets' must be an %s" % core.AssetsStep.__name__
        basename_index = assets.target.indexfilepattern.name
        cmd = ['%s' % self._execpath]
        cmd.append('index')
        cmd.extend(parameters)
        cmd.extend(['-p',
                    assets.target.indexfilepattern.name,
                    assets.source.reference.name])
        with open(os.devnull, "w") as fnull:
            logger.debug(subprocess.list2cmdline(cmd))
            returncode = subprocess.call(cmd, stdout = fnull, stderr = fnull)

        #FIXME: test return code ?
        return (cmd, returncode)


class SavedBowtie2Index(FilePattern):
    """ Set of Bowtie2 index files. """
    _extension = ('.bt2', '.bt21')


class Bowtie2Build(IndexerAbstract):
    """
    Bowtie2 as an indexer.
    
    Should the way it is run on the command line change, a child class
    will be implemented (and the relevant methods be overriden).
    """

    _name = 'bowtie2-build'
    _default_execpath = 'bowtie2-build'
    class Assets(AssetsIndexer):
        Target = core.assetfactory('Target',
                                   [core.AssetAttr('indexfilepattern', 
                                                   SavedBowtie2Index, '')])

    def __init__(self, executable=None):
        """        
        :param executable: path to the executable
        :type param: :class:`str`
        """
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
                tmp = subprocess.check_output(cmd)
            except OSError as ose:
                raise UnifexError("""Command: %s
                %s""" % (' '.join(cmd), ose))
            m = re.match(b'^.+? version ([^ \n]+).*', tmp)
            self._version = m.groups()[0]
        return self._version

    def run(self, assets, parameters=tuple()):
        """ The step on the assets given as parameters. 
        :param assets: the assets to run the step with.
        :param parameters: optional parameters as a sequence at the exception of the assets
        (e.g. ('-k 24', '--super-fast') ).
        :type parameters: a sequence. """
        assert environment.Executable.ispresent(self._execpath), "The executable %s could not be found in the PATH" % self._execpath
        assert isinstance(assets, core.AssetsStep), "The parameter 'assets' must be an %s" % core.AssetsStep.__name__
        basename_index = assets.target.indexfilepattern.name
        cmd = ['%s' % self._execpath]
        cmd.extend(parameters)
        cmd.extend(['-f',
                    assets.source.reference.name, 
                    assets.target.indexfilepattern.name])
        with open(os.devnull, "w") as fnull:
            logging.debug(cmd)
            returncode = subprocess.call(cmd, stdout = fnull, stderr = fnull)

        #FIXME: test return code ?
        return (cmd, returncode)


class AssetsAligner(core.AssetsStep):
    Source = core.assetfactory('Source', [core.AssetAttr('indexfilepattern', FilePattern, ''),
                                          core.AssetAttr('read1', FASTQPossiblyGzipCompressed, ''),
                                          core.AssetAttr('read2', FASTQPossiblyGzipCompressed, '', allownone=True)])
    Target = core.assetfactory('Target', [core.AssetAttr('alignment', BAMFile, '')])


class AlignerAbstract(core.StepAbstract):
    """

    """
    __metaclass__ = ABCMeta
    activities = (ACTIVITY.ALIGN, )
    Assets = AssetsAligner

class Bowtie2(AlignerAbstract):
    """
    Bowtie2 as an aligner.
    
    Should the way it is run on the command line change, a child class should
    be implemented (and the relevant methods be overriden).

    bt = Bowtie2(execpath)
    
    :param execpath: path to the executable

    """

    _name = 'bowtie2'
    _default_execpath = 'bowtie2'
    # vanilla aligner:
    # Source = assetfactory('Source', [core.AssetAttr('indexfilepattern', FilePattern, ''),
    #                                  core.AssetAttr('read1', SavedFASTQ, ''),
    #                                  core.AssetAttr('read2', SavedFASTQ, '', allownone=True)])
    class Assets(AssetsAligner):
        Source = core.assetfactory('Source', [core.AssetAttr('indexfilepattern', SavedBowtie2Index, ''),
                                              core.AssetAttr('read1', FASTQPossiblyGzipCompressed, ''),
                                              core.AssetAttr('read2', FASTQPossiblyGzipCompressed, '', allownone=True)])

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
                res = subprocess.check_output(cmd)
            except OSError as ose:
                raise UnifexError("""Command: %s
                %s""" % (' '.join(cmd), ose))

            m = re.match(b'^.+? version ([^ \n]+).*', res)
            if m is None:
                raise RuntimeError('Could not find the version number.')
            self._version = m.groups()[0]
        return self._version

    def run(self, assets, parameters=tuple()):
        # samtools is a dependency. In the absence of features in RRT that would handle secondary depencies,
        # there is a simple assert here to prevent cryptic failures if the software is not on the system.
        samtools = 'samtools'
        assert environment.Executable.ispresent(self._execpath)
        assert environment.Executable.ispresent(samtools)

        source = assets.source
        # FIXME: these checks should be moved to the Assets class !
        index_pattern = re.compile('^.+\.[1-9][0-9]*\.bt2$')
        indexfiles = [name for cls, name in source.indexfilepattern.iterlistfiles() if index_pattern.match(name)]
        if len(indexfiles) == 0:
            raise ValueError("No bowtie2 index files in %s" % source.indexfilepattern)

        # different strings for the command line
        if source.read1 is None or not isinstance(source.read1, core.SavedEntityAbstract):
            raise ValueError("Incorrect value %s for read1" % source.read1)

        if source.read2 is None:
            # single reads
            cmd_sub = ('-U %s' % source.read1.name, )
        else:
            cmd_sub = ('-1', source.read1.name, '-2', source.read2.name)


        unaligned_fn = None
        aligned_fn = None

        # for ext in sorted(assets.target.alignment._extension, key=len, reverse=True):
        #     if assets.target.alignment.name.endswith(ext):
        #         fn = assets.target.alignment.name
        #         fn = fn[:(len(fn)-len(ext))]
        #         unaligned_fn = fn + '_unaligned' + ext
        #         aligned_fn = fn + '_aligned' + ext
        #         break
        # if unaligned_fn is None:
        #     logger.error('Could not build a file name for unaligned reads.')
        #     #FIXME: abort here.

        # if source.read2 is None:
        #     if os.path.exists(unaligned_fn):
        #         logger.warn('File %s already existing.' % unaligned_fn)
        # else:
        #     # "ext" obtained from loop above
        #     for r_i in ('1', '2'):
        #         tmp_fn = re.sub('(.+)(\\%s)$' % ext, '\\1_%s\\2' % r_i, unaligned_fn)
        #         if os.path.exists(tmp_fn):
        #             logger.warn('File %s already existing.' % unaligned_fn)

        # build command line
        # notes:
        #    - q : FASTQ in input
        cmd_align = ['%s' % self._execpath]
        cmd_align.extend(parameters)
        #cmd_align.extend(('--un', unaligned_fn))
        cmd_align.extend(('-x', source.indexfilepattern.name, 
                    '-q'))
        cmd_align.extend(cmd_sub)

        cmd_sam2bam = (samtools, 'view', '-bS', '-')

        aligned_fn = assets.target.alignment.name
        cmd = cmd_align + ['|', ] + list(cmd_sam2bam) + ['>',] + [aligned_fn, ] 
        logger.debug(subprocess.list2cmdline(cmd))
        with open(os.devnull, "w") as fnull, \
             open(aligned_fn, "w") as fh_out, \
             core.Popen(cmd_align, stdout=subprocess.PIPE, stderr=fnull) as p_align:
            # we are done yet: the output should be BAM, not SAM
            with core.Popen(cmd_sam2bam, 
                            stdin=p_align.stdout, 
                            stdout=fh_out, 
                            stderr=fnull) as p_sam2bam:
                p_sam2bam.communicate()[0]
                returncode = p_sam2bam.returncode

        # now merge aligned and unaligned reads
        #cmd_mergebam = (samtools, 'merge', assets.target.alignment.name, aligned_fn, unaligned_fn)
        #try:
        #    subprocess.check_call(cmd_mergebam)
        #finally:
        #    os.unlink(aligned_fn)
        #    os.unlink(unaligned_fn)
        return (cmd, returncode)

class Bowtie(AlignerAbstract):
    """
    Bowtie as an aligner.
    
    Should the way it is run on the command line change, a child class should
    be implemented (and the relevant methods be overriden).

    bt = Bowtie(execpath)
    
    :param execpath: path to the executable

    """

    _name = 'bowtie'
    _default_execpath = 'bowtie'

    class Assets(AssetsAligner):
        Source = core.assetfactory('Source', [core.AssetAttr('indexfilepattern', SavedBowtieIndex, ''),
                                              core.AssetAttr('read1', FASTQPossiblyGzipCompressed, ''),
                                              core.AssetAttr('read2', FASTQPossiblyGzipCompressed, '', allownone=True)])

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
                res = subprocess.check_output(cmd)
            except OSError as ose:
                raise UnifexError("""Command: %s
                %s""" % (' '.join(cmd), ose))

            m = re.match(b'^.+? version ([^ \n]+).*', res)
            if m is None:
                raise RuntimeError('Could not find the version number.')
            self._version = m.groups()[0]
        return self._version

    def run(self, assets, parameters=tuple()):
        # samtools is a dependency. In the absence of features in RRT that would handle secondary depencies,
        # there is a simple assert here to prevent cryptic failures if the software is not on the system.
        samtools = 'samtools'
        assert environment.Executable.ispresent(self._execpath)
        assert environment.Executable.ispresent(samtools)
        
        source = assets.source
        indexfiles = tuple(name for cls, name in source.indexfilepattern.iterlistfiles())
        if len(indexfiles) == 0:
            raise ValueError("No bowtie index files in %s" % source.indexfilepattern)

        # different strings for the command line
        if source.read1 is None or not isinstance(source.read1, core.SavedEntityAbstract):
            raise ValueError("Incorrect value %s for read1" % source.read1)

        if source.read2 is None:
            # single reads
            cmd_sub = (source.read1.name, )
        else:
            cmd_sub = ('-1', source.read1.name, '-2', source.read2.name)

        # Gzip-compression not supported by bowtie1
        if source.read1.iscompressed or (source.read2 is not None and source.read2.iscompressed):
            raise NotImplementedError("Bowtie(1) does not support gzip-compressed FASTQ file. On-the-fly might eventually be implemented in the future.")
        # build command line
        # notes:
        #    - q : FASTQ in input
        cmd_align = ['%s' % self._execpath]
        cmd_align.extend(parameters)
        # unaligned_fn = None
        # aligned_fn = None
        # for ext in sorted(assets.target.alignment._extension, key=len, reverse=True):
        #     if assets.target.alignment.name.endswith(ext):
        #         fn = assets.target.alignment.name
        #         fn = fn[:(len(fn)-len(ext))]
        #         unaligned_fn = fn + '_unaligned' + ext
        #         aligned_fn = fn + '_aligned' + ext
        #         break
        # if unaligned_fn is None:
        #     logger.error('Could not build a file name for unaligned reads.')
        #     #FIXME: abort here.

        # if source.read2 is None:
        #     if os.path.exists(unaligned_fn):
        #         logger.warn('File %s already existing.' % unaligned_fn)
        # else:
        #     # "ext" obtained from loop above
        #     for r_i in ('1', '2'):
        #         tmp_fn = re.sub('(.+)(\\%s)$' % ext, '\\1_%s\\2' % r_i, unaligned_fn)
        #         if os.path.exists(tmp_fn):
        #             logger.warn('File %s already existing.' % unaligned_fn)

        # cmd_align.extend(('--un', unaligned_fn,))
        cmd_align.append('--sam') # otherwise a SAM-like format is produced...
        cmd_align.extend([source.indexfilepattern.name, 
                          '-q'])
        cmd_align.extend(cmd_sub)

        cmd_sam2bam = (samtools, 'view', '-bS', '-')

        aligned_fn = assets.target.alignment.name
        cmd = cmd_align + ['|', ] + list(cmd_sam2bam) + ['>',] + [aligned_fn, ] 
        logger.debug(subprocess.list2cmdline(cmd))
        with open(os.devnull, "w") as fnull, \
             open(aligned_fn, "w") as fh_out, \
             core.Popen(cmd_align, stdout=subprocess.PIPE, stderr=fnull) as p_align:
            # we are done yet: the output should be BAM, not SAM
            with core.Popen(cmd_sam2bam,
                            stdin=p_align.stdout,
                            stdout=fh_out,
                            stderr=fnull) as p_sam2bam:
                p_align.stdout.close()
                p_sam2bam.communicate()[0]
                returncode = p_sam2bam.returncode

        # # now merge aligned and unaligned reads
        # cmd_mergebam = (samtools, 'merge', assets.target.alignment.name, aligned_fn, unaligned_fn)
        # try:
        #     subprocess.check_call(cmd_mergebam)
        # finally:
        #     os.unlink(aligned_fn)
        #     os.unlink(unaligned_fn)

        return (cmd, returncode)


class BWA(AlignerAbstract):
    """
    BWA as an aligner.
    
    Should the way it is run on the command line change, a child class should
    be implemented (and the relevant methods be overriden).

    bt = Bowtie(execpath)
    
    :param execpath: path to the executable

    """

    _name = 'bwa-align'
    _default_execpath = 'bwa'

    class Assets(AssetsAligner):
        Source = core.assetfactory('Source', [core.AssetAttr('indexfilepattern', BWAIndexFiles, ''),
                                              core.AssetAttr('read1', FASTQPossiblyGzipCompressed, ''),
                                              core.AssetAttr('read2', FASTQPossiblyGzipCompressed, '', allownone=True)])

    def __init__(self, executable=None):
        if executable is None:
            executable = type(self)._default_execpath
        self._execpath = executable
        self._version = None

    @property
    def version(self):
        if self._version is None:
            self._version = _bwa_version(self._execpath)
        return self._version

    def run(self, assets, parameters=tuple()):
        # samtools is a dependency. In the absence of features in RRT that would handle secondary depencies,
        # there is a simple assert here to prevent cryptic failures if the software is not on the system.
        samtools = 'samtools'
        assert environment.Executable.ispresent(self._execpath)
        assert environment.Executable.ispresent(samtools)
        
        source = assets.source
        indexfiles = tuple(name for cls, name in source.indexfilepattern.iterlistfiles())
        if len(indexfiles) == 0:
            raise ValueError("No BWA index files in %s" % source.indexfilepattern)

        # different strings for the command line
        if source.read1 is None or not isinstance(source.read1, core.SavedEntityAbstract):
            raise ValueError("Incorrect value %s for read1" % source.read1)
        if source.read2 is None:
            # single reads
            cmd_sub = ('mem', (source.indexfilepattern.name, source.read1.name))
        else:
            cmd_sub = ('mem', (source.indexfilepattern.name, source.read1.name, source.read2.name))

        # Gzip-compression not supported by BWA 
        if source.read1.iscompressed or (source.read2 is not None and source.read2.iscompressed):
            raise NotImplementedError("BWA does not support gzip-compressed FASTQ file. On-the-fly might eventually be implemented in the future.")
        # build command line
        cmd_align = ['%s' % self._execpath]
        cmd_align.append(cmd_sub[0])
        cmd_align.extend(parameters)
        cmd_align.extend(cmd_sub[1])
        
        cmd_sam2bam = (samtools, 'view', '-bS', '-')

        cmd = cmd_align + ['|', ] + list(cmd_sam2bam) + ['>',] + [assets.target.alignment.name, ] 
        logger.debug(subprocess.list2cmdline(cmd))
        with open(os.devnull, "w") as fnull, \
             open(assets.target.alignment.name, "w") as fh_out, \
             core.Popen(cmd_align, stdout=subprocess.PIPE, stderr=fnull) as p_align:
            # we are done yet: the output should be BAM, not SAM
            with core.Popen(cmd_sam2bam,
                            stdin=p_align.stdout,
                            stdout=fh_out,
                            stderr=fnull) as p_sam2bam:
                p_sam2bam.communicate()[0]
                returncode = p_sam2bam.returncode
        return (cmd, returncode)

        cmd.append(assets.target.alignment.name)
        with open(os.devnull, "w") as fnull:
            logger.debug(subprocess.list2cmdline(cmd))
            returncode = subprocess.call(cmd, stdout = fnull, stderr = fnull)
        # at this point, the output is in SAM while when want BAM
        return (cmd, returncode)


def _star_version(execpath):
    # How to get the version number courtesy of Brant Peterson.
    d = tempfile.mkdtemp()
    cmd = [execpath, '--outFileNamePrefix', d + os.path.sep]
    with open(os.devnull, 'w') as fnull:
        try:
            logger.debug(subprocess.list2cmdline(cmd))
            returncode = subprocess.call(cmd,
                                         stderr=fnull,
                                         stdout=fnull )
        except OSError as ose:
            raise UnifexError("""Command: %s
            %s""" % (' '.join(cmd), ose))
    # STAR is returning 102 (woohoo....)
    assert returncode==102, "STAR should have been returning 102 (and it did not)"
    # now dig the information out
    with open(os.path.join(d, "Log.out"), "r") as logfile:
        row = logfile.readline()
    shutil.rmtree(d)
    m = re.match('STAR svn revision compiled=(.+)', row)
    if m is None:
        m = re.match('STAR version=(.+)', row)
        if m is None:
            raise ValueError('The version number could be determined.')
    return m.groups()[0]


# --starindexer-begin
class StarIndex(IndexerAbstract):
    """
    STAR as an indexer.    
    """

    _name = 'star-index'
    _default_execpath = 'STAR'

    def __init__(self, executable=None):
        if executable is None:
            executable = type(self)._default_execpath
        self._execpath = executable
        self._version = None

    @property
    def version(self):
        if self._version is None:
            self._version = _star_version(self._execpath)
        return self._version

    def run(self, assets, parameters=tuple()):
        """ The step on the assets given as parameters. """
        assert environment.Executable.ispresent(self._execpath)
        assert isinstance(assets, core.AssetsStep), "The parameter 'assets' should be an instance of %s" % core.AssetsStep.__name__

        #if not os.path.exists(assets.target.indexfilepattern.name):
        os.mkdir(assets.target.indexfilepattern.name)
        # 
        basename_index = assets.target.indexfilepattern.name
        # command line is like:
        # STAR --runMode genomeGenerate \
        #     --genomeDir /path/to/GenomeDir \
        #     --genomeFastaFiles /path/to/genome/fasta1 /path/to/genome/fasta2        
        cmd = ['%s' % self._execpath]
        cmd.extend(parameters)
        cmd.extend(['--runMode', 'genomeGenerate',
                    '--genomeDir',
                    assets.target.indexfilepattern.name,
                    '--genomeFastaFiles',
                    assets.source.reference.name])
        with open(os.devnull, "w") as fnull:
            logger.debug(subprocess.list2cmdline(cmd))
            returncode = subprocess.call(cmd, stdout = fnull, stderr = fnull)
        #FIXME: test return code ?
        return (cmd, returncode)
# --starindexer-end

class StarAlign(AlignerAbstract):
    """
    STAR as an aligner.    
    """

    _name = 'star-align'
    _default_execpath = 'STAR'

    class Assets(AssetsAligner):
        Target = core.assetfactory('Target', [core.AssetAttr('alignment', BAMFile, '')])

    def __init__(self, executable):
        if executable is None:
            executable = type(self)._default_execpath
        self._execpath = executable
        self._version = None

    @property
    def version(self):
        if self._version is None:
            self._version = _star_version(self._execpath)
        return self._version

    def run(self, assets, parameters=tuple()):
        # samtools is a dependency. In the absence of features in RRT that would handle secondary depencies,
        # there is a simple assert here to prevent cryptic failures if the software is not on the system.
        samtools = 'samtools'
        assert environment.Executable.ispresent(self._execpath)
        assert environment.Executable.ispresent(samtools)

        source = assets.source
        # different strings for the command line
        if source.read1 is None or not isinstance(source.read1, core.SavedEntityAbstract):
            raise ValueError("Incorrect value %s for read1" % source.read1)

        #--genomeDir /dir/STAR/Genome --readFilesIn /dir/*.fastq.gz --readFilesCommand zcat --outFileNamePrefix SampleA --runThreadN 20
        if source.read2 is None:
            # single reads
            cmd_sub = ['--readFilesIn', source.read1.name]
        else:
            cmd_sub = ['--readFilesIn', source.read1.name, source.read2.name]

        if source.read1.iscompressed:
            readfilescmd = 'zcat'
        else:
            readfilescmd = 'cat'

        cmd = [self._execpath]
        cmd.extend(parameters)
        cmd.extend(['--genomeDir', source.indexfilepattern.name, # reference index
                ])
        cmd.extend(cmd_sub)
        
        outdir = tempfile.mkdtemp(prefix='rrt_%s' % self._name,
                                  dir=os.path.dirname(assets.target.alignment.name))
        outfilename_prefix = os.path.join(outdir, 'star_output')
        cmd.extend([
            '--readFilesCommand', readfilescmd, #FIXME: check presence of zcat ?
            '--outFileNamePrefix', outfilename_prefix,
            '--outSAMunmapped', 'Within'
        ])
        with open(os.devnull, "w") as fnull:
            try:
                logger.debug(subprocess.list2cmdline(cmd))
                returncode = subprocess.call(cmd, stdout = fnull, stderr = fnull)
            except TypeError as te:
                print(cmd)
                raise te
        if returncode == 0:
            # success, but we want a specific name for the alignment file
            starname = os.path.join(outfilename_prefix + 'Aligned.out.sam')

            # with open(os.devnull, "w") as fnull, \
            #      open(aligned_fn, "w") as fh_out:
            #     p_align = core.Popen(cmd_align, stdout=subprocess.PIPE, stderr=fnull)
            #     # we are done yet: the output should be BAM, not SAM
            #     p_sam2bam = core.Popen(cmd_sam2bam, stdin=p_align.stdout, stdout=fh_out, stderr=fnull)
            #     p_align.stdout.close()
            #     p_sam2bam.communicate()[0]
            #     returncode = p_sam2bam.returncode

            Assets = SamtoolsSamToBam.Assets
            assets = Assets(Assets.Source(SAMFile(starname)),
                            Assets.Target(BAMFile(assets.target.alignment.name)))
            s2b = SamtoolsSamToBam(samtools)
            s2b_cmd, s2b_returncode = s2b.run(assets)
                #
            #assets.target.bam
                #
                
        return (cmd, returncode)


class GsnapIndex(IndexerAbstract):
    """
    Gsnap as an indexer.
    
    Should the way it is run on the command line change, a child class
    will be implemented (and the relevant methods be overriden).
    """

    _name = 'gsnap-index'
    _default_execpath = 'gmap_build'

    def __init__(self, executable=None):
        """        
        :param execpath: path to the executable
        :type param: :class:`str`
        """
        if executable is None:
            executable = type(self)._default_execpath
        self._execpath = executable
        self._version = None

    @property
    def version(self):
        if self._version is None:
            cmd = [self._execpath]
            with open(os.devnull, "w") as fnull:
                try:
                    logger.debug(subprocess.list2cmdline(cmd))
                    proc = core.Popen(cmd,
                                      stdout=subprocess.PIPE,
                                      stderr=fnull)
                except OSError as ose:
                    raise UnifexError("""Command: %s
                    %s""" % (' '.join(cmd), ose))
            m = None
            for row in proc.stdout:
                m = re.match(b'^.+? version ([0-9-]+)\.$', row)
                if m is not None:
                    break
            if m is None:
                raise Exception('The version number could not be extracted.')
            self._version = m.groups()[0]
        return self._version

    def run(self, assets, parameters=tuple()):
        """ Run the step on the assets given as parameters. 
        :param assets: the assets to run the step with.
        :param parameters: optional parameters as a sequence at the exception of the assets
        :type parameters: a sequence. """
        assert environment.Executable.ispresent(self._execpath)
        assert isinstance(assets, core.AssetsStep), "The parameter 'assets' must be an %s" % core.AssetsStep.__name__

        if not os.path.exists(assets.target.indexfilepattern.name):
            os.mkdir(assets.target.indexfilepattern.name)

        basename_index = assets.target.indexfilepattern.name
        #makefile = 'Makefile.reference'
        #makefile_path = os.path.join(basename_index, makefile)
        cmd = ['%s' % self._execpath]
        cmd.extend(parameters)
        cmd.extend((
            '-d', 'reference',
            '-D', assets.target.indexfilepattern.name,
            assets.source.reference.name))
        with open(os.devnull, "w") as fnull:
            logger.debug(subprocess.list2cmdline(cmd))
            returncode = subprocess.call(cmd, stdout = fnull, stderr = fnull)
        #FIXME: test return code ?
        return (cmd, returncode)


class GsnapAlign(AlignerAbstract):
    """
    Gsnap as an aligner.
    
    Should the way it is run on the command line change, a child class should
    be implemented (and the relevant methods be overriden).
    """

    _name = 'gsnap-align'
    _default_execpath = 'gsnap'

    def __init__(self, executable=None):
        if executable is None:
            executable = type(self)._default_execpath
        self._execpath = executable
        self._version = None

    @property
    def version(self):
        if self._version is None:
            with open(os.devnull, "w") as fnull:
                cmd = [self._execpath]
                try:
                    logger.debug(subprocess.list2cmdline(cmd))
                    proc = core.Popen(cmd,
                                      stdout=fnull,
                                      stderr=subprocess.PIPE)
                except OSError as ose:
                    raise UnifexError("""Command: %s
                    %s""" % (' '.join(cmd), ose))
            row = next(proc.stderr)
            m = re.match(b'^.+ version ([^ ]+) .+$', row)
            self._version = m.groups()[0]
        return self._version

    def run(self, assets, parameters=tuple()):
        # samtools is a dependency. In the absence of features in RRT that would handle secondary depencies,
        # there is a simple assert here to prevent cryptic failures if the software is not on the system.
        samtools = 'samtools'
        assert environment.Executable.ispresent(self._execpath)
        assert environment.Executable.ispresent(samtools)

        # directory for the index
        assert isinstance(assets, core.AssetsStep), "The parameter 'assets' must be an %s" % core.AssetsStep.__name__

        basename_index = assets.source.indexfilepattern.name

        # different strings for the command line
        gunzip = False
        if assets.source.read2 is None:
            # single reads
            cmd_sub = (assets.source.read1.name, )
            if assets.source.read1.iscompressed:
                gunzip = True
        else:
            cmd_sub = (assets.source.read1.name, assets.source.read2.name)
            if assets.source.read1.iscompressed:
                if assets.source.read2.iscompressed:
                    gunzip = True
                else:
                    raise ValueError('With GSNAP FASTQ can either be all gzip-compressed or all gzip-uncompressed.')
            else:
                # gunzip is already set to False
                pass

        # build command line
        cmd_align = [self._execpath,]
        cmd_align.extend(parameters)
        cmd_align.extend((
            '-d', 'reference', #FIXME: move to parameters + way to have default parameters clearly exposed to higher-level layers
            '-A', 'sam',
            '-D', basename_index))
        if gunzip:
            cmd_align.append('--gunzip')
        cmd_align.extend(cmd_sub)

        cmd_sam2bam = (samtools, 'view', '-bS', '-')

        cmd = cmd_align + ['|', ] + list(cmd_sam2bam) + ['>',] + [assets.target.alignment.name, ] 
        logger.debug(subprocess.list2cmdline(cmd))
        with open(os.devnull, "w") as fnull, \
             open(assets.target.alignment.name, "w") as fh_out, \
             core.Popen(cmd_align, stdout=subprocess.PIPE, stderr=fnull) as p_align:
            # we are done yet: the output should be BAM, not SAM
            with core.Popen(cmd_sam2bam,
                            stdin=p_align.stdout,
                            stdout=fh_out,
                            stderr=fnull) as p_sam2bam:
                p_align.stdout.close()
                p_sam2bam.communicate()[0]
                returncode = p_sam2bam.returncode
        return (cmd, returncode)

class TopHat(core.StepAbstract):
    """
    TopHat
    th = TopHat(execpath)
    
    :param execpath: path to the executable

    """

    _name = 'tophat'
    _default_execpath = 'tophat'
    activities = (ACTIVITY.ALIGN, )

    class Assets(AssetsAligner):
        Source = core.assetfactory('Source', [core.AssetAttr('indexfilepattern', SavedBowtieIndex, ''),
                                              core.AssetAttr('read1', FASTQPossiblyGzipCompressed, ''),
                                              core.AssetAttr('read2', FASTQPossiblyGzipCompressed, '', allownone=True)])
        Target = core.assetfactory('Target', [core.AssetAttr('alignment', BAMFile, '')])

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
                res = subprocess.check_output(cmd)
            except OSError as ose:
                raise UnifexError("""Command: %s
                %s""" % (' '.join(cmd), ose))

            m = re.match(b'^.+? v([0-9][^ \n]+).*', res)
            if m is None:
                raise RuntimeError('Could not find the version number.')
            self._version = m.groups()[0]
        return self._version

    def run(self, assets, parameters=tuple()):
        assert environment.Executable.ispresent(self._execpath)
        source = assets.source
        indexfiles = tuple(name for cls, name in source.indexfilepattern.iterlistfiles())
        if len(indexfiles) == 0:
            raise ValueError("No bowtie1 index files in %s" % source.indexfilepattern)

        # FIXME: checks already done at the asset level. Code clean up to consider...
        # different strings for the command line
        if source.read1 is None or not isinstance(source.read1, core.SavedEntityAbstract):
            raise ValueError("Incorrect value %s for read1" % source.read1)

        if source.read2 is None:
            # single reads
            cmd_sub = (source.read1.name, )
        else:
            cmd_sub = (source.read1.name, source.read2.name)

        # build command line
        # notes:

        cmd = ['%s' % self._execpath]
        cmd.extend(parameters)
        # TopHat saves results into a directory, and the alignment is called "accepted_hits.bam"
        output_dir = os.path.dirname(assets.target.alignment.name)
        if output_dir == '':
            output_dir = '.'
        cmd.extend(('-o', output_dir))
        cmd.append(source.indexfilepattern.name)
        cmd.extend(cmd_sub)
        with open(os.devnull, "w") as fnull:
            logger.debug(subprocess.list2cmdline(cmd))
            returncode = subprocess.check_call(cmd, stdout = fnull, stderr = fnull)        
        source_name = 'accepted_hits.bam'
        if source_name != assets.target.alignment.name:
            #FIXME: test for accidental overwrite ?
            #if os.path.exists(assets.target.alignment.name):
            shutil.move(os.path.join(output_dir, source_name), assets.target.alignment.name)
        return (cmd, returncode)

class TopHat2(TopHat):
    """
    TopHat2
    th = TopHat2(execpath)
    
    :param execpath: path to the executable

    """

    _name = 'tophat2'
    _default_execpath = 'tophat2'

    class Assets(TopHat.Assets):
        Source = core.assetfactory('Source', [core.AssetAttr('indexfilepattern', SavedBowtie2Index, ''),
                                              core.AssetAttr('read1', FASTQPossiblyGzipCompressed, ''),
                                              core.AssetAttr('read2', FASTQPossiblyGzipCompressed, '', allownone=True)])


    def run(self, assets, parameters=tuple()):
        assert environment.Executable.ispresent(self._execpath)
        source = assets.source
        #FIXME: TopHat2 can use bowtie or bowtie2
        indexfiles = tuple(name for cls, name in source.indexfilepattern.iterlistfiles())
        if len(indexfiles) == 0:
            raise ValueError("No bowtie2 index files in %s" % source.indexfilepattern)

        # FIXME: checks already done at the asset level. Code clean up to consider...
        # different strings for the command line
        if source.read1 is None or not isinstance(source.read1, core.SavedEntityAbstract):
            raise ValueError("Incorrect value %s for read1" % source.read1)

        if source.read2 is None:
            # single reads
            cmd_sub = (source.read1.name, )
        else:
            cmd_sub = ('%s,%s' % (source.read1.name, source.read2.name), )

        # build command line
        # notes:

        cmd = ['%s' % self._execpath]
        cmd.extend(parameters)
        # TopHat saves results into a directory, and the alignment is called "accepted_hits.bam"
        output_dir = os.path.dirname(assets.target.alignment.name)
        if output_dir == '':
            output_dir = '.'
        cmd.extend(('-o', output_dir))
        cmd.append(source.indexfilepattern.name)
        cmd.extend(cmd_sub)
        with open(os.devnull, "w") as fnull:
            logger.debug(subprocess.list2cmdline(cmd))
            returncode = subprocess.check_call(cmd, stdout = fnull, stderr = fnull)        
        source_name = 'accepted_hits.bam'
        if source_name != assets.target.alignment.name:
            #FIXME: test for accidental overwrite ?
            #if os.path.exists(assets.target.alignment.name):
            shutil.move(os.path.join(output_dir, source_name), assets.target.alignment.name)
        return (cmd, returncode)


class SamtoolsExtractUnaligned(core.StepAbstract):
    """
    Extract unaligned reads
    
    :param execpath: path to the executable

    """

    _name = 'samtools-extractunaligned'
    _default_execpath = 'samtools'


    class Assets(core.AssetsStep):
        """
        Assets for :class:`SamtoolsExtractUnaligned`
        
        """
        Source = core.assetfactory('Source', [core.AssetAttr('bamfile', BAMFile, '')])
        Target = core.assetfactory('Target', [core.AssetAttr('unaligned', BAMFile, '')])

    activities = (ACTIVITY.UTILITY, )

    def __init__(self, executable=None):
        warnings.warn('SamtoolsExtracUnaliged is deprecated, use SamtoolsFilter', DeprecationWarning)
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
        assert environment.Executable.ispresent(self._execpath)
        cmd = [self._execpath, 'view', '-b', '-f', '4', 
               '-o', assets.target.unaligned.name, assets.source.bamfile.name, ]
        with open(os.devnull, 'w') as fnull: 
            logging.debug(cmd)
            returncode = subprocess.check_call(cmd,
                                               stdout = fnull,
                                               stderr = fnull)
        return (cmd, returncode)
