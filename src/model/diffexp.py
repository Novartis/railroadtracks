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
Differential expression
"""

import os
import argparse
from railroadtracks import core
from railroadtracks import environment
from railroadtracks.model.files import CSVFile, File
from railroadtracks.unifex import _cmdfromuei

RSOURCES_DIR = os.path.dirname(__file__)

class ACTIVITY(core.Enum):
    """ Activities that can be performed by the different steps modeled.
    (note: steps can combine several activities - the most obvious example is
    existing monolithic pipelines) """
    __order__ = 'DIFFEXP'
    DIFFEXP = 'Differential Expression'

class AssetsDifferentialExpression(core.AssetsStep):
    Source = core.assetfactory('Source', [core.AssetAttr('counttable_fn', CSVFile, ''),
                                          core.AssetAttr('sampleinfo_fn', File, '')])
    Target = core.assetfactory('Target', [core.AssetAttr('diffexp_fn', File, '')])

class DifferentialExpressionMeasurer(core.StepAbstract):
    Assets = AssetsDifferentialExpression
    pass

class RDifferentialExpressionMeasurer(DifferentialExpressionMeasurer):

    activities = (ACTIVITY.DIFFEXP, )
    _default_execpath = 'R'
    parser = None

    def __init__(self, executable = None,
                 rsource_template = None):
        if executable is None:
            executable = type(self)._default_execpath
        if not isinstance(executable, environment.R):
            executable = environment.R(executable)
        self.r = executable
        if rsource_template is None:
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

    def run(self, assets, parameters=(), magicvariable = 'railroadtracks_import'):
        """ 
        :param sources: named tuple 
        :param targets: name tuple
        """
        assert isinstance(assets, AssetsDifferentialExpression), "The parameter 'assets' must be an %s" % AssetsDifferentialExpression.__name__
            
        var_in = {'counttable_fn': assets.source.counttable_fn.name,
                  'sampleinfo_fn': assets.source.sampleinfo_fn.name,
                  'diffexp_fn': assets.target.diffexp_fn.name}

        if self.parser is not None:
            options, unknown = self.parser.parse_known_args(parameters)
            #FIXME: just ignore unknown arguments ?
            for k,v in options.__dict__.items():
                if k in var_in:
                    raise ValueError("The parameter '%s' is already used internally.")
                var_in[k] = v

        with open(self.rsource_template) as template:
            code = template.read()

        code = os.linesep.join([code, 'run(%s)' % magicvariable])
        returncode = self.r.run_snippet(code, var_in=var_in)
        # FIXME: should return the unified execution command line
        uei = core.UnifiedExecInfo(self.r.path, self._name, assets.source, assets.target, parameters, 
                                   None, None # logging_file and logging_level
        )
        cmd = _cmdfromuei(uei)
        return (cmd, returncode)

# -- note-R-differential-expression-begin

class EdgeR(RDifferentialExpressionMeasurer):
    _rscript_name = 'edger.R'
    _rpackagename = 'edgeR'
    _name = 'edger'

class DESeq(RDifferentialExpressionMeasurer):
    _rscript_name = 'deseq.R'
    _rpackagename = 'DESeq'
    _name = 'deseq'
    parser = argparse.ArgumentParser(_name)
    parser.add_argument('--dispersion-fittype',
                        dest = 'dispersion_fittype',
                        choices = ("parametric", "local"),
                        default = "parametric")

class DESeq2(DESeq):
    _rscript_name = 'deseq2.R'
    _rpackagename = 'DESeq2'
    _name = 'deseq2'

class LimmaVoom(RDifferentialExpressionMeasurer):
    _rscript_name = 'limmavoom.R'
    _rpackagename = 'limma'
    _name = 'limma-voom'

# -- note-R-differential-expression-end

class EBSeq(RDifferentialExpressionMeasurer):
    _rscript_name = 'ebseq.R'
    _rpackagename = 'EBSeq'
    _name = 'ebseq'
    parser = argparse.ArgumentParser(_name)
    parser.add_argument('--maxround',
                        dest = 'maxround',
                        help = 'Parameter maxround in EBSeq::EBTest(). Default: %(default)s',
                        type = float,
                        default = 5)
    parser.add_argument('--qtrm',
                        dest = 'qtrm',
                        help = 'Parameter Qtrm in EBSeq::EBTest(). Default: %(default)s',
                        type = float,
                        default = 0.75)
    parser.add_argument('--qtrmcut',
                        dest = 'qtrmcut',
                        help = 'Parameter QtrmCut in EBSeq::EBTest(). Default: %(default)s',
                        type = int,
                        default = 10)
