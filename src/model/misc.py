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

import os
import subprocess
import logging
logger = logging.getLogger(__name__)

from railroadtracks import core

from railroadtracks.model.files import FilePattern

class ACTIVITY(core.Enum):
    """ Activities that can be performed by the different steps modeled.
    (note: steps can combine several activities - the most obvious example is
    existing monolithic pipelines) """
    __order__ = 'MISC'
    MISC = 'Miscellaneous'

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
    activities = (ACTIVITY.MISC, )

    def __init__(self, executable=None):
        if executable is None:
            executable = type(self)._default_execpath
        self._execpath = executable
        # no third-party executable involved
        if self._execpath is not None:
            raise ValueError('No third-party executable involved. The executable cannot be None')

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
