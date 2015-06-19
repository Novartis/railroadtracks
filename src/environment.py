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
Information about the environment this is running on 
"""

import sys, os, subprocess
import tempfile, json
import logging
import re
from railroadtracks import core

logger = logging.getLogger(__name__)


class MissingSoftware(ValueError):
    pass


def _is_executable(path):
        return os.path.isfile(path) and os.access(path, os.X_OK)

def _find_executable(executable):
    path, fname = os.path.split(executable)
    if path:
        if _is_executable(executable):
            return executable
    else:
        for env_path in os.environ["PATH"].split(os.pathsep):
            path = env_path.strip('"')
            candidate_path = os.path.join(path, executable)
            if _is_executable(candidate_path):
                return candidate_path
    return None
class Executable(object):

    def __init__(self, executable):
        """ Name of executable, or full path to executable. """
        if executable.startswith('/'):
            path = executable
        else:
            path = _find_executable(executable)
        if (path is None) or (not os.path.exists(path)):
            raise MissingSoftware('%s cannot be found' % executable)
        self.path = path

    def __repr__(self):
        s = [super(Executable, self).__repr__(),
             '  absolute path: %s' % self.path]
        return '\n'.join(s)

    @staticmethod
    def ispresent(name):
        try:
            e = Executable(name)
        except MissingSoftware:
            return False
        return True

class R(Executable):

    __version = None
    rjsonversion = None

    def getversion(self):
        if self.__version is None:
            with core.Popen((self.path, '--version'), \
                            stdout=subprocess.PIPE) as proc:
                version = proc.stdout.readline()
                self.__version = re.sub(b'R version (.+)', b'\\1', version)
        logger.info('R version is %s' % self.__version)
        return self.__version    
    version = property(getversion, None, None, 'R version')


    def packageversion(self, name):
        """ Check whether an R package is installed and return the version number.
        If the package cannot be loaded (generally because it is not installed),
        an exception :class:`MissingSoftware` is raised.

        :param name: name of the package

        """
        code = """res <- suppressMessages(suppressWarnings(require("%s", quietly=TRUE))); cat(res)"""
        cmd = (self.path, '--slave', '--no-restore', '-e', code % name)

        rsays = subprocess.check_output(cmd).rstrip()
        if rsays == b'FALSE':
            raise MissingSoftware("The R package '%s' is not installed" % name)

        code = 'res <- suppressMessages(packageVersion("%s")); cat(as.character(res))'
        cmd = (self.path, '--slave', '-e', code % name)
        with open(os.devnull, "w") as fnull,\
             core.Popen(cmd, stdout=subprocess.PIPE,
                        stderr=fnull) as proc:
            version = proc.stdout.readline().rstrip()
            tmp = proc.stdout.read()

        logger.info('R package "%s" version is %s' % (name, version))
        return version

    def packageversion_or_none(self, name):
        """ Return the package version, or `None` if the R package cannot be loaded.
        This method is a wrapper around the method `packageversion`.
        """
        try:
            version = self.packageversion(name)
        except MissingSoftware:
            version = None
        return version

    def run_snippet(self, code, var_in,
                    magicvariable = 'railroadtracks_import'):
        """
        :param code: snippet of R code as text
        :param var_in: dict of variables to import into R's global environment before the snippet is evaluated.
                       this is a rather ad-hoc system, and significant performance problems will occur
                       with large data, or complex data structures.

        For example:

        >>> code = 'print(x+1)'
        >>> d = {'x': 2}
        >>> r = R()
        >>> r.run_snippet(code, d)
        [1] 3
        """
        if self.rjsonversion is None:
            # package rjson used to pass context data
            self.rjsonversion = self.packageversion('rjson')

        code_json_extract = """
        # paranoid check
        if (exists("%(magicvariable)s", envir=.GlobalEnv)) {
        cat("*** The variable name '%(magicvariable)s' is reserved !")
        print()
        }
        
        railroadtracks_import <- local({
            suppressWarnings(require("rjson", quietly=TRUE))
            con <- file("%(filename)s", "r")
            res <- rjson::fromJSON(file=con)
            close(con)
            res
        })
        """ 

        if logger.getEffectiveLevel() == logging.DEBUG:
            delete_tempfile = False
        else:
            delete_tempfile = True
        with tempfile.NamedTemporaryFile(delete=delete_tempfile,
                                         mode='wt') as fh_out:
            var_in_json = json.dump(var_in, fh_out)
            fh_out.flush()
            if not delete_tempfile:
                logger.debug('File "%s" not deleted (because logging level is DEBUG)' % fh_out.name)
            #FIXME: test that rjson is installed earlier ?
            code = os.linesep.join([code_json_extract % {'magicvariable': magicvariable, 
                                                         'filename': fh_out.name},
                                    code])
            logger.debug(os.linesep.join(('Running R code snippet:',
                                          code)))
            returncode = self.run(code)
            if returncode != 0:
                raise Exception("R returned non-zero exit status %i" % returncode)
        return returncode

    def run(self, code,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE):
        """ R arbitrary code in a subprocess. 

        :param code: a string with R code

        Returns the return code for the child process.
        """
        cmd = (self.path, '--slave',)
        with open(os.devnull, "w") as fnull, \
             core.Popen(cmd, stdin=stdin,
                        stdout=stdout,
                        stderr=stderr) as proc:
            stdout, stderr = proc.communicate(input=code.encode('ascii'))
            if proc.returncode != 0:
                logger.warning(stderr)
            return proc.returncode
        

    def __repr__(self):
        s = [super(R, self).__repr__(), '  version: %s' % self.version]
        return '\n'.join(s)


