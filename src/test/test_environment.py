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

import unittest
import tempfile
from railroadtracks import environment

has_R = environment.Executable.ispresent('R')

class EnvironmentTestCase(unittest.TestCase):

    def test_executable(self):
        ls_exec = environment.Executable('ls')
        #FIXME: technically no test here

    @unittest.skipIf(not has_R,
                     'R is missing')
    def test_R(self):
        r_exec = environment.R('R')
        r_version = r_exec.version
        # missing package
        self.assertRaises(ValueError, r_exec.packageversion, 'foobarbaz')

        version = r_exec.packageversion('stats')
        self.assertTrue(r_version.startswith(version))

    @unittest.skipIf(not (has_R and \
                          environment.R('R').packageversion_or_none('rjson') is not None),
                     'R and its package "rjson" must be present.')
    def test_R_run_snippet(self):
        r_exec = environment.R('R')
        # run a snippet
        magicvariable = 'railroadtracks_import'
        code = """        
        run <- function(p) {
          res <- 1 + p$x
          # output file
          fh_out <- file(p$output_fn, "w")
          write(res, file = fh_out)
          flush(fh_out)
        }
        run(%s)
        """ % magicvariable
        with tempfile.NamedTemporaryFile() as fh_out:
            var_in = {'x': 3,
                      'output_fn': fh_out.name}
            r_exec.run_snippet(code, var_in)
            res = fh_out.read()
        self.assertEqual(b'4', res.rstrip())
