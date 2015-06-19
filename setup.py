#!/usr/bin/env python

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
if (sys.version_info[0] == 2) and (sys.version_info[1] < 7):
    print("Python>=2.7 is required.")
    sys.exit(1)

try:
    from setuptools import setup
except ImportError:
    print("The Python package 'setuptools' is required in order to install this package.")
    sys.exit(1)

from src import __version__

packagename='railroadtracks'
setup(name=packagename,
      version=__version__,
      description='Framework to connect computional steps with an emphasis on RNA-Seq',
      author='Laurent Gautier',
      author_email='laurent.gautier@novartis.com',
      license='License :: OSI Approved :: Apache Software License',
      url='https://github.com/Novartis/railroadtracks',
      package_dir={packagename: 'src'},
      package_data={packagename: ['src/cache.sql',
                                  'src/model/featurecount.R',
                                  'src/model/edge.R',
                                  'src/model/deseq.R',
                                  'src/model/deseq2.R',
                                  'src/model/limmavoom.R',
                                  'src/EF204940.fa',
                                  'src/ef204940.gff',
                                  'src/ef204940.gtf',
                                  'src/templates/*.html']},
      include_package_data=True,
      packages=[packagename, 
                packagename + '.easy',
                packagename + '.model',
                packagename + '.test'],
      install_requires=['flufl.enum', 'networkx', 'jinja2', 'six', 'enum34'] #FIXME: enum34 not needed if Python >= 3.4
     )
