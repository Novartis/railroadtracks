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
from railroadtracks import core
from railroadtracks import unifex


class UnifexTestCase(unittest.TestCase):
    def test__extract_argdict(self):
        # invalid one arg
        args = ('/bar/baz.fq', )
        self.assertRaises(ValueError, unifex._extract_argdict, args)

        # one arg
        arglist = ('foo=/bar/baz.fq', )
        d = unifex._extract_argdict(arglist)
        self.assertEqual(set(('foo', )), set(d.keys()))

        # two args
        arglist = ('foo=/bar/baz.fq', 'bar=/baz/foo.bam')
        d = unifex._extract_argdict(arglist)
        self.assertEqual(set(('foo', 'bar')), set(d.keys()))

        # three args
        arglist = (('a=123', 'b=abc', 'a=456'))
        argdict = unifex._extract_argdict(arglist)
        # 
        self.assertEqual(set(argdict.keys()), set(('a','b')))
        self.assertEqual(argdict['a'], ['123','456'])
        self.assertEqual(argdict['b'], ['abc',])
        # test round trip
        arglist2 = unifex._extract_arglist(argdict)
        self.assertEqual(list(arglist).sort(), list(arglist2).sort())

    def test_build_AssetSet(self):
        values = (('abc',), ('123','456'))
        AssetSet = core.assetfactory('AssetSet', [core.AssetAttr('afile', core.File,''),
                                                  core.AssetAttr('asequenceoffiles', core.FileSequence, '')])
        assetset = unifex.build_AssetSet(AssetSet, values)
        #FIXME: where is the test ?
        #FIXME: how are allownone=True assets handled ?

        # values not allowed to be None cannot be None
        values = (None, ('123','456'))
        self.assertRaises(AssertionError, unifex.build_AssetSet, AssetSet, values)

        values = (('abc',), ('123','456'))
        AssetSet = core.assetfactory('AssetSet', [core.AssetAttr('afile', core.File,'', allownone=True),
                                                  core.AssetAttr('asequenceoffiles', core.FileSequence, '')])
        assetset = unifex.build_AssetSet(AssetSet, values)
        values = (None, ('123','456'))
        assetset = unifex.build_AssetSet(AssetSet, values)


    def test_AssetSet_repr(self):
        values = (('abc',), ('123','456'))
        AssetSet = core.assetfactory('AssetSet', [core.AssetAttr('afile', core.File,''),
                                                  core.AssetAttr('asequenceoffiles', core.FileSequence, '')])
        assetset = unifex.build_AssetSet(AssetSet, values)
        res = str(assetset)
        self.assertTrue(isinstance(res, str))
