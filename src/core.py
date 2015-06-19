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

import collections
import os
import abc
from operator import attrgetter
import sys
import subprocess

from six import with_metaclass

if sys.version_info[0] < 3:
    class Popen(subprocess.Popen):
        def __enter__(self):
            return self
        def __exit__(self, exc_type, exc_val, exc_tb):
            if self.stdout:
                self.stdout.close()
            if self.stderr:
                self.stderr.close()
            if self.stdin:
                self.stdin.close()
            self.wait()
else:
    # Python 3 !
    unicode = str
    Popen = subprocess.Popen

from flufl.enum import Enum


# Error message 
_NOTIMPLEMENTED_ABSTRACT = "This is an interface, and should be implemented by child classses."


class SavedEntityAbstract(with_metaclass(abc.ABCMeta, object)):
    """
    Represent a file, or group of files, whether already existing on the disk or not.

    Instances of the class are iterable in order to have a unified handling for single
    files and sequences of files.
    """

    _defined = False
    _name = False
    _extension = None

    def __init__(self, *args, **kwargs):
        raise TypeError('The class %s does not have a constructor.' % str(type(self)))

    def isnewerthan(self, savedentityabstract):
        """ """

        assert self._defined

        t_self = self.lastmodified()
        t = savedentityabstract.lastmodified()
        # special cases with None
        if t_self is None:
            # This is a problematic situation
            # savedentityabstract cannot be newer
            # than "self", that is still be materialized
            return False
        elif t is None:
            # t_self cannot be None (already tested above)
            # self is the newest then
            return True
        return t_self > t

    @staticmethod
    def _hash_components(cls, name):
        clsname = cls.__module__ + '.' + cls.__name__
        return (clsname, name)

    @property
    def hashdb(self):
        if not self._defined:
            raise ValueError("The hashdb for an undefined saved entity cannot be computed.")
        return hash(type(self)._hash_components(type(self), self._name))

    @abc.abstractmethod
    def __iter__(self):
        raise NotImplementedError(_NOTIMPLEMENTED_ABSTRACT)

    @abc.abstractmethod
    def __len__(self):
        raise NotImplementedError(_NOTIMPLEMENTED_ABSTRACT)

    # FIXME: does it make any sense to keep this method at the top abstract level,
    #        or should it be cleaner to only have this for containers (sequences) ?
    @abc.abstractmethod
    def iteritems(self):
        raise NotImplementedError(_NOTIMPLEMENTED_ABSTRACT)


def _check_extension(name, extension):
    if (name is not None) and (extension is not None):
        ok = any(name.endswith(x) for x in extension)
        assert ok, \
            "The file name %s should have one of the following extensions: %s" % (name, repr(extension))

class File(SavedEntityAbstract):
    """
    File name, or prospective file name.

    This differs from regular files because it may represent a file that is not yet existing on the disk.
    The reason for this is to able to model the case where the input, or output, from a step is not yet
    existing.

    The class implements __iter__ to allow duck typing with instances of this class and with
    instances of the sibling class FilePattern.
    """

    def __init__(self,
                 name=None):
        """
        :param name: a string, or a sequence of strings with only one string
        """
        if name is not None:
            self._defined = True
            if not (isinstance(name, str) or isinstance(name, unicode)):
                if len(name) == 1:
                    name = next(iter(name))
                else:
                    raise ValueError('"name" should be a string.')
            _check_extension(name, self._extension)
        self._name = name

    def _getname(self):
        assert self._defined
        return self._name

    def _setname(self, name):
        assert not self._defined
        _check_extension(name, self._extension)
        self._name = name
        self._defined = True

    name = property(_getname, _setname, None,
                    'Path to the file. When not defined, this can be None and accessing it will raise an AssertionError (it can only be set then).')

    def lastmodified(self):
        assert self._defined
        if os.path.exists(self.name):
            return os.path.getmtime(self.name)
        else:
            return None

    def __len__(self):
        return 1

    def __iter__(self):
        """ Iterate through _the_ filename. """
        assert self._defined
        yield self.name

    def iteritems(self):
        """ """
        assert self._defined
        yield(type(self), self.name)

class FileSequenceAbstract(SavedEntityAbstract):
    """ Sequence of files """
    _type = SavedEntityAbstract

    @property
    def hashdb(self):
        return hash((super(FileSequenceAbstract, self).hashdb, 
                     tuple(x for x in self)))

    def __init__(self, savedentities):
        savesentities = tuple(savedentities)
        if savedentities is not None:
            # sanity checks
            savedentities = tuple(savedentities)
            for x in savedentities:
                # FIXME: check that there is the _same_ type for all items ? 
                assert isinstance(x, self._type), 'Expected type "%s" but got "%s"' % (str(self._type), str(type(x)))
                assert not isinstance(x, type(self)) # no recursive FileSequence allowed
            self._defined = True
        else:
            pass
            #FIXME: do we allow undefined FileSequences ???
        self._savedentities = savedentities

    def __len__(self):
        return len(self._savedentities)

    def __iter__(self):
        """ Iterate through files in the pattern """
        assert self._defined
        return iter(x.name for x in self._savedentities)

    def iteritems(self):
        assert self._defined
        for se in self._savedentities:
            yield (type(se), se.name)



class FileSequence(FileSequenceAbstract):
    _type = File

#FIXME: why isn't this in unifex.py ? circular reference ?
UnifiedExecInfo = collections.namedtuple('UnifiedExecInfo', 
                                         'executable model source target parameters logging_file logging_level')


# -- stepabstract-begin
class StepAbstract(with_metaclass(abc.ABCMeta, object)):
    """ Abstract parent for steps. """

    # monicker under which the step will be known (must be unique within a step list)
    _name = abc.abstractproperty()
    # default name for executable associated with the class
    _default_execpath = abc.abstractproperty()
    
    # class of assets for the step (must be a child of :class:`AssetsStep`)
    Assets = abc.abstractproperty()

    activities = abc.abstractproperty(None, None, None, 
                                      "Activities associated with the step.")

    version = abc.abstractproperty(None, None, None,
                                   "Version of the executable associated with the step.")


    @property
    def hashdb(self):
        return hash((type(self), self._default_execpath, self.version, self.activities))

    @abc.abstractmethod
    def run(self, assets, parameters=tuple()):
        """ 
        :param assets: Assets to run the step with
        :type assets: :class:`AssetsStep`
        :param parameters: optional parameters
        """
        raise NotImplementedError(_NOTIMPLEMENTED_ABSTRACT)

    def uei(self, assets, parameters = tuple()):
        # FIXME: should return the unified execution command line
        uei = UnifiedExecInfo(self._execpath, self._name,
                              assets.source, assets.target, parameters,
                              None, None) #logging_file and logging_level
        return uei
# -- stepabstract-end

# iterate over the specification of sources
# name: name of the attribute
# cls: expected class for the attribute
# pattern: regular expression
class AssetAttr(object):
    """ Describe what an asset (as a attribute in an :class:`core.AssetsStep`) should be like. """
    __slots__ = ('name', 'cls', 'pattern', 'allownone')
    def __init__(self, name, cls, pattern, allownone=False):
        self.name = name
        self.cls = cls
        self.pattern = pattern
        self.allownone = allownone


#FIXME: _sources is not a good variable name (because an AssetSet is either called Source or Target) 
class AssetMetaReserved(Enum):
    """ Reserved attribute names for AssetMeta type
    (made as a central reference to facilitate refactoring in the choice
    of names proves unlucky later on)"""
    __order__ = 'SOURCES FIELDS'
    SOURCES = '_sources'
    FIELDS = '_fields' # don't touch this one in any case - it is providing compatibility with namedtuples

def _pset_factory(name):
    def func(self, value):
        if getattr(self, name) is not None:
            raise AttributeError("'%s' object attribute '%s' can only be set if uninitialized." % (type(self), name))
        privatename = '_' + name
        setattr(self, privatename, value)

class AssetMeta(type):
    """ Meta-class looking for a property "AssetMetaReserved.SOURCES",
    expected to be a sequence of objects with the attributes

    :param name: name for the asset
    :type name: class:`str`
    :param cls: class inheriting from SavedEntityAbstract
    :type cls: type
    :param pattern: pattern for file name
    :type pattern: class:`str`
    """
    def __new__(cls, name, parents, dct):
        _sources = dct.get(AssetMetaReserved.SOURCES.value)
        if _sources is None:
            _sources = collections.OrderedDict()
            dct[AssetMetaReserved.SOURCES.value] = _sources
        
        reserved = set(AssetMetaReserved.__dict__.keys())
        _fields = list()
        for aa in _sources:
            assert isinstance(aa, AssetAttr)
            privatename = '_'+aa.name
            # sanity check
            if (aa.name in reserved) or (privatename in reserved):
                raise ValueError('The attribute name "%s" is reserved and cannot be used.' % aa.name)
            elif (aa.name in dct) or (privatename in dct):
                # duplicated name
                raise ValueError("The name %s is defined more than once" % aa.name)
            assert isinstance(aa.name, str)
            assert issubclass(aa.cls, SavedEntityAbstract)
            assert isinstance(aa.pattern, str)
            # add to the class for tab-completion (whenever a user is looking for hints),
            # # but do not allow the wrong type to be set
            dct[privatename] =  None
            dct[aa.name] = property(attrgetter(privatename),
                                    None, None)
            # add to _fields (so it is looking like a named tuple)
            _fields.append(aa.name)
        #
        dct[AssetMetaReserved.FIELDS.value] = _fields
        # complete the initialization
        return super(AssetMeta, cls).__new__(cls, name, parents, dct)

    def __repr__(cls):
        res = [super(AssetMeta, cls).__repr__(),
               'The fields for the set of assets are:']
        for f, v in zip(cls._fields, getattr(cls, AssetMetaReserved.SOURCES.value)):
            res.append('    - %s (%s)' % (f, v.cls.__name__))
        return os.linesep.join(res)

#FIXME: implement a deferred check (to allow unspecified asset elements to exist)
class AssetSet(with_metaclass(AssetMeta, object)):
    """ Ordered set of assets """

    @property
    def hashdb(self):
        lst = list()
        for x in self._sources:
            val = getattr(self, x.name)
            if val is None:
                if not x.allownone:
                    raise ValueError("Invalid asset")
                hashsb = val
            else:
                hashdb = val.hashdb
            lst.append((x.name, hashdb))
        return hash(tuple(lst))

    def __init__(self, *args):
        sources = getattr(self, AssetMetaReserved.SOURCES.value)
        assert len(args) == len(sources), 'The following parameter(s) must be specified, and in that order: %s' % \
            str(tuple(x.name for x in sources))
        for s,value in zip(sources, args):
            if (value is None) and (s.allownone):
                 continue
            assert isinstance(value, s.cls), 'The object "%s" was expected to be of type %s (but is a %s)' % \
                (repr(value), repr(s.cls), type(value))
            privatename = '_'+s.name
            setattr(self, privatename, value)
        pass

    def __str__(self):
        l = len(self)
        res = [super(AssetSet, self).__repr__(),
               '  with %i element(s) (-: defined, *: undefined):' % l]
        for i in range(min(5,l)):
            if getattr(self, self._sources[i].name)._defined:
                bullet = '-'
            else:
                bullet = '*'
            res.append('  %s `%s` (a %s)' % (bullet, self._sources[i].name, self._sources[i].cls))
        if l > 5:
            res.append('    ...')
        return os.linesep.join(res)
        
    def __len__(self):
        """
        Number of assets in the set 
        :rtype: :class:`int`
        """
        return len(self._sources)

    def __iter__(self):
        for s in self._sources:
            yield getattr(self, s.name)

    @classmethod
    def createundefined(cls):
        """
        Create an unspecified AssetsSet-inheriting instance.
        This is an helper function for the situations where
        the AssetSet cannot be fully specified (for example, it is
        not yet known what will the target files be.
        
        :param obj: object such as an instance of Source or Target
        :param classname: class name for the tuple
        """
        # loop over AssetAttr elements
        args = list()
        for aa in getattr(cls, AssetMetaReserved.SOURCES.value):
            # aa.name
            # aa.cls
            # aa.pattern
            args.append(aa.cls(None))
        res = cls(*args)
        return res

    @classmethod
    def getassetattr(cls, attrname):
        """ Get the 'AssetAttr' for an attribute name. """
        try:
            i = cls._fields.index(attrname)
        except ValueError as ve:
            raise ValueError('The attribute "%s" is not in "%s".' % (attrname, cls.__name__))
        return cls._sources[i]

def assetfactory(name, savedentities):
    cls = AssetMeta(name, (AssetSet, ), {AssetMetaReserved.SOURCES.value: savedentities})
    return cls

class AssetsStep(object):
    """
    Assets used by any step are split into sources,
    and targets.
    """
    Source = assetfactory('Source', [])
    Target = assetfactory('Target', [])

    @property
    def hashdb(self):
        return hash((self.source.hashdb, self.target.hashdb))

    def __init__(self, source, target=None):
        """
        :param source: 
        :type source: :class:`AssetSet`
        :param target:
        :type target: :class:`AssetSet`
        """
        # sanity check before storing into the instance
        self.source = self._check(source, self.Source)
        if target is None:
            # allow unspecified targets to exists.
            # build the asset as undefined
            self.target = self.Target.createundefined()
        else:
            self.target = self._check(target, self.Target)

    def _check(self, obj, cls):
        """
        :param obj: object such as an instance of Source or Target
        :param cls: class inheriting from :class:`AssetSet`
        """
        # try:
        #     getattr(cls, AssetMetaReserved.SOURCES)
        # except:
        #     import pdb; pdb.set_trace()
        for assetattr, field in zip(getattr(cls, AssetMetaReserved.SOURCES.value), cls._fields):
            if assetattr.allownone:
                continue
            assert hasattr(obj, field), '%s is missing "%s"' % (str(cls), field)
            field_attr = getattr(obj, field)
            assert isinstance(field_attr,
                              SavedEntityAbstract), \
                '%s is not an instance of %s for class %s (got %s)' % \
                (field, SavedEntityAbstract,
                 type(self), type(field_attr))
        return obj

    def __str__(self):
        return os.linesep.join((super(AssetsStep, self).__str__(), 
                                str(self.source), 
                                '---', 
                                str(self.target)))


def steplist(model):
    """ Return the list of steps (as a list of classes inheriting from core.StepAbstract) 
    present in a a model (a namespace such as a module).
    
    If the name _STEPLIST_CLASSES (a list of classe) is present in the model, it will be 
    used as a definition of the exported classes. If not the content of the namespace will
    be inspected and any class inheriting from StepAbstract will be returned.

    """
    if hasattr(model, '_STEPLIST_CLASSES'):
        classlist = model._STEPLIST_CLASSES
    else:
        # catch all steps
        import inspect
        def is_Step(x):
            return inspect.isclass(x) \
                and (not inspect.isabstract(x)) \
                and issubclass(x, (StepAbstract,))

        classlist = inspect.getmembers(model, 
                                       is_Step)
        classlist = ([x[1] for x in classlist])
    return classlist

class DEFAULT_ACTIVITY(Enum):
    MISC = 'Misc.'

class Model(object):
    """ Namespace for a model. """
    def __init__(self, 
                 steplist_classes = tuple(),
                 activity = DEFAULT_ACTIVITY):
        self._STEPLIST_CLASSES = tuple(steplist_classes)
        self.ACTIVITY = activity
        
