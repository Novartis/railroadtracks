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
Unified execution layer.

One general way to run things on the command line.
"""

import sys, os, argparse
import collections
import logging
logger = logging.getLogger(__name__)
import railroadtracks
from importlib import import_module
from . import environment
from . import core


LOG_DEBUG = 'DEBUG'
LOGGING_LEVELS = ('INFO', LOG_DEBUG)

class UnifexError(Exception):
    pass



#FIXME: several classes 'Call' in rrt, I think. Delete/rename others
class Call(object):
    """ Unified call, turning a step + assets + parameters into a task. """
    def __init__(self, step, assets, parameters):
        self._step = step
        self._assets = assets
        self._parameters = parameters

    def execute(self):
        """ Execute the task. """
        return self._step.run(self._assets, self._parameters)

    @property
    def step(self):
        return self._step
    @property
    def assets(self):
        return self._assets
    @property
    def parameters(self):
        return self._parameters

    def __str__(self):
        res = (super(Call, self).__str__(),
               '  step: %s' % str(self.step),
               '  assets: %s' % str(self.assets),
               '  parameters: %s' % str(self.parameters))
        return os.linesep.join(res)

def _set_logging(args):
    if args.logging_file is not None:
        logging.basicConfig(filename=args.logging_file,
                            level=getattr(logging, args.logging_level))

def _make_stepdict(module):
    """ 
    Make a dict of classes (key=step name, value=class).
    Duplicate step names will raise a ValueError.
    :param classlist: module 
    :rtype: :class:`namespace`
    """
    classlist = core.steplist(module)
    d = dict()
    for cls in classlist:
        if not issubclass(cls, (core.StepAbstract,)):
            raise ValueError("Classes used as step must inherit from core.StepAbstract.")
        stepname = cls._name
        if stepname is None:
            raise ValueError('The step name for class "%s" is not defined.' % cls)
        elif stepname in d:
            raise ValueError('The step name "%s" is defined twice.' % stepname)
        d[stepname] = cls
    return d


def _extract_argdict(arglist):
    """ Build dictionary from list of parameters.
    """
    argdict = collections.defaultdict(list)
    if arglist is not None:
         for src in arglist:
            key, value = src.split('=', 1)
            argdict[key].append(value)
    return argdict

def _extract_arglist(argdict):
    """ Build list of parameters from a mapping of parameters.
    """
    arglist = list()
    if argdict is not None:
        for key, value in argdict.items():
            for v in value:
                src = '='.join((key, v))
                arglist.append(src)
    return tuple(arglist)


#FIXME: moved to unifex.py ?
def _cmdfromuei(uei):
    cmd = [sys.executable, '-m', 'railroadtracks.unifex', 'run', uei.model]
    # In the unified command line, the executable does not have to be specified,
    # and when parsing command line arguments with argparse the associated
    # parameter is set to None.
    # When building a unified command line, we skip the executable if it is
    # set to None.
    if uei.executable is not None:
        cmd.append(uei.executable)
    cmd.append('-s')
    cmd.extend('%s=%s' % (x,y.name) for x,y in zip(uei.source._fields, uei.source))
    cmd.append('-t')
    cmd.extend('%s=%s' % (x,y.name) for x,y in zip(uei.target._fields, uei.target))
    if len(uei.parameters) > 0:
        cmd.append('-p')
        cmd.extend(" '%s'" % x for x in uei.parameters)
    return cmd


#FIXME: needed ?
# def wrapper_run(model_cls, sources, targets,
#                 parameters):
#     executable = model_cls(args.executable, parameters = args.parameters)
#     res = executable.run(sources, targets)
#     return res


def _model_instance(args, steplist):
    """
    :param args: arguments, such as the ones returned by :meth:`argparse.ArgumentParser.parse`.
    
    This should contain 2 attributes:
    - :attr:`model`, the name of the model used for execution,
    - :attr:`executable` the name (or path) to the executable used.

    :param steplist: sequence of known steps
    :rtype: :class:``
    """
    model_cls = steplist[args.model]
    try:
        model_instance = model_cls(args.executable)
    except environment.MissingSoftware:
        raise UnifexError("""The executable '%s' to use with the model class '%s' is not in the ${PATH}.
Use either the full path, or add the executable to your PATH.
        """ % (args.executable, args.model))
    except Exception as e:
        raise UnifexError("""Internal error while creating a %s using %s:
%s
""" % (model_cls,
       args.executable,
       e))
    return model_instance


def build_AssetSet(AssetSet, 
                   values):
    """
    :param AssetSet:
    :param values: values to create instances in the AssetSet
    """
    # classes expected for the assets
    assetset_cls = getattr(AssetSet, core.AssetMetaReserved.SOURCES.value)
    assets = list()
    for x, val in zip(assetset_cls, values):
        if val is None:
            # if None, whether it is allowed (allownone True/False) is pushed to the construction
            # of the AssetSet below
            #FIXME: is this really the most transparent way to do it ?
            a = val
        else:
            if issubclass(x.cls, core.FileSequence):
                a = x.cls(x.cls._type(z) for z in val)
            else:
                a = x.cls(val)
        assets.append(a)
    assetset = AssetSet(*assets)
    return assetset

def unified_exec_run(args, steplist, msg=[]):
    """
    Run a command.
    :param args: arguments in a class such as the one returned by :meth:`argparse.ArgumentParser.parse`
    :param steplist: sequence of known steps. The `args` will be matched against this to find the model class.
    :type steplist: sequence of :class:`core.StepAbstract`-inherting instances
    :param msg: list with (eventual) messages
    """
    #FIXME: deprecate msg (use logging)
    _set_logging(args)
    model = _model_instance(args, steplist)
    Assets = model.Assets

    # extract the sources from the args
    sources_dict = _extract_argdict(args.source)
    tmp = ['Source parameters:', ]
    tmp.extend('%s: %s%s' % (k, str(v), os.linesep) for k,v in sources_dict.items())
    logger.debug(os.linesep.join(tmp))
    sources = tuple(sources_dict.get(f) for f in Assets.Source._fields)
    sources_undefined = list()
    for src, assetattr, assetname in zip(sources, Assets.Source._sources, Assets.Source._fields):
        if src is None:
            if not assetattr.allownone:
                sources_undefined.append(assetname)
    if len(sources_undefined) > 0:
        msg.append('The following sources must be defined (and are missing):\n'+\
                   '\n'.join('- %s' % sources_undefined))

    # extract the targets from the args
    targets_dict = _extract_argdict(args.target)
    targets = tuple(targets_dict.get(f) for f in Assets.Target._fields)

    if None in targets:
        msg.append('The following targets must be defined (and are missing):\n'+\
                   '\n'.join('- %s' % (y) for x,y in enumerate(Assets.Target._fields) if targets[x] is None))

    # exit early if missing parameters
    if len(msg) > 0:
        raise ValueError('\n'.join(msg))

    # build the assets
    assets = Assets(build_AssetSet(Assets.Source, sources),
                    build_AssetSet(Assets.Target, targets))
    cleanparameters = tuple(x.lstrip() for x in args.parameters)
    cmd, returncode = model.run(assets, cleanparameters)

    msg = 'return code: %i' % returncode
    logger.debug(msg)
    if returncode != 0:
        sys.stderr.write(msg)
        sys.stderr.write('\n')
    return returncode

def unified_exec_version(args, steplist):
    _set_logging(args)
    executable = _model_instance(args, steplist)
    print(executable.version)
    return 0

def unified_exec_activities(args, steplist):
    _set_logging(args)
    executable = _model_instance(args, steplist)
    print(executable.activities)
    return 0

def unified_exec_model(args, steplist):
    _set_logging(args)
    print(' '.join(str(x) for x in steplist.keys()))
    return 0

def unified_exec():

    def _steplist(module_name):
        module = import_module(module_name)
        steplist = _make_stepdict(module)
        return steplist

    def _unified_exec_run(args):
        msg = []
        try:
            returncode = unified_exec_run(args, _steplist(args.module), msg=msg)
        except ValueError as ve:
            logger.error(msg)
            raise ve
        return returncode

    def _unified_exec_version(args):
        return unified_exec_version(args, _steplist(args.module))
    def _unified_exec_activities(args):
        return unified_exec_activities(args, _steplist(args.module))
    def _unified_exec_model(args):
        return unified_exec_model(args, _steplist(args.module))

    parser = argparse.ArgumentParser(description='Unified model for RNA-Seq steps')
    subparsers = parser.add_subparsers(help='Action to perform. `run` will call the executable specified in the next parameter, `version` will return the version number, `activities` will list the RNA-Seq activities the executable is associated with')

    #FIXME: This means that one will be able to point to different R version, and packages will be coming
    #       from that R version


    # run
    parser_run = subparsers.add_parser('run', help="Run a step")
    #FIXME: This means that one will be able to point to different R version, and packages will be coming
    #       from that R version
    parser_run.add_argument('model',
                            help='Class name used to model that executable. The classes are defined in the model (see -m/--module).')
    parser_run.add_argument('executable',
                            nargs='?', # this is optional
                            help='Name, or full path, for an executable. '
                            'Note that the executable can be R, '
                            'and in that case the class name (argument "model") '
                            'will wrap an R script to be run with that R version.')
    parser_run.add_argument('-m', '--module',
                            default='railroadtracks.rnaseq',
                            help='module defining all models (default: %(default)s)')
    parser_run.add_argument('-s', '--source', nargs='*',
                            help='Source file indicated as <label>=<filename>')
    parser_run.add_argument('-t', '--target', nargs='*',
                            help='target file indicated as <label>=<filename>')
    parser_run.add_argument('-p', '--parameters', nargs='*',
                            default = (),
                            help='Additional parameter(s) for the wrapped executable')
    parser_run.add_argument('--logging-file',
                            help='Name of log file.')
    parser_run.add_argument('--logging-level',
                            choices = LOGGING_LEVELS,
                            default = 'INFO',
                            help='Level of logging information. DEBUG is quite verbose (default: %(default)s).')
    parser_run.set_defaults(func = _unified_exec_run)

    # version
    parser_version = subparsers.add_parser('version',
                                           help="Query version number")
    parser_version.add_argument('model',
                                help='Class name used to model that executable. The classes are defined in the model (see -m/--module).')
    parser_version.add_argument('executable', 
                                nargs='?', # this is optional
                                help='Name, or full path, for an executable. '
                                'Note that the executable can be R, '
                                'and in that case the class name (argument "model") '
                                'will wrap an R script to be run with that R version.')
    parser_version.add_argument('-m', '--module',
                                default='railroadtracks.rnaseq',
                                help='module defining all models (default: %(default)s)')
    parser_version.add_argument('--logging-file',
                                help='Name of log file.')
    parser_version.add_argument('--logging-level',
                                choices = ('INFO', 'DEBUG'),
                                default = 'INFO',
                                help='Level of logging information. DEBUG is quite verbose (default: %(default)s).')

    parser_version.set_defaults(func = _unified_exec_version)

    # activities
    parser_activities = subparsers.add_parser('activities',
                                              help="Query activities")
    parser_activities.add_argument('model',
                                   help='Class name used to model that executable. The classes are defined in the model (see -m/--module).')
    parser_activities.add_argument('-m', '--module',
                                   default='railroadtracks.rnaseq',
                                   help='module defining all models (default: %(default)s)')
    parser_activities.add_argument('executable',
                                   nargs='?', # this is optional
                                   help='Name, or full path, for an executable. '
                                   'Note that the executable can be R, '
                                   'and in that case the class name (argument "model") '
                                   'will wrap an R script to be run with that R version.')
    parser_activities.add_argument('--logging-file',
                                   help='Name of log file.')
    parser_activities.add_argument('--logging-level',
                                   choices = ('INFO', 'DEBUG'),
                                   default = 'INFO',
                                   help='Level of logging information. DEBUG is quite verbose (default: %(default)s).')

    parser_activities.set_defaults(func = _unified_exec_activities)

    # model
    parser_model = subparsers.add_parser('model',
                                         help="Query on a model module")
    parser_model.add_argument('-m', '--module',
                              default='railroadtracks.rnaseq',
                              help='module defining all models (default: %(default)s)')
    parser_model.add_argument('--logging-file',
                              help='Name of log file.')
    parser_model.add_argument('--logging-level',
                              choices = ('INFO', 'DEBUG'),
                              default = 'INFO',
                              help='Level of logging information. DEBUG is quite verbose (default: %(default)s).')

    parser_model.set_defaults(func = _unified_exec_model)
    
    args = parser.parse_args()
    returncode = args.func(args)


if __name__ == '__main__':
    import sys
    returncode = unified_exec()
    sys.exit(returncode)
