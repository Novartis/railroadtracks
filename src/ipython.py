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

from __future__ import division
import os
from railroadtracks import core, easy
from railroadtracks.hortator import _TASK_DONE
from collections import namedtuple
from jinja2 import Environment, PackageLoader
from railroadtracks import core, hortator, easy
import railroadtracks.easy.taskgraph
import railroadtracks.easy.tasksetgraph
import networkx
try:
    import pygraphviz
    has_pygraphviz = True
except:
    has_pygraphviz = False

env = Environment(loader=PackageLoader('railroadtracks'),
                  extensions = ['jinja2.ext.with_'])



def format_largenumber(value):
    return "{:,.2f}".format(value)
env.filters['format_largenumber'] = format_largenumber


def html_project_view_storage(project, cssclass='rsrsproject'):
    d = easy.dict_project_view_storage(project)
    d['cssclass'] = cssclass
    template = env.get_template('project_view_storage.html')
    res = template.render(d)
    return res
    

def html_project_view_activities(project, done=_TASK_DONE,
                                 cssclass='rsrsproject'):
    template = env.get_template('project_view_activities.html')
    d = easy.dict_project_view_activities(project, done=done)
    d['cssclass'] = cssclass
    res = template.render(d)
    return res


def html_project_view_results(project,
                              cssclass='rsrsproject'):
    d = easy.dict_project_view_results(project)
    d['cssclass'] = cssclass
    template = env.get_template('project_view_results.html')
    res = template.render(d)
    return res



rsrsproject_css = """
table.rsrsproject {
  width: 100%;
  max-width: 40em;
  border-collapse:separate;
  border:solid black 2px;
  -moz-border-radius:6px;
  -webkit-border-radius:6px;
  border-radius:6px;
}
table.rsrsproject thead {
  background-color:rgb(235,235,235);
}
table.rsrsproject th:first-child {
  -moz-border-radius:6px 0 0 0;
  -webkit-border-radius:6px 0 0 0;
  border-radius:6px 0 0 0;
}
table.rsrsproject th:last-child {
  -moz-border-radius: 0 6px 0 0;
  -webkit-border-radius: 0 6px 0 0;
  border-radius: 0 6px 0 0;
}
table.rsrsproject th:only-child {
  -moz-border-radius: 6px 6px 0 0;
  -webkit-border-radius: 6px 6px 0 0;
  border-radius: 6px 6px 0 0;
}
table.rsrsproject .inverse {
  background-color: rgb(50,50,50);
  color: rgb(230,230,230);
}

progress.rsrsproject {
  -webkit-appearance:none;
  -moz-appearance:none;
  appearance:none;
}
progress[value].rsrsproject::-webkit-progress-bar {
  background-color: rgb(230,230,230);
  box-shadow: 0 2px 5px rgba(0,0,0,0.25) inset;
}
progress[value].rsrsproject::-webkit-progress-value {
  background-color: rgb(150,150,250);
}
progress[value].rsrsproject {
  border: none;
  box-shadow: 0 2px 5px rgba(0,0,0,0.25) inset;
}
progress[value].rsrsproject::-moz-progress-bar {
  background-color: rgb(150,150,250);
}
"""

def html_project_view(project,
                      cssclass='rsrsproject'):
    res = ['<style>', rsrsproject_css, '</style>']
    res.extend((x(project, cssclass=cssclass) for x in (html_project_view_storage,
                                                        html_project_view_activities,
                                                        html_project_view_results)))
    return os.linesep.join(res)


def html_assetset_view(assetset,
                       cssclass='rsrsproject'):
    d = easy.dict_assetset_view(assetset)
    d['cssclass'] = cssclass    
    template = env.get_template('assetset_view.html')
    res = template.render(d)
    return res


def html_task_view(task,
                   cssclass='rsrsproject'):
    d = dict()
    d['source_elements'] = easy.dict_assetset_view(task.call.assets.source)['elements']
    d['target_elements'] = easy.dict_assetset_view(task.call.assets.target)['elements']
    d['parameters'] = task.call.parameters
    d['task'] = task
    d['task_class'] = type(task.call.step).__name__
    d['cssclass'] = cssclass    
    template = env.get_template('task_view.html')
    res = template.render(d)
    return res

def node_label(instance):
    if isinstance(instance, hortator.StoredEntity):
        label_ext = "."+instance.entityname.split(".")[-1]
        label_cls = instance.clsname
        if len(label_ext) < len(label_cls):
            label = label_ext
        else:
            label = label_cls
        label = 't:%i - %s' % (instance.id, label)
    elif isinstance(instance, hortator.StoredSequence):
        label_cls = instance.clsname
        label = label_cls
        label = 's:%i - %s' % (instance.id, label)
    else:
        label = 'a:%i - %s' % (instance.task_id, instance.call.step._name)
    return label

def basic_tasknode_func(ag, inst, nodelabel):
    node = ag.get_node(nodelabel)
    node.attr['shape'] = 'box'
    node.attr['tooltip'] = 'parameters: %s' % ' '.join(str(x) for x in inst.call.parameters)
    return node

def default_tasknode_func(ag, inst, nodelabel):
    node = basic_tasknode_func(ag, inst, nodelabel)
    if inst.status == easy._TASK_FAILED:
        node.attr.update(style='filled', fillcolor='black', fontcolor="white")
    elif inst.status == easy._TASK_DONE:
        node.attr.update(style='filled', fillcolor='white', fontcolor="black", penwidth=1.5)
    else:
        node.attr.update(style='filled', fillcolor='white', color="grey", fontcolor="grey")
    return node

def default_assetnode(ag, inst, nodelabel, font_size, pred_done):
    node = ag.get_node(nodelabel)
    node.attr['fontsize'] = int(font_size*.75)
    node.attr['tooltip'] = inst.entityname
    if not pred_done:
        node.attr['fontcolor'] = "grey"
        node.attr['color'] = "grey"
    return node

def default_assetseqnode(ag, inst, nodelabel, font_size, pred_done):
    node = ag.get_node(nodelabel)
    node.attr['fontsize'] = int(font_size*.75)
    if not pred_done:
        node.attr['fontcolor'] = "grey"
        node.attr['color'] = "grey"
    return node

def _get_pred_done(node, digraph):
    pred_done = True
    for pnode in digraph.predecessors(node):
        instance2 = digraph.node[pnode]['instance']
        if not hasattr(instance2, 'status'):
            pred_done = _get_pred_done(pnode, digraph)
        elif instance2.status != hortator._TASK_DONE:
            pred_done = False
            break
    return pred_done


def agraph_fromdigraph(digraph, layout = "dot",
                       nodelabel_func = node_label,
                       tasknode_func = default_tasknode_func,
                       assetnode_func = default_assetnode,
                       assetseqnode_func = default_assetseqnode,
                       graph_size = None,
                       font_size = 14):
    """
    Plot a :class:`DiGraph` from networkx.

    :param digraph: a :class:`DiGraph`
    :param format: the output format
    :param layout: the layout algorithm (from graphviz)
    :param nodelabel_func:
    :param tasknode_func:
    :param assetnode_func:
    :param assetseqnode_func:
    :param graph_size: attribute "size" for the graph (e.g., "3,2")
    :param font_size: font size
    """
    if not has_pygraphviz:
        raise RuntimeError('The Python package "pygraphviz" is required but could not be imported.')
    ag = pygraphviz.AGraph(directed=True)
    if graph_size is not None:
        ag.graph_attr.update(size = graph_size)
    ag.graph_attr.update(font_size = font_size)
    for node in networkx.topological_sort(digraph):
        instance = digraph.node[node]['instance']
        nodelabel = nodelabel_func(instance)
        for node2 in digraph.successors(node):
            instance2 = digraph.node[node2]['instance']
            if _get_pred_done(node2, digraph):
                # if all predecessors are done, default edge display
                ag.add_edge(nodelabel, nodelabel_func(instance2))
            else:
                # if not all predecessors done, lighter edge
                ag.add_edge(nodelabel, nodelabel_func(instance2), color="grey")
        if isinstance(instance, easy.Task):
            tasknode_func(ag, instance, nodelabel)
            pred_done = instance.status == hortator._TASK_DONE
        else:
            pred_done = _get_pred_done(node, digraph)
            if isinstance(instance, hortator.StoredSequence):
                assetseqnode_func(ag, instance, nodelabel, font_size, pred_done)
            else:
                assetnode_func(ag, instance, nodelabel, font_size, pred_done)
    ag.layout(prog=layout)
    return ag


def _plot_digraph(digraph, format, graph_size=None):
    ag = agraph_fromdigraph(digraph, graph_size=graph_size)
    string = ag.draw(format=format)
    return string

def svg_digraph(digraph, graph_size=None):
    return _plot_digraph(digraph, 'svg', graph_size=graph_size)

def svg_taskgraph_view(taskgraph, graph_size=None):
    svg = svg_digraph(taskgraph.digraph)
    return svg

def agraph_fromtasksetgraph(tasksetgraph, taskgraph,
                            graph_size=None,
                            nodelabel_func = node_label, layout="dot"):
    ag = agraph_fromdigraph(taskgraph.digraph, graph_size=graph_size, layout=layout)
    el = tasksetgraph.execution_list()
    knownlabels = dict()
    for tasksetbundle in el:
        if len(tasksetbundle.taskset) == 1:
            continue
        nodebundle = list()
        label = tasksetbundle.label
        if label is None:
            a_set = set()
            for task in tasksetbundle.taskset:
                a_set.update(set(a.value for a in task.activities)),
            label = '/'.join(a_set)

        if label not in knownlabels:
            knownlabels[label] = 1
        else:
            knownlabels[label] += 1
            label = "%s - %i" % (label, knownlabels[label])

        for task in tasksetbundle.taskset:
            #node = easy.taskgraph.Node('Task', task.task_id)
            nodelabel = nodelabel_func(task)
            nodebundle.append(nodelabel)

        ag.add_subgraph(nodebundle, name='cluster_' + label, label=label, color="grey")
    # likely bug in pygraphviz - this is needed
    ag.layout(prog=layout)
    return ag
        
def svg_tasksetgraph_view(tasksetgraph, taskgraph,
                          graph_size=None,
                          nodelabel_func = node_label, layout="dot"):
    ag = agraph_fromtasksetgraph(tasksetgraph, taskgraph,
                                 graph_size=graph_size,
                                 nodelabel_func = node_label, layout="dot")
    svg = ag.draw(format='svg')
    return svg


def _project_status_curses(stdscr, project, checkinterval=3):
    import time, curses, datetime
    time_lastupdate = os.path.getmtime(project._db_fn)
    while True:
        pos_x = 0
        pos_y = 0
        stdscr.addstr(pos_x, pos_y, 'Fetching status...', curses.A_REVERSE)
        stdscr.refresh()
        str_storage = easy.str_project_view_storage(project, width=15)
        str_activities = easy.str_project_view_activities(project, width=15)
        str_results = easy.str_project_view_results(project, width=15)        
        for row in str_storage.split(os.linesep):
            stdscr.addstr(pos_x, pos_y, row)
            pos_x += 1
        stdscr.addstr(pos_x, pos_y, '    Activities    ', curses.A_REVERSE)
        pos_x += 1
        for row in str_activities.split(os.linesep):
            stdscr.addstr(pos_x, pos_y, row)
            pos_x += 1
        stdscr.addstr(pos_x, pos_y, '    Results    ', curses.A_REVERSE)
        pos_x += 1
        for row in str_results.split(os.linesep):
            stdscr.addstr(pos_x, pos_y, row)
            pos_x += 1

        pos_x += 2
        pos_update_x = pos_x
        tdelta = datetime.timedelta(0, 
                                    int(time.time()-time_lastupdate), # seconds 
                                    0, #microseconds 
                                    0, #milliseconds
                                    0)
        stdscr.addstr(pos_update_x, pos_y, '    [Last status change was %s ago]  ' % tdelta,
                      curses.A_DIM)
        pos_x += 1
        stdscr.addstr(pos_x, pos_y, 
                      '    [Checking updates every %.1f seconds]' % (checkinterval), 
                      curses.A_DIM)
        pos_x += 1
        stdscr.addstr(pos_x, pos_y, '    [Ctrl-c to exit]  ',
                      curses.A_DIM)
        stdscr.refresh()
        time_lastupdate = os.path.getmtime(project._db_fn)
        while time_lastupdate == os.path.getmtime(project._db_fn):
            time.sleep(checkinterval)
            tdelta = datetime.timedelta(0, 
                                        int(time.time()-time_lastupdate), # seconds 
                                        0, #microseconds 
                                        0, #milliseconds
                                        0)
            stdscr.addstr(pos_update_x, pos_y, '    [Last status change was %s ago]  ' % tdelta,
                          curses.A_DIM)
            stdscr.refresh()

def _project_status_str(args):
    import importlib
    model = importlib.import_module(args.model)
    project = easy.Project(model, args.wd, db_fn=args.db_fn)
    if args.follow:
        import curses
        import sys
        try:
            curses.wrapper(_project_status_curses, project, checkinterval=args.checkinterval)
        except KeyboardInterrupt:
            sys.exit(1)
    else:
        print(easy.str_project_view(project))


# Custom display
def init_printing():
    from railroadtracks import easy, core
    ip = get_ipython()
    html_f = ip.display_formatter.formatters['text/html']
    html_f.for_type(easy.Project, html_project_view)
    html_f.for_type(core.AssetSet, html_assetset_view)
    html_f.for_type(easy.Task, html_task_view)
    svg_f = ip.display_formatter.formatters['image/svg+xml']
    svg_f.for_type(easy.taskgraph.TaskGraph, svg_taskgraph_view)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='rsrs project tools')
    subparsers = parser.add_subparsers(help='Action to perform. `status` will report the status of the project.')

    parser_status = subparsers.add_parser('status', help="Report status")

    parser_status.add_argument('wd',
                               help='Working directory')
    parser_status.add_argument('model',
                               help='Model (as a Python module)')
    parser_status.add_argument('--db-fn', dest='db_fn',
                               default = None,
                               help='File name for the database')
    parser_status.add_argument('-f', '--follow', dest='follow',
                               action='store_true',
                               help='Follow the status')
    parser_status.add_argument('-t', '--time',
                               dest='checkinterval',
                               type=float,
                               default=3,
                               help='Time interval to update status (in seconds)')
    parser_status.set_defaults(func = _project_status_str)
    
    args = parser.parse_args()
    args.func(args)
