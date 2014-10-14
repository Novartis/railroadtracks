# Copyright 2014 Novartis Institutes for Biomedical Research

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
from railroadtracks import core
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



from railroadtracks.hortator import GRAPHDISPLAY
    
def _svg_graph(project, stored_entities,
               func,
               display=GRAPHDISPLAY.STEPLEVEL):
    import networkx as nx
    dag_list = list()
    for stored_entity in stored_entities:
        dag = func(stored_entity,
                   display=display)
        dag_list.append(dag)
    final_dag = nx.compose_all(dag_list)
    import tempfile
    with tempfile.NamedTemporaryFile(suffix = '.svg') as fh_out:
        adag = nx.to_agraph(final_dag)
        adag.layout(prog=display.value['layout'],
                    args=display.value['layoutargs'])
        adag.draw(path = fh_out.name)
        svg = fh_out.read()
    return svg

def svg_provenancegraph(project, stored_entity, display=GRAPHDISPLAY.STEPLEVEL):
    return _svg_graph(project, stored_entity,
                      project.todo.provenancegraph_storedentity,
                      display = display)

def svg_destinationgraph_stepconcrete(project, step_concrete, display=GRAPHDISPLAY.STEPLEVEL):
    return _svg_graph(project, step_concrete,
                      project.todo.destinationgraph_stepconcrete,
                      display = display)

def svg_destinationgraph_storedentity(project, stored_entity, display=GRAPHDISPLAY.STEPLEVEL):
    return _svg_graph(project, stored_entity,
                      project.todo.destinationgraph_storedentity,
                      display = display)

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
    import railroadtracks.easy as easy
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
