## Creating a project

Railroadtracks is a toolkit designed to help the handling of computational tasks
when:
- the tasks have a relatively high rate of failure because running the underlying software
  is brittle (resources such as RAM or time not known beforehand, spurious crash)
- the tasks can be iteratively added (the results of tasks can decide what is next)
- the programmatic construction of pipelines can be desirable

This ipython notebook is a self-contained tutorial to show some its features.

With railroadtracks, projects are typically contained in a directory on the file system: this
will contain a database describing the dependency graph as well as the results for all tasks.

We will use the following directory for this tutorial (it can be changed to a different location
if needed):

```python
tmpdir = '/tmp/rrt_tutorial' # this directory will be erased !
```

The "!" will run a shell command. With this we make sure that we can start cleanly.

```python
! rm -rf $tmpdir
```

In addition to the directory, a project will require a model module. This is a python module in which
classes to represent the tasks are declared or imported from other module. Here we use the model
for RNA-Seq included with `railroadtracks`.

```python
# import the model for RNA-Seq
from railroadtracks import rnaseq
```

The package is designed with interactive (REPL) programming in mind.
Error messages should be relatively explicit, and lead the user to a quick resolution of problems.

For example, trying to create a project now will end with an error.

```python
import os
from railroadtracks import easy

model = rnaseq
# try/except to allow a fully automated evaluation of the notebook
try:
    project = easy.Project(model, 
                           wd = os.path.join(tmpdir, 'project'))
except Exception as e:
    print('Error: %s' % str(e))
```

The issue can be fixed by creating the project directory.

```python
! mkdir -p $tmpdir/project
```

This time the creation of a project is working.

```python
project = easy.Project(model,
                       wd = os.path.join(tmpdir, 'project'))
```

## Displaying a project

Railroadtracks objects often have a custom string representation. For example, we can print our project:

```python
print(project)
```

HTML display is also available for some of the objects, and this can be used in the IPython notebook.

```python
from IPython.display import display
import railroadtracks.ipython as rsip
rsip.init_printing()
```

Our project can now be "displayed" in the notebook.

```python
display(project)
```

## Working with model environments

A model environment is a convenience class that exposes the content of a model (see the creation of a project above) in a form convenient for writing recipes (that is sequence of steps), with a special attention to interactive use (for example with ipython).

Creating a working environment from a model is meant to be as easy as creating a project.

<div>Note

The attentive reader will note that an other model could be written and used (either as a module like `rnaseq` currently is, or as an object)

</div>


```python
env = easy.Environment(rnaseq)
```

In this example, we will build
- a bowtie index for a reference genome
- align reads in FASTQ files using that index

The package is shipping with the public sequence of a phage, and will use this as a conveniently small
toy genome.

```python
import railroadtracks.model.simulate
reference_fn = railroadtracks.model.simulate.PHAGEFASTA
```

We also need an indexer (indexing part of an aligner).

The working environment we created contains an attribute activities that acts a namespace of activiyies declared in the model. Writing env.activities and hitting `<tab>` on ipython would list the options available.

Here we select bowtie2.

```python
bowtie2index = env.activities.INDEX.bowtie2build
```


## Aligning reads using bowtie2

The object `bowtie2index` is a "step" that can perform the indexing of reference genome.
The embodiment of that step into a task, that is the application of the step to a specific
set of input file requires the definition of the associated assets.

The assets are basically the input files (source assets) and the output files (target assets).

Each step contains a nested class `Assets`,
itself with a nested class `Source` and a nested
class `Target`. These should be used to construct our assets.

```python
try:
    source = bowtie2index.Assets.Source(reference_fn)
except Exception as e:
    print('Error: %s' % str(e))
```

Again, error messages are meant to be rather explicit. Here the model is typing the input assets (which allows early consistency checks). The error disappears when the required class is used.

```python
reference_fasta = rnaseq.FASTAFile(reference_fn)
source = bowtie2index.Assets.Source(reference_fasta)
display(source)
```

When thinking about a task, we let one focus on the input for the task while deferring
the definition of the output (where the files will be stored) to a later stage
(and even have it handled for you).

Target assets can be created as "undefined".

```python
target = bowtie2index.Assets.Target.createundefined()

display(target)
```

Building the complete set of assets is done by combining the source and target assets.

```python
assets_index = bowtie2index.Assets(source, target)
```

Whenever the target assets are undefined, one can implicitly mean it by using Python's `None` or
omitting the target altogether:

```python
assets_index = bowtie2index.Assets(source, target=None)
# or
assets_index = bowtie2index.Assets(source)
```

With both a step and assets, we can now create a task:

Note that the task is labelled as "To do"...

```python
task = project.add_task(bowtie2index, assets_index)
display(task)
```

```python
task.info
```

```python
task.task_id
```

Note that the target is now fully defined. This was performed automatically by the framework.

```python
display(assets_index)
```

Thanks to the persistent layer, adding tasks can be done without having to explicitly
keep track of the operations already performed. Identical tasks will be identified as
such and a unique ID is assigned to them.

To demonstrate this, we add the same task again (with an undefined target).
The result is the same ID (1), and the attribute new tells that this ID is not new.

```python
assets_again = bowtie2index.Assets(bowtie2index.Assets.Source(rnaseq.FASTAFile(reference_fn)))
task = project.add_task(bowtie2index, assets_again)
display(task)

```

Since our model provides a relatively unified interface for steps, building
the indexes for several tools can be done easily in a loop:

```python
listofindexers = (env.activities.INDEX.bwaindex, 
                  env.activities.INDEX.bowtie2build)
# list of tasks
listoftasks = list()
for indexer in listofindexers:
    # create assets for the indexer
    assets_index = indexer.Assets(indexer.Assets.Source(reference_fasta))
    # add to the list of tasks
    listoftasks.append(project.add_task(indexer, assets_index))
```

Note that our little project is now growing:
```python
display(project)
```


## Task execution

Executing tasks can be done manually with:
```python
for task in listoftasks:
    try:
        task.execute()
        # set the status to "done"
        task.status = easy._TASK_DONE
    except:
        # set the status to "failed"
        task.status = easy._TASK_FAILED
display(project)
```

```python
for task in listoftasks:
    task.reset()
display(project)
```

```python
ts_index = easy.TaskSet()
for task in listoftasks:
    ts_index.add(task)
```

However, railroadtracks is also providing a level of abstraction to execute sets
of tasks, separating the declaration of tasks from where they will be running.

Task mappers are objects that implement a method `map` that can take a `TaskSet`
and run all the tasks it contains.

One such mapper provided with `railroadtracks` is a naive mapper that executes tasks
iteratively and in the current process.

```python
engine = easy.IterativeExecution()
engine.map(ts_index)
```

An other mapper provided is using multiple processes:

```python
for task in ts_index:
    task.reset()
# multiprocessing with 2 parallel processes
NPROCS = 2
engine = easy.MultiprocessingExecution(NPROCS)
engine.map(ts_index)
```

We are also providing a mapper for SGE's qsub and have implemented an in-house one for an
Hadoop cluster. The key point is that the description of the workflow is independent of the
type of infrastructure the tasks will be executed on.


## The graph in the dependency graph

Our dependency graph can be turned into a `networkx` `DiGraph` object.

```python
from railroadtracks.easy.taskgraph import TaskGraph

taskgraph = TaskGraph(project)

```

We are also providing a utility to display the dependency graph.

```python
display(taskgraph)
```

Now let's add more tasks: the aligment of FASTQ files and the counts for the annotated features.

We assume that the file names for the sequencing data are contained in a table:

```python
import pandas

dataf_fq = pandas.DataFrame({'sample_n': tuple('sample_%i' % i for i in range(1, 4)),
                             'read1_fn': tuple('sample_%i_1.fq' % i for i in range(1, 4)),
                             'read2_fn': tuple('sample_%i_2.fq' % i for i in range(1, 4))})
```

```python

import os, tempfile
tmpdir = tempfile.mkdtemp()
with open(reference_fn) as fasta_fh:
    reference = next(railroadtracks.model.simulate.readfasta_iter(fasta_fh))
for idx, row in dataf_fq.iterrows():
    with open(os.path.join(tmpdir, row['read1_fn']), 'w') as read1_fh, \
        open(os.path.join(tmpdir, row['read2_fn']), 'w') as read2_fh:
        railroadtracks.model.simulate.randomPEreads(read1_fh, read2_fh, reference)
```

Adding aligment tasks is then as simple as looping through them:

```python
task_index = project.add_task(bowtie2index,
                              bowtie2index.Assets(bowtie2index.Assets.Source(reference_fasta)))
bowtie2 = env.activities.ALIGN.bowtie2
Assets = bowtie2.Assets
FASTQ = rnaseq.FASTQPossiblyGzipCompressed
ts_align = easy.TaskSet()
for idx, row in dataf_fq.iterrows():
    assets_align = Assets(Assets.Source(task_index.assets.target.indexfilepattern,
                                        FASTQ(os.path.join(tmpdir, row['read1_fn'])),
                                        FASTQ(os.path.join(tmpdir, row['read2_fn']))))
    task = project.add_task(bowtie2, assets_align)
    ts_align.add(task)

```

Adding the quantification tasks is equally easy:

```python
htseqcount = env.activities.QUANTIFY.htseqcount
Assets = htseqcount.Assets
ts_quantify = easy.TaskSet()
for task_align in ts_align:
    assets_quantify = Assets(Assets.Source(task_align.assets.target.alignment,
                                           rnaseq.GFFFile(railroadtracks.model.simulate.PHAGEGFF)))
    task_quantify = project.add_task(htseqcount, assets_quantify,
                                     parameters = tuple(['--will-fail']))
    ts_quantify.add(task_quantify)
```

```python
columnmerger = env.activities.UTILITY.columnmerger
Assets = columnmerger.Assets
filestomerge = tuple(rnaseq.CSVFile(x.assets.target.counts.name) for x in ts_quantify)
assets = Assets(Assets.Source(rnaseq.CSVFileSequence(filestomerge)))
task_merge_into_table = project.add_task(columnmerger, assets, parameters=(0,1))
```


```python
taskgraph = TaskGraph(project)
from IPython.display import SVG
display(SVG(rsip.svg_digraph(taskgraph.digraph, graph_size="8,6")))
```

Having declared tasks and grouped them in sets of tasks that can be parallelized is all that is
needed:


```python
from railroadtracks.easy.tasksetgraph import TaskSetGraph
engine = easy.IterativeExecution()
tsg = TaskSetGraph(project=project, defaultmapper=engine)
# add them in reverse order to show that the execution plan
# does not depend on the order they are entered.
tsg.add(ts_quantify)
tsg.add(ts_align)
# individual tasks can also be added
tsg.add(task_merge_into_table)

# execute all tasks, parallelized with task sets and using
# the default mapper (execution engine)
tsg.execute()
```

Unfortunately, while the alignment tasks were successful the quantification
tasks were not:
```python
display(project)
```

We can visualize this on a graph:
```python
taskgraph = TaskGraph(project)
display(SVG(rsip.svg_digraph(taskgraph.digraph,
                             graph_size="8,6")))
```

Investigating the source of error is straightforward. Here the all quantification
have failed, so we just pick one in the set and try to execute it.
```python
task = next(iter(ts_quantify))
# we could also look at the task number (ID) on the graph and do:
# taskid = 10
# task = project.get_task(task_id)

try:
    task.execute()
except Exception as e:
    print(e)
```


The cause of the failure is an invalid parameter, which we correct in a new task set.

```python
htseqcount = env.activities.QUANTIFY.htseqcount
Assets = htseqcount.Assets
ts_quantify = easy.TaskSet()
for task_align in ts_align:
    assets_quantify = Assets(Assets.Source(task_align.assets.target.alignment,
                                           rnaseq.GFFFile(railroadtracks.model.simulate.PHAGEGFF)))
    task_quantify = project.add_task(htseqcount, assets_quantify,
                                     parameters = htseqcount._noexons_parameters)
    ts_quantify.add(task_quantify)
Assets = columnmerger.Assets
filestomerge = tuple(rnaseq.CSVFile(x.assets.target.counts.name) for x in ts_quantify)
assets = Assets(Assets.Source(rnaseq.CSVFileSequence(filestomerge)))
task_merge_into_table = project.add_task(columnmerger, assets, parameters=(0,1))
```

The task set is then added to our task set graph.

```python
tsg.add(ts_quantify)
tsg.add(task_merge_into_table)
tsg.execute()
```

The new parameters allow us to succesfully run all tasks and obtain a table of counts.

This time we are using the additional grouping of tasks in the `TaskSetGraph` tsg when
displaying the dependency graph:

```python
taskgraph = TaskGraph(project)
display(SVG(rsip.svg_tasksetgraph_view(tsg, taskgraph, graph_size="8,6")))
```

More task can be added, such as for example extracting unaligned reads from the BAM files and
producing FASTQ files containing them.

```python
SamtoolsExtractUnaligned = railroadtracks.model.aligners.SamtoolsExtractUnaligned
taskset_unaligned = easy.TaskSet()

bamtofastq = railroadtracks.model.files.BedtoolsBamToFastqPE()

taskset_backtofq = easy.TaskSet()

for task_align in ts_align:
    # extract unaligned reads
    Assets = SamtoolsExtractUnaligned.Assets
    assets = Assets(Assets.Source(task_align.assets.target.alignment))
    task_extract = project.add_task(SamtoolsExtractUnaligned(), assets)
    taskset_unaligned.add(task_extract)

    # convert to FASTQ
    Assets = bamtofastq.Assets
    assets = Assets(Assets.Source(task_extract.assets.target.unaligned))
    taskset_backtofq.add(project.add_task(bamtofastq, assets))

tsg.add(taskset_unaligned)
tsg.add(taskset_backtofq)
taskgraph = TaskGraph(project)
display(SVG(rsip.svg_tasksetgraph_view(tsg, taskgraph, graph_size="8,6")))
```

## Building a full pipeline

Adding additional alternative processing is then only a matter of nesting loops.
If we want to expand on our initial consideration to use BWA, say for comparison purposes,
and build the same tasks for it is done in very few additional lines of code.

### Building blocks for the pipeline

Four functions can be written.

#### 1/4: Align the reads

```python

def align_reads(task_index, aligner, dataf_fq, project):
    # align the reads
    taskset_align = easy.TaskSet()
    Assets = aligner.Assets
    FASTQ = rnaseq.FASTQPossiblyGzipCompressed
    for idx, row in dataf_fq.iterrows():
        assets_align = Assets(Assets.Source(task_index.assets.target.indexfilepattern,
                                            FASTQ(os.path.join(tmpdir, row['read1_fn'])),
                                            FASTQ(os.path.join(tmpdir, row['read2_fn']))))
        task = project.add_task(aligner, assets_align)
        taskset_align.add(task)
    return taskset_align
```

#### 2/4: Quantify the alignments

```python
def quantify_reads(taskset_align, project, additional_parameters = ()):
    # quantify the alignments
    htseqcount = env.activities.QUANTIFY.htseqcount
    parameters = list(htseqcount._noexons_parameters)
    parameters.extend(additional_parameters)
    Assets = htseqcount.Assets
    taskset_quantify = easy.TaskSet()
    for task_align in taskset_align:
        assets_quantify = Assets(Assets.Source(task_align.assets.target.alignment,
                                               rnaseq.GFFFile(railroadtracks.model.simulate.PHAGEGFF)))
        task_quantify = project.add_task(htseqcount, assets_quantify,
                                         parameters = parameters)
        taskset_quantify.add(task_quantify)
    return taskset_quantify
```

#### 3/4: Build a table of counts

```python
def build_table_of_counts(taskset_quantify, project):
    # merge the alignments into a table of counts
    Assets = columnmerger.Assets
    filestomerge = tuple(rnaseq.CSVFile(x.assets.target.counts.name) for x in taskset_quantify)
    assets = Assets(Assets.Source(rnaseq.CSVFileSequence(filestomerge)))
    task_merge = project.add_task(columnmerger, assets, parameters=(0,1))
    return task_merge
```

#### 4/4: Save unaligned reads into FASTQ files

```python
def unaligned_to_FASTQ(taskset_align, project):
    SamtoolsExtractUnaligned = railroadtracks.model.aligners.SamtoolsExtractUnaligned
    bamtofastq = railroadtracks.model.files.BedtoolsBamToFastqPE()
    taskset_unaligned = easy.TaskSet()
    taskset_backtofq = easy.TaskSet()
    for task_align in taskset_align:
        # extract unaligned reads
        Assets = SamtoolsExtractUnaligned.Assets
        assets = Assets(Assets.Source(task_align.assets.target.alignment))
        task_extract = project.add_task(SamtoolsExtractUnaligned(), assets)
        taskset_unaligned.add(task_extract)

        # convert to FASTQ
        Assets = bamtofastq.Assets
        assets = Assets(Assets.Source(task_extract.assets.target.unaligned))
        taskset_backtofq.add(project.add_task(bamtofastq, assets))
    return (taskset_unaligned, taskset_backtofq)
```

### Assembly of the pipeline

With the four short functions defined above, a workflow that aligns reads with either bowtie2 or BWA and
for each aligner saves the unaligned reads into FASTQ files and quantifies the aligments
can be written in very few lines of code.

One can note that the workflow is also quite generic. It is
only depending on `reference_fasta`, the FASTA file
with the reference genome, and on `dataf_fq`, a pandas `DataFrame` that contains
the file names for the FASTQ files.

```python
os.mkdir(os.path.join(tmpdir, 'project2'))
project2 = easy.Project(rnaseq,
                        wd = os.path.join(tmpdir, 'project2'))

# pairs of indexer/aligner
listofindexers = ((env.activities.INDEX.bwaindex, env.activities.ALIGN.bwa),
                   (env.activities.INDEX.bowtie2build, env.activities.ALIGN.bowtie2))

# additional parameters for the quantifier
quant_params = (('-m', 'union'),
                ('-m', 'intersection-strict'),
                ('-m', 'intersection-nonempty'))

tsg2 = TaskSetGraph(defaultmapper=engine)

taskset_index = easy.TaskSet()
# loop through the indexers
for indexer, aligner in listofindexers:

    # create assets for the indexer
    assets_index = indexer.Assets(indexer.Assets.Source(reference_fasta))
    # add to the list of tasks
    task_index = project2.add_task(indexer, assets_index)
    taskset_index.add(task_index)

    # align the reads
    taskset_align = align_reads(task_index, aligner, dataf_fq, project2)
    tsg2.add(taskset_align)

    for additional_parameters in quant_params:
        # quantify the alignments
        taskset_quantify = quantify_reads(taskset_align, project2,
                                          additional_parameters=additional_parameters)
        tsg2.add(taskset_quantify)
        # merge the counts into a table
        task_merge = build_table_of_counts(taskset_quantify, project2)
        tsg2.add(task_merge)

	
    # extract the unaligned read into FASTQ files
    taskset_unaligned, taskset_backtofq = unaligned_to_FASTQ(taskset_align, project2)	
    tsg2.add(taskset_unaligned)
    tsg2.add(taskset_backtofq)

# add the set of indexing tasks
# (order of addition to the TaskSetGraph does not matter)
tsg2.add(taskset_index)
```

Running the pipeline is done as simply as before:

```python
tsg2.execute()
```

Visualizing the pipeline, together with the grouping of tasks for parallel execution
is also done as for simpler pipelines:

```python
taskgraph = TaskGraph(project2)
display(SVG(rsip.svg_tasksetgraph_view(tsg2, taskgraph, graph_size="10,8")))
```




## Handling results (target assets)

### Results from one task

Let's assume that the table of count is our primary interest.

The `Task` "task_merge_into_table" we created earlier for our first project is all we need.

```python
dataf_counts = pandas.read_csv(task_merge_into_table.assets.target.counts.name, engine="python")
display(dataf_counts)
```

It is also possible to look for tasks with the help of helping functions

```python
all_merge = project.get_tasksoftype(rnaseq.ColumnMerger)
all_success_merge = all_merge.filter_on_status(easy.TASK_DONE)
# one task left. That's our task_merge
for task in all_success_merge:
    display(task)
```

Pre-analysis settings such as fixing sample name is, again, straightforward
(the order of is the one in the data frame)

```python
new_names = dict(zip(dataf_counts.columns.values[1:],
                     tuple(dataf_fq['sample_n'])))
dataf_counts.rename(columns = new_names)
dataf_counts.columns.values[1:] = dataf_fq['sample_n']
display(dataf_counts)
```


Flattening the provenance is relatively easy with this example.
This can be helpful for building a table with information about how any given result
was obtained.

```python
# starting with the result of merging the result of quantification into a table of counts
result = tuple(task_merge_into_table.iter_targetassets())[0]
tg = TaskGraph(project)
tasklist = tg.digraph.predecessors(easy.taskgraph.Node('Asset', result.id.id))
taskset = easy.TaskSet()
while len(tasklist) > 0:
    t = tasklist.pop()
    tasklist.extend(tt for tt in tg.predecessor_tasks(t))
    taskset.add(tg.digraph.node[t]['instance'])

set(task.call.step._name for task in taskset)
```

If individual information about the intermediate files is wanted, flattening the graph
into a table could be done with the following function:

```python
def provenance_table(taskgraph, task_merge):
    result = tuple(task_merge.iter_sourceassets())[0]
    node = taskgraph.digraph.node[easy.taskgraph.Node('AssetSequence', result.id.id)]
    assetsequence_ids = tuple(node['instance'])
    ld = list()
    for idx, asset_id in enumerate(assetsequence_ids):
        tasklist = taskgraph.digraph.predecessors(easy.taskgraph.Node('Asset', asset_id))
        taskset = easy.TaskSet()
        d = dict()
        d['idx'] = idx
        while len(tasklist) > 0:
            t = tasklist.pop()
            tasklist.extend(tt for tt in taskgraph.predecessor_tasks(t))
            this_task = taskgraph.digraph.node[t]['instance']
            taskset.add(this_task)
        for task in taskset:
            label = ' '.join(x.value for x in task.call.step.activities)
            d[label + '-taskid'] = task.task_id
            d[label] = task.call.step._name
            d[label + '-params'] = ' '.join(task.call.parameters)
        ld.append(d)
        
    dataf_provenance = pandas.DataFrame.from_records(ld,)
    return dataf_provenance
```


```python
tg = TaskGraph(project)
dataf_provenance = provenance_table(tg, task_merge_into_table)
display(dataf_provenance)
```

Since the order of sample used from the begining is the one in `dataf_fq`, making
a table with sample data and provenance is straightforward:

```python

display(dataf_provenance.join(dataf_fq))

```

### Handling several alternative results

Extracting results from alternative processing routes integrated into one pipeline is also quite easy.

First we execute the tasks in our second project (in case we have not done it before - executing
the tasks will not be attempted if previously successful).

```
tsg2.execute()
```

Build a data frame that contains the results from all alternatives is straightforward:

```python
res = list()
for task in project2.get_tasksoftype(rnaseq.ColumnMerger):
    dataf_counts = pandas.read_csv(task.assets.target.counts.name, engine="python")
	# add the task ID that produced that table. the task ID can be used to recover anything about the
	# task or its provenance.
    dataf_counts['task_id'] = pandas.Series(task.task_id, index=dataf_counts.index)
    res.append(dataf_counts)
dataf_counts = pandas.concat(res)
```

Here we observe the difference on the number of reads not mapping to any feature across the alternative
pipelines:

```python
criterion = dataf_counts['ID'].map(lambda x: x == '__no_feature')
display(dataf_counts.loc[criterion])
```

Linking a task ID to its provenance seen above in the example with one alternative would be similar in
this case.


```python
tg = TaskGraph(project2)
res = list()
for task in project2.get_tasksoftype(rnaseq.ColumnMerger):
    dataf_provenance = provenance_table(tg, task)
    dataf_provenance['task_id'] = pandas.Series(task.task_id, index=dataf_provenance.index)
    res.append(dataf_provenance)
dataf_provenance = pandas.concat(res)

dataf_provenance = dataf_provenance[['Align', 'Quantify', "Quantify-params", 'task_id']]
dataf_provenance = dataf_provenance.drop_duplicates()
dataf = pandas.merge(dataf_counts.loc[criterion], dataf_provenance, on='task_id', how='inner')
display(dataf)
```
