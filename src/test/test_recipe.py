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


import unittest, tempfile, csv, shutil, logging

from railroadtracks import core, hortator, rnaseq, easy, environment

from railroadtracks.model.simulate import (PHAGEGFF,
                                           PHAGEGTF,
                                           PHAGEFASTA)
import railroadtracks.model.simulate

from railroadtracks.test.test_model import (_build_UpToAlign,
                                            _build_StepIndex,
                                            _build_StepAlign,
                                            _build_StepQuantify)


# Test the writing of recipes
from railroadtracks import easy

class RecipeTestCase(unittest.TestCase):

    def setUp(self):
        # -- recipe-init-begin
        # -- initialization boiler plate code
        wd = tempfile.mkdtemp()
        project = easy.Project(rnaseq, wd=wd)

        # declare the 3rd-party command-line tools we will use
        env = easy.Environment(rnaseq)
        # -- recipe-init-end

        # -- recipe-data-begin
        # Phage genome shipped with the package for testing purposes

        PHAGEFASTA = railroadtracks.model.simulate.PHAGEFASTA
        PHAGEGFF = railroadtracks.model.simulate.PHAGEGFF

        # create random data for 6 samples (just testing here)
        nsamples = 6
        samplereads = list()
        with open(PHAGEFASTA) as fasta_fh:
            reference = next(railroadtracks.model.simulate.readfasta_iter(fasta_fh))
        for sample_i in range(nsamples):
            read1_fh = tempfile.NamedTemporaryFile(prefix='read1', suffix='.fq')
            read2_fh = tempfile.NamedTemporaryFile(prefix='read2', suffix='.fq')
            read1_fh, read2_fh = railroadtracks.model.simulate.randomPEreads(read1_fh,
                                                                             read2_fh, 
                                                                             reference)
            samplereads.append((read1_fh, read2_fh))

        sampleinfo_fh = tempfile.NamedTemporaryFile(suffix='.csv', mode='w+')
        csv_w = csv.writer(sampleinfo_fh)
        csv_w.writerow(['sample_id', 'group'])
        for i in range(6):
            csv_w.writerow([str(i), ('A','B')[i%2]])
        sampleinfo_fh.flush()
        referenceannotation = rnaseq.GFFFile(PHAGEGFF)
        # -- recipe-data-end

        self._wd = wd
        self.project = project
        self.reference_fn = PHAGEFASTA
        self.env = env
        self.nsamples = nsamples
        self.samplereads = samplereads
        self.sampleinfo_fh = sampleinfo_fh
        self.referenceannotation = referenceannotation
        self._PHAGEFASTA = PHAGEFASTA
        self._PHAGEGFF = PHAGEGFF

    def tearDown(self):
        samplereads = self.samplereads
        # -- recipe-teardown-begin
        for read1_fh, read2_fh in self.samplereads:
            read1_fh.close()
            read2_fh.close()
        # FIXME: delete the temporary directory
        shutil.rmtree(self.project.wd)
        # -- recipe-teardown-end

    def test_File(self):
        #FIXME: rather test it in the model ?
        reference = core.File(self.reference_fn)

    @unittest.skipIf(not (environment.Executable.ispresent('bowtie2-build') and \
                          environment.Executable.ispresent('htseq-count') and \
                          environment.Executable.ispresent('R') and \
                          environment.R('R').packageversion_or_none('edgeR') is not None),
                     'bowtie2 and/or htseq-count is not in the PATH')
    def test_RecipeSimpleIncremental(self):
        project = self.project
        env = self.env
        nsamples = self.nsamples
        samplereads = self.samplereads
        sampleinfo_fh = self.sampleinfo_fh
        reference_fn = self.reference_fn
        referenceannotation = self.referenceannotation
        PHAGEFASTA = self._PHAGEFASTA
        PHAGEGFF = self._PHAGEGFF

        # steps used
        bowtie2index = env.activities.INDEX.bowtie2build
        bowtie2align = env.activities.ALIGN.bowtie2
        htseqcount = env.activities.QUANTIFY.htseqcount
        merge = env.activities.UTILITY.columnmerger
        edger = env.activities.DIFFEXP.edger

        from railroadtracks import easy

        # sequence of tasks to run
        torun = list()

        # index for alignment
        Assets = bowtie2index.Assets
        assets = Assets(Assets.Source(rnaseq.FASTAFile(reference_fn)),
                        Assets.Target.createundefined())
        task_index = project.add_task(bowtie2index, 
                                      assets)
        # the step is not done
        self.assertEqual(hortator._TASK_TODO, task_index.info[1])
        torun.append(task_index)
        # run the tasks
        for task in torun:
            # run only if not done
            if task.info[1] != hortator._TASK_DONE:
                task.execute()
                task.status = hortator._TASK_DONE

        self.assertEqual(1, project.persistent_graph.nconcrete_steps)
        # now that the tasks have run let's open the same project
        project_same = easy.Project(project.model, wd=project.wd)

        # index for alignment
        Assets = bowtie2index.Assets
        assets = Assets(Assets.Source(rnaseq.FASTAFile(reference_fn)),
                        Assets.Target.createundefined())
        task_index_same = project_same.add_task(bowtie2index, 
                                                assets)

        self.assertNotEqual(task_index, task_index_same)
        self.assertNotEqual(task_index.call.assets, task_index_same.call.assets)
        self.assertListEqual(list(task_index.call.assets.source.reference), 
                             list(task_index_same.call.assets.source.reference))
        self.assertListEqual(list(task_index.call.assets.target.indexfilepattern), 
                             list(task_index_same.call.assets.target.indexfilepattern))
        self.assertEqual(hortator._TASK_DONE, task_index_same.info[1])
        self.assertEqual(1, project.persistent_graph.nconcrete_steps)


    def _recipesimpleincremental(self, runtasks):
        project = self.project
        env = self.env
        nsamples = self.nsamples
        samplereads = self.samplereads
        sampleinfo_fh = self.sampleinfo_fh
        reference_fn = self.reference_fn
        referenceannotation = self.referenceannotation
        PHAGEFASTA = self._PHAGEFASTA
        PHAGEGFF = self._PHAGEGFF
        
        # steps used
        bowtie2index = env.activities.INDEX.bowtie2build
        bowtie2align = env.activities.ALIGN.bowtie2
        htseqcount = env.activities.QUANTIFY.htseqcount
        merge = env.activities.UTILITY.columnmerger
        edger = env.activities.DIFFEXP.edger

        for iteration in range(5):
            nextiteration = False
            # sequence of tasks to run
            torun = list()

            # index for alignment
            Assets = bowtie2index.Assets
            assets = Assets(Assets.Source(rnaseq.FASTAFile(reference_fn)),
                            Assets.Target.createundefined())
            task_index = project.add_task(bowtie2index, assets)
            torun.append(task_index)
            if iteration < 1:
                nextiteration = True
                runtasks(torun)
                self.assertEqual(1, project.persistent_graph.nconcrete_steps)
                continue
            # process all samples
            sample_counts = list()
            for sample_i, (read1_fh, read2_fh) in enumerate(samplereads):
                # align
                Assets = bowtie2align.Assets
                assets = Assets(Assets.Source(task_index.call.assets.target.indexfilepattern, 
                                              rnaseq.FASTQPossiblyGzipCompressed(read1_fh.name),
                                              rnaseq.FASTQPossiblyGzipCompressed(read2_fh.name)),
                                Assets.Target.createundefined())
                task_align = project.add_task(bowtie2align, assets)
                torun.append(task_align)
                if iteration < 2:
                    nextiteration = True
                    runtasks(torun)
                    self.assertEqual(1+(sample_i+1), project.persistent_graph.nconcrete_steps)
                    continue

                # quantify
                # (non-default parameters to fit our demo GFF)
                params = rnaseq.HTSeqCount._noexons_parameters
                Assets = htseqcount.Assets
                assets = Assets(Assets.Source(task_align.call.assets.target.alignment,
                                              rnaseq.GFFFile(referenceannotation)),
                                Assets.Target.createundefined())
                task_quantify = project.add_task(htseqcount,
                                                 assets,
                                                 parameters=params)
                torun.append(task_quantify)
                if iteration < 3:
                    nextiteration = True
                    runtasks(torun)
                    self.assertEqual(1+len(samplereads)+(sample_i+1), 
                                     project.persistent_graph.nconcrete_steps)
                    continue

                # keep a pointer to the counts, as we will use it in the merge step
                sample_counts.append(task_quantify.call.assets)

            if nextiteration:
                continue
            # merge the sample data into a table (so differential expression can be computed)
            Assets = merge.Assets
            counts = tuple(x.target.counts for x in sample_counts)
            assets = Assets(Assets.Source(rnaseq.CSVFileSequence(counts)),
                            merge.Assets.Target.createundefined())

            task_merge = project.add_task(merge,
                                          assets,
                                          parameters=("0", "1"))
            torun.append(task_merge)
            if iteration < 4:
                nextiteration = True
                runtasks(torun)
                self.assertEqual(1+2*len(samplereads)+1, 
                                 project.persistent_graph.nconcrete_steps)
                continue

            # differential expression with edgeR
            Assets = edger.Assets
            assets = Assets(Assets.Source(task_merge.call.assets.target.counts,
                                          rnaseq.CSVFile(sampleinfo_fh.name)),
                            Assets.Target.createundefined())
            task_de = project.add_task(edger,
                                       assets)
            if iteration < 5:
                nextiteration = True
                runtasks(torun)
                self.assertEqual(1+2*len(samplereads)+2, # 1 index + 2 FASTQ per sample + 1 merge + 1 differential expression
                                 project.persistent_graph.nconcrete_steps)
                continue

    @unittest.skipIf(not (environment.Executable.ispresent('bowtie2-build') and \
                          environment.Executable.ispresent('htseq-count') and \
                          environment.Executable.ispresent('R') and \
                          environment.R('R').packageversion_or_none('edgeR') is not None),
                     'bowtie2 and/or htseq-count is not in the PATH')
    def test_RecipeSimpleIncrementalComplete(self):
        def runtasks(torun):
            # run the tasks
            for task in torun:
                # run only if not done
                if task.info[1] != hortator._TASK_DONE:
                    task.execute()
        self._recipesimpleincremental(runtasks)

    @unittest.skipIf(not (environment.Executable.ispresent('bowtie2-build') and \
                          environment.Executable.ispresent('htseq-count') and \
                          environment.Executable.ispresent('R') and \
                          environment.R('R').packageversion_or_none('edgeR') is not None),
                     'bowtie2, htseq-count, R (with package "edgeR") must be in the PATH')
    def test_RecipeSimpleIncrementalCompleteNoRun(self):
        def runtasks(torun):
            # do nothing
            pass
        self._recipesimpleincremental(runtasks)


    @unittest.skipIf(not (environment.Executable.ispresent('bowtie2-build') and \
                          environment.Executable.ispresent('htseq-count') and \
                          environment.Executable.ispresent('R') and \
                          environment.R('R').packageversion_or_none('edgeR') is not None),
                     'bowtie2, htseq-count, R (with package "edgeR") must be in the PATH')
    def test_RecipeSimple(self):
        project = self.project
        env = self.env
        nsamples = self.nsamples
        samplereads = self.samplereads
        sampleinfo_fh = self.sampleinfo_fh
        reference_fn = self.reference_fn
        referenceannotation = self.referenceannotation
        PHAGEFASTA = self._PHAGEFASTA
        PHAGEGFF = self._PHAGEGFF
        
        # -- recipesimple-test-begin

        # steps used
        bowtie2index = env.activities.INDEX.bowtie2build
        bowtie2align = env.activities.ALIGN.bowtie2
        htseqcount = env.activities.QUANTIFY.htseqcount
        merge = env.activities.UTILITY.columnmerger
        edger = env.activities.DIFFEXP.edger

        from railroadtracks import easy

        # sequence of tasks to run
        torun = list()
                            
        # index for alignment
        Assets = bowtie2index.Assets
        assets = Assets(Assets.Source(rnaseq.FASTAFile(reference_fn)),
                        Assets.Target.createundefined())
        task_index = project.add_task(bowtie2index, assets)
        torun.append(task_index)

        # process all samples
        sample_counts = list()
        for read1_fh, read2_fh in samplereads:
            # align
            Assets = bowtie2align.Assets
            assets = Assets(Assets.Source(task_index.call.assets.target.indexfilepattern, 
                                          rnaseq.FASTQPossiblyGzipCompressed(read1_fh.name),
                                          rnaseq.FASTQPossiblyGzipCompressed(read2_fh.name)),
                            Assets.Target.createundefined())
            task_align = project.add_task(bowtie2align, assets)
            torun.append(task_align)

            # quantify
            # (non-default parameters to fit our demo GFF)
            params = rnaseq.HTSeqCount._noexons_parameters
            Assets = htseqcount.Assets
            assets = Assets(Assets.Source(task_align.call.assets.target.alignment,
                                          rnaseq.GFFFile(referenceannotation)),
                            Assets.Target.createundefined())
            task_quantify = project.add_task(htseqcount,
                                             assets,
                                             parameters=params)
            torun.append(task_quantify)
            # keep a pointer to the counts,
            # as we will use them in the merge step
            sample_counts.append(task_quantify.call.assets)

        # merge the sample data into a table
        # (so differential expression can be computed)
        Assets = merge.Assets
        counts = tuple(x.target.counts for x in sample_counts)
        assets = Assets(Assets.Source(rnaseq.CSVFileSequence(counts)),
                        merge.Assets.Target.createundefined())
        task_merge = project.add_task(merge,
                                      assets,
                                      parameters=("0","1"))
        torun.append(task_merge)

        # differential expression with edgeR
        Assets = edger.Assets
        assets = Assets(Assets.Source(task_merge.call.assets.target.counts,
                                      rnaseq.CSVFile(sampleinfo_fh.name)),
                        Assets.Target.createundefined())
        task_de = project.add_task(edger,
                                   assets)

        # run the tasks
        for task in torun:
            # run only if not done
            if task.info[1] != hortator._TASK_DONE:
                task.execute()

        # get results
        final_storedentities = project.get_targetsofactivity(rnaseq.ACTIVITY.DIFFEXP)

        # get the step that created the results files
        final_steps = list()
        for stored_entity in final_storedentities:
            final_steps.append(project.persistent_graph.get_parenttask_of_storedentity(stored_entity))
        
        # -- recipesimple-test-end
        
        self.assertEqual(1, len(final_storedentities))
        self.assertEqual(core.File.__name__, final_storedentities[0].clsname)
        self.assertEqual('railroadtracks.model.diffexp.EdgeR', final_steps[0].clsname)

        # FIXME: not yet implemented
        # now that we have all steps, we "only" have to run them
        #steps = todo.stepcrawler()
        #for s in steps:
        #    print('%s' % (s.unifiedname))
        #    s.run()

    @unittest.skipIf(not (environment.Executable.ispresent('bowtie2-build') and \
                          environment.Executable.ispresent('bowtie-build') and \
                          environment.Executable.ispresent('STAR') and \
                          environment.Executable.ispresent('R') and \
                          environment.R('R').packageversion_or_none('edgeR') is not None and \
                          environment.R('R').packageversion_or_none('DESeq') is not None and \
                          environment.R('R').packageversion_or_none('DESeq2') is not None and \
                          environment.R('R').packageversion_or_none('limma') is not None),
                     'bowtie2, bowtie, STAR, TopHat2, and R (with packages "edgeR", "DESeq", "DESeq2", "limma") must be in the PATH')
    def test_RecipeLoop(self):
        project = self.project
        env = self.env
        nsamples = self.nsamples
        samplereads = self.samplereads
        sampleinfo_fh = self.sampleinfo_fh
        reference_fn = self.reference_fn
        referenceannotation = self.referenceannotation
        PHAGEFASTA = self._PHAGEFASTA
        PHAGEGFF = self._PHAGEGFF

        # -- recipeloop-test-begin
        from railroadtracks import easy

        torun = list()

        # bowtie
        bowtie1index = env.activities.INDEX.bowtiebuild
        bowtie1align = env.activities.ALIGN.bowtie
        Assets = bowtie1index.Assets
        fa_file = rnaseq.FASTAFile(reference_fn)
        task_index_bowtie1 = project.add_task(bowtie1index, 
                                              Assets(Assets.Source(fa_file),
                                                     None))
        torun.append(task_index_bowtie1)

        # bowtie2
        bowtie2index = env.activities.INDEX.bowtie2build
        bowtie2align = env.activities.ALIGN.bowtie2
        Assets = bowtie2index.Assets
        fa_file = rnaseq.FASTAFile(reference_fn)
        task_index_bowtie2 = project.add_task(bowtie2index,
                                              Assets(Assets.Source(fa_file),
                                                     None))
        torun.append(task_index_bowtie2)

        # STAR
        starindex = env.activities.INDEX.starindex
        staralign = env.activities.ALIGN.staralign
        Assets = starindex.Assets
        fa_file = rnaseq.FASTAFile(reference_fn)
        task_index_star = project.add_task(starindex, 
                                           Assets(Assets.Source(fa_file),
                                                  None))
        torun.append(task_index_star)

        # TopHat2
        # (index from bowtie2 used)
        #tophat2 = env.activities.ALIGN.tophat2

        # featureCount
        featurecount = env.activities.QUANTIFY.featurecount

        # Merge columns (obtained from counting)
        merge = env.activities.UTILITY.columnmerger

        # EdgeR, DESeq, DESeq2, and LIMMA voom
        edger = env.activities.DIFFEXP.edger
        deseq = env.activities.DIFFEXP.deseq
        deseq2 = env.activities.DIFFEXP.deseq2
        voom = env.activities.DIFFEXP.limmavoom
        

        # Now explore the different alignment presets in bowtie2, and vanilla star
        from itertools import cycle
        from collections import namedtuple
        Options = namedtuple('Options', 'aligner assets_index parameters')
        # Try various presets for bowtie2
        bowtie2_parameters = (('--very-fast', ), ('--fast', ), 
                              ('--sensitive', ), ('--very-sensitive', ))
        options = [Options(*x) for x in zip(cycle((bowtie2align,)),
                                            cycle((task_index_bowtie2.call.assets.target,)),
                                            bowtie2_parameters)]

        # add bowtie
        options.append(Options(bowtie1align, task_index_bowtie1.call.assets.target, tuple()))
        # add STAR (vanilla, no specific options beside the size of index k-mers)
        options.append(Options(staralign, 
                               task_index_star.call.assets.target, 
                               ('--genomeChrBinNbits', '12')))
        # add TopHat2
        #options.append(Options(tophat2, task_index_bowtie2.call.assets.target, tuple()))

        # loop over the options
        for option in options:
            sample_counts = list()
            # loop over the samples
            for sample_i in range(nsamples):
                read1_fh, read2_fh = samplereads[sample_i]
                # align
                Assets = option.aligner.Assets
                assets = Assets(Assets.Source(option.assets_index.indexfilepattern,
                                              rnaseq.FASTQPossiblyGzipCompressed(read1_fh.name), 
                                              rnaseq.FASTQPossiblyGzipCompressed(read2_fh.name)),
                                Assets.Target.createundefined())
                task_align = project.add_task(option.aligner,
                                              assets,
                                              parameters=option.parameters)
                torun.append(task_align)

                # quantify
                # (non-default parameters to fit our demo GFF)
                Assets = featurecount.Assets
                assets = Assets(Assets.Source(task_align.call.assets.target.alignment,
                                              rnaseq.GFFFile(referenceannotation)),
                                Assets.Target.createundefined())
                task_quantify = project.add_task(featurecount,
                                                 assets,
                                                 parameters = ('--gtf-featuretype', 'CDS',
                                                               '--gtf-attrtype', 'ID'))
                torun.append(task_quantify)

                # keep a pointer to the counts, as we will use it in the merge step
                sample_counts.append(task_quantify.call.assets)

            # merge the sample data into a table (so differential expression can be computed)
            Assets = merge.Assets
            source = Assets.Source(rnaseq.CSVFileSequence(tuple(x.target.counts\
                                                                for x in sample_counts)))
            assets_merge = Assets(source,
                                  Assets.Target.createundefined())
            task_merge = project.add_task(merge,
                                          assets_merge,
                                          parameters=("0","1"))
            torun.append(task_merge)

            # differential expression with edgeR, deseq2, and voom
            # (deseq is too whimsical for tests)
            for diffexp, params in ((edger, ()),
                                    (deseq, ('--dispersion-fittype=local', )), 
                                    (deseq2, ()),
                                    (voom, ())):
                Assets = diffexp.Assets
                assets = Assets(Assets.Source(task_merge.call.assets.target.counts,
                                              core.File(sampleinfo_fh.name)),
                                Assets.Target.createundefined())
                task_de = project.add_task(diffexp,assets)
                torun.append(task_de)

        # run the tasks
        # (this is an integration test rather than a unit test - the 
        # 3rd-party tools are often brittle and we want to keep the noise level down)
        env_log_level = environment.logger.level
        environment.logger.level = logging.ERROR
        try:
            for task in torun:
                if task.info[1] != hortator._TASK_DONE:
                    try:
                        task.execute()
                        status = easy.hortator._TASK_DONE
                    except:
                        status = easy.hortator._TASK_FAILED
                project.persistent_graph.step_concrete_state(hortator.DbID(task.task_id, False),
                                                             easy.hortator._TASK_STATUS_LIST[status])
        finally:
            environment.logger.level = env_log_level
        # -- recipeloop-test-end


if __name__ == '__main__':
    unittest.main()
