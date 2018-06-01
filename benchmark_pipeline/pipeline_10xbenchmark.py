##############################################################################
#
#   MRC FGU CGAT
#
#   $Id$
#
#   Copyright (C) 2018 Andreas Heger
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
###############################################################################
"""===========================
10X benchmark pipeline
===========================

:Author: Ian Sudbery and Tom Smith
:Release: 0.1
:Date: 7th March 2018
:Tags: Python

Overviews
========

This pipeline accompanies the first Alevin preprint. It compares Alevin performance on 10X datasets with CellRanger, UMI-Tools.

Usage
=====

Generate default configuration files with::

python /path/to/repo/pipeline_10xbenchmark.py config

Run pipeline with:

python /path/to/repo/pipeline_10xbenchmark.py make full


Configuration
-------------
Configuration may be neccesary to run accross a cluster.


Input files
-----------

None

Requirements
------------

Requirements:

* CGATPipelines
* CGAT
* CellRanger == 2.1.1
* UMI-Tools == 0.5.3
* STAR == 2.4.2a
* featureCounts == 1.4.6
* wget
* samtools


Pipeline output
===============

.. Describe output files of the pipeline here

Glossary
========

.. glossary::

To Do:
- ASdd

Code
====

"""
from ruffus import *
from ruffus.combinatorics import *

import sys
import os
import re
import sqlite3
import collections
import itertools
import glob
import shutil
import json
import numpy as np

#try:
#    import CGAT.IOTools as IOTools
#    import CGAT.Experiment as E
#    import CGATPipelines.Pipeline as P
#    import CGATPipelines.PipelineMapping as PipelineMapping
#except:
import CGATCore.IOTools as IOTools
import CGATCore.Experiment as E
import CGATCore.Pipeline as P

import Benchmark
REPO_PATH = os.path.dirname(__file__)
# load options from the config file
PARAMS = P.getParameters(
    ["%s/pipeline.ini" % REPO_PATH,
     "pipeline.ini"])

# ---------------------------------------------------
# Specific pipeline tasks

################################################################################
# Download and prepare input data and reference transcriptomes
################################################################################

TENX2INFO = collections.defaultdict(lambda: collections.defaultdict())
TENX_DATASETS = set()
with IOTools.openFile(os.path.join(REPO_PATH, "sample_info"), "r") as inf:
    header = next(inf)
    for line in inf:
        if len(line.strip()) == 0:
            break
        if line.startswith("#"):
            continue

        sample_name, r_version, n_cells, species, seq_sat, chem = line.strip().split(",")
        #if chem == "v2":  # can't currently handle the v1 data (Where are the UMIs?!)
        TENX_DATASETS.add(sample_name)
        TENX2INFO[sample_name]["version"] = r_version
        TENX2INFO[sample_name]["n_cells"] = n_cells
        TENX2INFO[sample_name]["species"] = species
        TENX2INFO[sample_name]["chem"] = chem

if PARAMS['limit_samples']:
    TENX_DATASETS = [x for x in TENX_DATASETS if x in P.asList(PARAMS['limit_samples'])]

# restrict for testing (hgmm_12K > 75GB fastqs!!)
TENX_DATASETSV2 = [x for x in TENX_DATASETS if TENX2INFO[x]['chem']=="v2"]
TENX_DATASETSV1 = [x for x in TENX_DATASETS if TENX2INFO[x]['chem']=="v1"]

# ---------------------------------------------------
@mkdir('raw', 'raw/10X_fastqs/')
@originate(['raw/10X_fastqs/%s.fastq.1.gz' % x for x in TENX_DATASETSV2])
def download10x(outfile):
    ''' Download the 10X data, untar, and CellRanger do not want
    catted files, so do leave these present..'''

    sample_name = P.snip(os.path.basename(outfile), ".fastq.1.gz")

    ranger_version = TENX2INFO[sample_name]["version"]

    tar_file = "%s_fastqs.tar" % sample_name
    if sample_name in ["t_4k", "t_3k", "neuron_9k"]:
        url = "http://s3-us-west-2.amazonaws.com/10x.files/samples/cell-exp/%s/%s/%s" % (
        ranger_version, sample_name, tar_file)
    else:
        url = "http://cf.10xgenomics.com/samples/cell-exp/%s/%s/%s" % (
        ranger_version, sample_name, tar_file)

    outfile2 = outfile.replace(".fastq.1.gz", ".fastq.2.gz")

    tmp_dir = sample_name + ".dir"
    if not os.path.isdir(tmp_dir):
        os.mkdir(tmp_dir)

    statement = '''
    wget %(url)s;
    tar -xf %(tar_file)s -C %(tmp_dir)s;
    mv %(tmp_dir)s/*/* %(tmp_dir)s;
    cat %(tmp_dir)s/%(sample_name)s_S1_*_R1_001.fastq.gz
    > %(outfile)s;
    cat %(tmp_dir)s/%(sample_name)s_S1_*_R2_001.fastq.gz
    > %(outfile2)s;
    rm -rf %(tar_file)s
    '''
    P.run()


@mkdir('raw', 'raw/10X_fastqs/')
@originate(['raw/10X_fastqs/%s.fastq.1.gz' % x for x in TENX_DATASETSV1])
def download10xV1(outfile):
    ''' Download the 10X data, untar. and CellRanger do not want
    catted files, so leave these present..'''

    sample_name = P.snip(os.path.basename(outfile), ".fastq.1.gz")

    ranger_version = TENX2INFO[sample_name]["version"]

    tar_file = "%s_fastqs.tar" % sample_name

    url = "http://s3-us-west-2.amazonaws.com/10x.files/samples/cell-exp/%s/%s/%s" % (
        ranger_version, sample_name, tar_file)

    outfile2 = outfile.replace(".fastq.1.gz", ".fastq.2.gz")

    tmp_dir = os.path.join("raw/10X_fastqs", sample_name + ".dir")

    if not os.path.isdir(tmp_dir):
        os.mkdir(tmp_dir)

    statement = '''
    wget %(url)s; checkpoint ;
    tar -xf %(tar_file)s -C %(tmp_dir)s;
    checkpoint;
    mv %(tmp_dir)s/*/* %(tmp_dir)s;
    checkpoint;
    cat %(tmp_dir)s/read-I1*.fastq.gz
    > %(outfile)s; checkpoint ;
    cat %(tmp_dir)s/read-RA*.fastq.gz
    > %(outfile2)s; checkpoint ;
    rm -rf %(tar_file)s
    '''
    P.run()


@mkdir('raw', 'raw/10X_fastqs/')
@originate(['raw/10X_fastqs/%s_whitelist.txt' % x for x in TENX_DATASETS])
def download10xWhitelists(outfile):
    ''' Download the whitelists according to 10X'''

    sample_name = P.snip(os.path.basename(outfile), "_whitelist.txt")

    ranger_version = TENX2INFO[sample_name]["version"]

    tar_file = "%s_filtered_gene_bc_matrices.tar.gz" % sample_name

    url = "http://cf.10xgenomics.com/samples/cell-exp/%s/%s/%s" % (
        ranger_version, sample_name, tar_file)

    outfile2 = outfile.replace(".fastq.1.gz", ".fastq.2.gz")

    tmp_dir = P.getTempDir("./")

    statement = '''
    wget %(url)s; checkpoint ;
    tar -xf %(tar_file)s -C %(tmp_dir)s ; checkpoint ;
    cat %(tmp_dir)s/filtered_gene_bc_matrices/*/barcodes.tsv |
    sed 's/-1//g' > %(outfile)s ; checkpoint ;
    rm -rf %(tmp_dir)s %(tar_file)s
    '''
    P.run()




@follows(mkdir("nreads.dir"))
@transform((download10x, download10xV1),
           regex("raw/10X_fastqs/(\S+).fastq.1.gz"),
           r"nreads.dir/\1.nreads")
def countReads(infile, outfile):

    # hacky way to check statement command
    statement = Benchmark.makeCounterStatement(infile, outfile)
    P.run()


# TS: We can't use variable 'outfile' for this task. I think because
# P.run() specifically looks for outfile in the locals and expands the
# glob or something? I haven't actually got to the bottom of it!
@follows(countReads)
@split("raw/10X_fastqs/%s.fastq.1.gz" % PARAMS['subset_sample'],
       "raw/10X_fastqs/%s_subset*.fastq.1.gz" % PARAMS['subset_sample'])
def subsetRange(infile, subset_outfile):
    '''subset sample to 10%-100% depth.

    In order to make the output compatible with cellranger and alevin
    and to make sure we get a properly random set of reads, we
    randomly subsample the individual fastqs and then concatenate
    them.

    Note, this random sampling holds the all fastq files (R1, R2, I1)
    in memory so they can be randomly sorted together. This probably
    isn't the most efficient method'''

    sample_name = P.snip(os.path.basename(infile), ".fastq.1.gz")

    # making a hard assumption about the naming format here:
    # expect e.g:
    # neuron_9k_S1_L008_R2_001.fastq.gz
    # neuron_9k_S1_L008_R1_001.fastq.gz
    # neuron_9k_S1_L008_I1_001.fastq.gz
    # neuron_9k_S1_L007_I1_001.fastq.gz
    # neuron_9k_S1_L007_R1_001.fastq.gz
    # neuron_9k_S1_L007_R2_001.fastq.gz
    fastq_files = glob.glob(sample_name + ".dir/*R1_001.fastq.gz")

    # need to count the reads
    #tmp_nreads = P.getTempFilename(shared=True)
    #statement = Benchmark.makeCounterStatement(fastq_files[0], tmp_nreads)
    #P.run()

    #nreads = None
    #for line in IOTools.openFile(tmp_nreads, "r"):
    #    if not line.startswith("nreads"):
    #        continue
    #    nreads = int(line[:-1].split("\t")[1])
    #    break
    #if not nreads:
    #    raise ValueError("can't read nreads: %s" % tmp_nreads)

    #os.unlink(tmp_nreads)

    subset_depths = list(np.arange(0.1, 1.1, 0.1))
    #limits = [int(nreads / (100.0 / int(depth)))
    #          for depth in subset_depths]

    # build a big paste|sort -R|awk command to randomly subset triplicate fastqs
    for infile_1 in fastq_files:

        infile_2 = infile_1.replace("R1_001.fastq.gz", "R2_001.fastq.gz")
        infile_3 = infile_1.replace("R1_001.fastq.gz", "I1_001.fastq.gz")

        awk_cmd = "{rfloat=rand()}"
        for sub_ix, limit in enumerate(subset_depths):
            sub_outdir = infile.replace(".fastq.1.gz", "_subset%s.dir" % sub_ix)
            if not os.path.isdir(sub_outdir):
                os.mkdir(sub_outdir)

            # TS: needlessly repetitive but couldn't get the behaviour
            # I wanted otherwise. Please feel free to simplify the awk lines below
            awk_cmd += """
            {if (rfloat<=%s) print $1,$4,$7,$10 | "gzip > %s";
             if (rfloat<=%s) print $2,$5,$8,$11 | "gzip > %s";
             if (rfloat<=%s) print $3,$6,$9,$12 | "gzip > %s"}
            """ % (limit, os.path.join(sub_outdir, os.path.basename(infile_1)),
                   limit, os.path.join(sub_outdir, os.path.basename(infile_2)),
                   limit, os.path.join(sub_outdir, os.path.basename(infile_3)))

        statement = """
        paste <(zcat %(infile_1)s) <(zcat %(infile_2)s) <(zcat %(infile_3)s)|
        paste - - - - |
        awk -F'\\t' 'BEGIN{OFS="\\n";} %(awk_cmd)s' """

        P.run()

        outdir = os.path.dirname(infile)
        for sub_ix, limit in enumerate(subset_depths):
            sub_outdir = infile.replace(".fastq.1.gz", "_subset%s.dir" % sub_ix)

            outfile1 = os.path.join(
                outdir, "%s_subset%s.fastq.1.gz" % (sample_name, sub_ix))
            outfile2 = os.path.join(
                outdir, "%s_subset%s.fastq.2.gz" % (sample_name, sub_ix))
            statement = '''
            cat %(sub_outdir)s/*_R1_001.fastq.gz > %(outfile1)s; checkpoint ;
            cat %(sub_outdir)s/*_R2_001.fastq.gz > %(outfile2)s;
            '''

            P.run()

# ---------------------------------------------------
@mkdir("references")
@originate(["references/hg_transcriptome.gtf",
            "references/mm_transcriptome.gtf"])
def DownloadGTFs(outfile):
    ''' Make a merged GTF for human & mouse'''

    species = os.path.basename(outfile).split("_")[0]

    link = PARAMS["%s_geneset_address" % species]

    statement = '''
    wget %(link)s -O - | zcat > %(outfile)s
    '''

    P.run()

@mkdir("references")
@originate("references/hgbasic_transcriptome.gtf")
def DownloadBasicGTF(outfile):
    '''Download the basic geneset for human'''

    species = os.path.basename(outfile).split("_")[0]

    link = PARAMS["%s_geneset_address" % species]

    statement = '''
    wget %(link)s -O - | zcat > %(outfile)s
    '''

    P.run()


# ---------------------------------------------------
@mkdir("references")
@merge(DownloadGTFs,
       "references/hgmm_transcriptome.gtf")
def MakeMergedGTF(infiles, outfile):
    ''' Make a merged GTF for human & mouse'''

    hg_infile = [x for x in infiles if "hg_" in os.path.basename(x)][0]
    mm_infile = [x for x in infiles if "mm_" in os.path.basename(x)][0]

    statement = '''
    awk '$3=="exon"' %(hg_infile)s
    | grep -v '^!'
    | sed  's/^chr/hg_chr/'
    > %(outfile)s ; checkpoint;

    awk '$3=="exon"'  %(mm_infile)s
    | grep -v '^!'
    | sed  's/^chr/mm_chr/'
    >> %(outfile)s
    '''

    P.run()


# ---------------------------------------------------
@mkdir("references")
@originate(["references/hg_transcriptome.fasta",
            "references/mm_transcriptome.fasta"])
def DownloadTranscriptomes(outfile):
    ''' Download the transcriptomes'''

    species = os.path.basename(outfile).split("_")[0]

    link = PARAMS["%s_transcriptome_address" % species]

    statement = '''
    wget %(link)s -O - | zcat > %(outfile)s
    '''

    P.run()


# ---------------------------------------------------
@mkdir("references")
@merge(DownloadTranscriptomes,
       "references/hgmm_transcriptome.fasta")
def MakeMergedTranscriptome(infiles, outfile):
    ''' Make a merged transcriptome fasta for human & mouse'''

    infiles = " ".join(infiles)
    statement = '''
    cat %(infiles)s > %(outfile)s '''

    P.run()

# ---------------------------------------------------


@mkdir("references")
@merge([DownloadBasicGTF, DownloadTranscriptomes],
       "references/hgbasic_transcriptome.fasta")
def MakeBasicTranscriptome(infiles, outfile):
    '''Subset the fasta to only the genes in the basic GTF'''
    
    gtf = infiles[0]
    hg_fasta = [x for x in infiles[1:] if
                os.path.basename(x) == "hg_transcriptome.fasta"][0]
    
    transcripts = set()
    with open(gtf, "r") as inf:
        for line in inf:
            if line.startswith("#"):
                continue
            attributes = line.strip().split("\t")[8]
            if "transcript_id" in attributes:
                attributes = attributes.split("; ")
                attributes = {x.split(" ")[0]:x.split(" ")[1] for x in attributes}
                if "transcript_id" in attributes:
                    transcripts.add(attributes["transcript_id"].replace('"', ''))

    print(len(transcripts))
    
    transcripts_seen=set()
    with open(outfile, "w") as outf:
        write_lines = False
        for line in open(hg_fasta, "r"):
            if line.startswith("#"):
                continue
            if line.startswith(">"):
                entry_id = line.split("|")[0][1:]
                if entry_id in transcripts:
                    transcripts_seen.add(entry_id)
                    write_lines = True
                else:
                    write_lines = False
            if write_lines:
                outf.write(line)

    assert len(transcripts.difference(transcripts_seen))==0, (
        "not all transcripts in gtf were found in the transcriptome fasta! %s, %s" %(
            len(transcripts), len(transcripts_seen)))

# ---------------------------------------------------

@transform(MakeMergedGTF,
           formatter(),
           "references/tx2gene.tsv")
def get_txp_to_gene(infile, outfile):
    ''' Make a map of trancript -> gene for use later'''

    import CGAT.GTF as GTF
    tx2gene = set()

    for entry in GTF.iterator(IOTools.openFile(infile)):
        tx2gene.add((entry.transcript_id, entry.gene_id))

    IOTools.writeLines(outfile, tx2gene)


# ---------------------------------------------------
@transform(MakeMergedGTF,
           formatter(),
           "references/mito_genes.tsv")
def get_mito_genes(infile, outfile):
    ''' Identify the Mitochondrial genes using the contig ID'''

    import CGAT.GTF as GTF
    mt_genes = set()

    for entry in GTF.iterator(IOTools.openFile(infile)):
        if entry.contig in ("hg_chrMT", "hg_chrM", "mm_chrM", "mm_chrM"):
            mt_genes.add(entry.gene_id)

    with IOTools.openFile(outfile, "w") as outf:
        outf.write("\n".join(mt_genes))


# ---------------------------------------------------
@transform(MakeMergedGTF,
           formatter(),
           "references/rrna_genes.tsv")
def get_rrna_genes(infile, outfile):
    ''' Identify the rRNA genes using gene biotype'''

    import CGAT.GTF as GTF
    rrna_genes = set()

    for entry in GTF.iterator(IOTools.openFile(infile)):
        if entry.gene_type == "rRNA":
            rrna_genes.add(entry.gene_id)

    with IOTools.openFile(outfile, "w") as outf:
        outf.write("\n".join(rrna_genes))


# ---------------------------------------------------
@mkdir("references")
@originate(["references/mm_genome.fasta",
            "references/hg_genome.fasta"])
def DownloadGenomes(outfile):
    ''' Download the individual genome fasta'''
    species = os.path.basename(outfile).split("_")[0]
    link = PARAMS["%s_genome_address" % species]

    statement = ''' wget %(link)s -O - | zcat > %(outfile)s '''

    P.run()


# ---------------------------------------------------
@mkdir("references")
@merge(DownloadGenomes, "references/hgmm_genome.fasta")
def MakeMergedGenomes(infiles, outfile):

    hg_infile = [x for x in infiles if "hg_" in os.path.basename(x)][0]
    mm_infile = [x for x in infiles if "mm_" in os.path.basename(x)][0]

    statement = '''
    sed 's/>chr/>hg_chr/' %(hg_infile)s
    > %(outfile)s;

    checkpoint ;

    sed 's/>chr/>mm_chr/' %(mm_infile)s
    >> %(outfile)s;
    '''

    P.run()


# ---------------------------------------------------
@jobs_limit(1)
@follows(MakeMergedGTF)
@transform((DownloadGenomes, MakeMergedGenomes),
           regex("references/(\S+)_genome.fasta"),
           add_inputs(r"references/\1_transcriptome.gtf"),
           r"references/\1/\1.star_index")
def MakeStarIndex(infiles, outfile):
    ''' Index the reference genome for STAR'''

    fasta, gtf = infiles
    index_dir = os.path.dirname(outfile)
    fasta_name = os.path.basename(fasta)

    if not os.path.isdir(os.path.dirname(outfile)):
        os.mkdir(os.path.dirname(outfile))

    job_threads = 30
    job_memory = "4G"

    statement = '''
    STAR --runMode genomeGenerate
    --genomeDir %(index_dir)s
    --genomeFastaFiles %(fasta)s
    --genomeSAindexNbases 14
    --genomeChrBinNbits 18
    --genomeSAsparseD 2
     --sjdbGTFfile %(gtf)s
    --sjdbOverhang %(star_tx_overhang)s
    --runThreadN %(job_threads)s
    --limitGenomeGenerateRAM 57993269973
    --outFileNamePrefix %(outfile)s.
    > %(outfile)s 2>&1 ;
    '''

    P.run()

@jobs_limit(1)
@follows(DownloadBasicGTF)
@transform(DownloadGenomes,
           regex("references/hg_genome.fasta"),
           add_inputs(DownloadBasicGTF),
           r"references/hgbasic/hgbasic.star_index")
def MakeBasicStarIndex(infiles, outfile):
    ''' Index the reference genome for STAR'''

    fasta, gtf = infiles
    index_dir = os.path.dirname(outfile)
    fasta_name = os.path.basename(fasta)

    if not os.path.isdir(os.path.dirname(outfile)):
        os.mkdir(os.path.dirname(outfile))

    job_threads = 30
    job_memory = "4G"

    statement = '''
    STAR --runMode genomeGenerate
    --genomeDir %(index_dir)s
    --genomeFastaFiles %(fasta)s
    --genomeSAindexNbases 14
    --genomeChrBinNbits 18
    --genomeSAsparseD 2
     --sjdbGTFfile %(gtf)s
    --sjdbOverhang %(star_tx_overhang)s
    --runThreadN %(job_threads)s
    --limitGenomeGenerateRAM 57993269973
    --outFileNamePrefix %(outfile)s.
    > %(outfile)s 2>&1 ;
    '''

    P.run()

# ---------------------------------------------------
@transform((DownloadTranscriptomes, MakeMergedTranscriptome),
           suffix(".fasta"),
           ".salmon.index")
def buildSalmonIndex(infile, outfile):
    '''Generate a salmon index of the merged transcritome
    for use with alevin. This is our first timed step'''

    job_options = PARAMS["timed_job_options"]
    job_threads = PARAMS["timed_job_threads"]
    job_memory = "2G"

    timing_file = P.snip(outfile,".salmon.index") + "_quant.time"

    statement = '''
    /usr/bin/time salmon index
    --gencode
    -p %(job_threads)s
    -t %(infile)s
    -i %(outfile)s
    -k %(salmon_kmer)s
    2> %(timing_file)s
    '''

    P.run()

@transform(MakeBasicTranscriptome,
           suffix(".fasta"),
           ".salmon.index")
def buildBasicSalmonIndex(infile, outfile):
    '''Generate a salmon index for the basic geneset'''

    job_options = PARAMS["timed_job_options"]
    job_threads = PARAMS["timed_job_threads"]
    job_memory = "2G"

    timing_file = P.snip(outfile,".salmon.index") + "_quant.time"

    statement = '''
    /usr/bin/time salmon index
    --gencode
    -p %(job_threads)s
    -t %(infile)s
    -i %(outfile)s
    -k %(salmon_kmer)s
    2> %(timing_file)s
    '''

    P.run()


# ---------------------------------------------------
@follows(download10xWhitelists, buildSalmonIndex)
@mkdir("processed", "processed/10X_fastqs/")
@transform(download10x,
           regex("raw/10X_fastqs/(\S+).fastq.1.gz"),
           add_inputs(r"raw/10X_fastqs/\1_whitelist.txt",
                      buildSalmonIndex, get_txp_to_gene),
           r"processed/10X_fastqs/\1_all_reads.fastq")
def ExtractAndCorrectReads(infiles, outfile):

    # TS: remove txp2gene when salmon updated to not require for fq dumping
    infile_1, whitelist, index, txp2gene = infiles
    infile_2 = infile_1.replace(".1.gz", ".2.gz")

    statement = '''
    %(alevin_bin)s alevin
    -lISR
    -1 %(infile_1)s
    -2 %(infile_2)s
    --chromium
    -o %(outfile)s_log.dir
    -i %(index)s
    --whitelist %(whitelist)s
    -p 10 --dumpfq --nosoftmap --noquant
    --tgMap %(txp2gene)s >
    %(outfile)s'''

    P.run()


@follows(download10x,
         MakeMergedGenomes,
         MakeMergedTranscriptome,
         MakeStarIndex,
         get_txp_to_gene,
         get_mito_genes,
         get_rrna_genes)
def PrepareInputs():
    pass


################################################################################
################################################################################
# Run Alevin
################################################################################
################################################################################

# TS: Definition of alevin option files is getting a bit
# long winded. Hence handling outside ruffus task
def makeAlevinFiles():

    alevin_options_outfile = [
        "options/alevin/%s-%s" % (x, PARAMS['timed_job_threads'])
        for x in P.asList(PARAMS['alevin_params_var'])] # add thread tasks

    alevin_options_outfile += [
        "options/alevin/%s-%s" % ("wo_whitelist", x)
        for x in P.asList(PARAMS['alevin_threads_var'])] # add alevin option tasks

    alevin_options_outfile = list(set(alevin_options_outfile)) # avoid duplicate jobs

    alevin_options_outfile = [
        x + "-%s.option" % rep for (x,rep) in itertools.product(
            list(set(alevin_options_outfile)),
            range(0, PARAMS['alevin_replicates']))] # replicate each job

    return(set(alevin_options_outfile))


@mkdir("options", "options/alevin")
@originate(makeAlevinFiles())
def make_options_files_alevin(outfile):
    ''' To run alevin with multiple sets of options, we are
    creating empty files to define the "types" of alevin run'''

    P.touch(outfile)


@mkdir("alevin")
@jobs_limit(3, "alevin")
@follows(download10xWhitelists)
@product(download10x,
         formatter(".+/(?P<RAW>\S+).fastq.1.gz"),
         make_options_files_alevin,
         formatter(".+/(?P<OPTIONS>\S+)-(?P<THREADS>\d+)-(?P<REPS>\d+).option"),
         add_inputs(buildSalmonIndex,
                    get_txp_to_gene,
                    get_mito_genes, get_rrna_genes,
                    "raw/10X_fastqs/{RAW[0][0]}_whitelist.txt"),
         "alevin/{RAW[0][0]}-{OPTIONS[1][0]}-{THREADS[1][0]}-{REPS[1][0]}/quant.time",
         "{OPTIONS[1][0]}", "{THREADS[1][0]}", "{RAW[0][0]}")
def run_alevin_options(infiles, timing_file, options, job_threads, sample_name):
    ''' run alevin with various options switched on/off and diff. # threads '''

    job_threads = int(job_threads)

    infiles, ix1, ix2, ix3, tx2gene, mito_genes, rrna_genes, whitelist = infiles
    bc_file, run_type_file = infiles

    reads_file = re.sub(".fastq.1.gz", ".fastq.2.gz", bc_file)

    outdir = os.path.dirname(timing_file)

    # very hacky way to identify the correct index!
    # note that buildSalmonIndex gets expanded into three index file names above
    species = TENX2INFO[sample_name]["species"]
    index = [x for x in (ix1, ix2, ix3) if os.path.basename(x).split("_")[0] == species][0]

    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    run_type2type_options = {
        "wo_whitelist": "",
        "w_whitelist": "--whitelist %s" % whitelist,
        "naive": "--naive --dumpcsvcounts",
        "usecor": "--usecorrelation",
        "dumpeq": "--dumpbarcodeeq",
        "dumpcsv":"--dumpcsvcounts",
        "dumpfeatures":"--dumpfeatures"}

    type_options = run_type2type_options[options]

    job_options = PARAMS["timed_job_options"]

    # at least 40Gb total, at least 1Gb per thread
    job_memory = str(max(40/job_threads, 1)) + "G"

    statement = '''
    /usr/bin/time -v %(alevin_bin)s alevin
    --chromium
    -lISR
    -1 %(bc_file)s
    -2 %(reads_file)s
    --out %(outdir)s
    --index %(index)s
    -p %(job_threads)s
    --tgMap %(tx2gene)s
    --mrna %(mito_genes)s
    --rrna %(rrna_genes)s
    %(type_options)s
    2> %(timing_file)s'''
    P.run()

    P.touch(timing_file)


@mkdir("alevin")
@jobs_limit(6, "alevin_subsets")
@follows(download10xWhitelists)
@transform(subsetRange,
           regex("raw/10X_fastqs/(\S+).fastq.1.gz"),
           add_inputs(buildSalmonIndex,
                      get_txp_to_gene,
                      get_mito_genes, get_rrna_genes,
                      "raw/10X_fastqs/%s_whitelist.txt" % PARAMS['subset_sample']),
           r"alevin/\1-w_whitelist-%s-0/quant.time" % (PARAMS['timed_job_threads']))
def run_alevin_subsets(infiles, timing_file):
    ''' run alevin on the subsampled data'''

    bc_file, ix1, ix2, ix3, tx2gene, mito_genes, rrna_genes, whitelist = infiles
    reads_file = re.sub(".fastq.1.gz", ".fastq.2.gz", bc_file)

    outdir = os.path.dirname(timing_file)

    # very hacky way to identify the correct index!
    # note that buildSalmonIndex gets expanded into three index file names above
    species = TENX2INFO[PARAMS['subset_sample']]["species"]
    index = [x for x in (ix1, ix2, ix3) if os.path.basename(x).split("_")[0] == species][0]

    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    job_options = PARAMS["timed_job_options"]
    job_threads = int(PARAMS['timed_job_threads'])

    # at least 40Gb total, at least 1Gb per thread
    job_memory = str(max(40/job_threads, 1)) + "G"

    #statement = '''
    #/usr/bin/time -v %(alevin_bin)s alevin

    statement = '''
    /usr/bin/time -v
    %(alevin_bin)s alevin
    --chromium
    -lISR
    -1 %(bc_file)s
    -2 %(reads_file)s
    --out %(outdir)s
    --index %(index)s
    -p %(job_threads)s
    --tgMap %(tx2gene)s
    --mrna %(mito_genes)s
    --rrna %(rrna_genes)s
    --whitelist %(whitelist)s
    2> %(timing_file)s'''
    P.run()

    #second run with dump options included
    outdir_2 = outdir.replace("w_whitelist", "dump_quants")
    if not os.path.isdir(outdir_2):
        os.mkdir(outdir_2)

    statement = '''
    /usr/bin/time -v
    %(alevin_bin)s alevin
    --chromium
    -lISR
    -1 %(bc_file)s
    -2 %(reads_file)s
    --out %(outdir_2)s
    --index %(index)s
    -p %(job_threads)s
    --tgMap %(tx2gene)s
    --mrna %(mito_genes)s
    --rrna %(rrna_genes)s
    --whitelist %(whitelist)s
    --dumpbarcodeeq
    --dumpcsvcounts
    2> %(timing_file)s_dump_quants'''
    P.run()

    #third run with noem option to derive "unique" quantification
    outdir_3 = outdir.replace("w_whitelist", "noem")
    if not os.path.isdir(outdir_3):
        os.mkdir(outdir_3)

    statement = '''
    /usr/bin/time -v
    %(alevin_bin)s alevin
    --chromium
    -lISR
    -1 %(bc_file)s
    -2 %(reads_file)s
    --out %(outdir_3)s
    --index %(index)s
    -p %(job_threads)s
    --tgMap %(tx2gene)s
    --mrna %(mito_genes)s
    --rrna %(rrna_genes)s
    --dumpbarcodeeq
    --dumpcsvcounts
    --noem
    2> %(timing_file)s_noem'''
    #P.run()



@mkdir("alevin")
@jobs_limit(3, "alevin")
@follows(download10xWhitelists)
@transform("raw/10X_fastqs/%s.fastq.1.gz" % PARAMS['subset_sample'],
           regex("raw/10X_fastqs/(\S+).fastq.1.gz"),
           add_inputs(buildBasicSalmonIndex,
                      get_txp_to_gene,
                      get_mito_genes, get_rrna_genes),
           r"alevin/\1_basic/quant.time")
def run_alevin_basic(infiles, timing_file):
    ''' run alevin with the base geneset'''

    bc_file, index, tx2gene, mito_genes, rrna_genes = infiles

    reads_file = re.sub(".fastq.1.gz", ".fastq.2.gz", bc_file)

    outdir = os.path.dirname(timing_file)

    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    
    job_threads = PARAMS["timed_job_threads"]
    job_options = PARAMS["timed_job_options"]

    # at least 40Gb total, at least 1Gb per thread
    job_memory = str(max(40/job_threads, 1)) + "G"

    statement = '''
    /usr/bin/time -v %(alevin_bin)s alevin
    --chromium
    -lISR
    -1 %(bc_file)s
    -2 %(reads_file)s
    --out %(outdir)s
    --index %(index)s
    -p %(job_threads)s
    --tgMap %(tx2gene)s
    --mrna %(mito_genes)s
    --rrna %(rrna_genes)s
    --dumpcsvcounts
    2> %(timing_file)s'''
    P.run()

    P.touch(timing_file)

################################################################################
################################################################################
# Run CellRanger
################################################################################
################################################################################

@merge((DownloadGenomes, DownloadGTFs), ["hg38/.sentinel",
                                         "mm10/.sentinel",
                                         "hg38_and_mm10/.sentinel"])
def make_cellranger_transcriptome_dir(infiles, sentinels):
    ''' cell ranger requires the input files in a single directory
    with specific structure '''

    gtfs = [x for x in infiles if x.endswith(".gtf")]
    fas = [x for x in infiles if x.endswith(".fasta")]

    hg_fa = [x for x in fas if "hg_" in x][0]
    mm_fa = [x for x in fas if "mm_" in x][0]
    hg_gtf = [x for x in gtfs if "hg_" in x][0]
    mm_gtf = [x for x in gtfs if "mm_" in x][0]

    job_memory = "4G"
    job_threads = 30

    statements = []

    statements.append('''
    cellranger mkref --nthreads=%(job_threads)s --memgb=40
    --genome=hg38 --fasta=%(hg_fa)s --genes=%(hg_gtf)s
    ''' % locals())

    statements.append('''
    cellranger mkref --nthreads=%(job_threads)s --memgb=40
    --genome=mm10 --fasta=%(mm_fa)s --genes=%(mm_gtf)s
    ''' % locals())

    statements.append('''
    cellranger mkref --nthreads=%(job_threads)s --memgb=40
    --genome=hg38 --fasta=%(hg_fa)s --genes=%(hg_gtf)s
    --genome=mm10 --fasta=%(mm_fa)s --genes=%(mm_gtf)s
    ''' % locals())

    P.run()

    for s in sentinels:
        P.touch(s)


@merge((DownloadGenomes, DownloadBasicGTF),
       "hg38basic/.sentinel")
def make_basic_cellranger_transcriptome_dir(infiles, sentinel):
    ''' cell ranger requires the input files in a single directory
    with specific structure '''

    hg_gtf = [x for x in infiles if x.endswith(".gtf")][0]
    hg_fa = [x for x in infiles if os.path.basename(x)=="hg_genome.fasta"][0]

    job_memory = "4G"
    job_threads = 30

    statements = []

    statements.append('''
    cellranger mkref --nthreads=%(job_threads)s --memgb=40
    --genome=hg38basic --fasta=%(hg_fa)s --genes=%(hg_gtf)s
    ''' % locals())

    P.run()
    P.touch(sentinel)


# TS: Definition of cellranger option files is getting a bit
# long winded. Hence handling outside ruffus task
def makeCellrangerFiles():

    # add thread tasks
    cellranger_options_outfile = [
        "options/cellranger/%s-%s" % (x, PARAMS['timed_job_threads'])
        for x in P.asList(PARAMS['cellranger_params_var'])]

    # add alevin option tasks
    cellranger_options_outfile += [
        "options/cellranger/%s-%s" % ("wo_whitelist", x)
        for x in P.asList(PARAMS['cellranger_threads_var'])]

    # avoid duplicate jobs
    cellranger_options_outfile = list(set(cellranger_options_outfile))

    # replicate each job
    cellranger_options_outfile = [
        x + "-%s.option" % rep for (x,rep) in itertools.product(
            list(set(cellranger_options_outfile)),
            range(0, PARAMS['cellranger_replicates']))]

    return(set(cellranger_options_outfile))


@mkdir("options", "options/cellranger")
@originate(makeCellrangerFiles())
def make_cellranger_option_files(outfile):
    ''' In order to run cell ranger with multiple sets of options, here we
    create empty files to define the "types" of cellranger run. These
    option 'files' are used with each input dataset in all
    combinations in run_cellranger task '''
    outdir = os.path.dirname(outfile)

    P.touch(outfile)


@jobs_limit(2, "cellranger")
@follows(download10xWhitelists)
@mkdir("cellranger")
@product(download10x,
         formatter(".+/(?P<RAW>\S+).fastq.1.gz"),
         make_cellranger_option_files,
         formatter(".+/(?P<OPTIONS>\S+)-(?P<THREADS>\d+)-(?P<REPS>\d+).option"),
         add_inputs(make_cellranger_transcriptome_dir,
                    r"raw/10X_fastqs/{RAW[0][0]}_whitelist.txt",),
         "cellranger/{RAW[0][0]}-{OPTIONS[1][0]}-{THREADS[1][0]}-{REPS[1][0]}/quant.time",
         "{OPTIONS[1][0]}", "{THREADS[1][0]}", "{RAW[0][0]}")
def run_cellranger_options(infiles, timing_file, options, job_threads, sample_name):

    job_threads = int(job_threads)

    infiles, transcriptome_sentinels, whitelist = infiles
    infile, run_type_file = infiles

    # ugly hacky way to identify the correct index!
    species = TENX2INFO[sample_name]["species"]
    if species in ["hg", "mm"]:
        transcriptome_sentinel = [
            x for x in transcriptome_sentinels if
            os.path.basename(os.path.dirname(x))[0:2] == species][0]
    else:
        transcriptome_sentinel = [
            x for x in transcriptome_sentinels if
            os.path.basename(os.path.dirname(x)) == "hg38_and_mm10"][0]

    transcriptome_dir = os.path.dirname(transcriptome_sentinel)

    tmpdir = os.path.basename(P.getTempDir("."))
    outdir = os.path.dirname(timing_file)

    # empty the directory if it exists
    if os.path.isdir(outdir):
        shutil.rmtree(outdir)

    os.mkdir(outdir)

    n_cells = 0
    for line in open(whitelist, "r"):
        n_cells += 1

    sample_name = P.snip(os.path.basename(infile), ".fastq.1.gz")
    fastq_dir = sample_name + ".dir"

    run_type2type_options = {
        "wo_whitelist": "--expect-cells %s" % n_cells,
        "w_whitelist": "--force-cells %s" % n_cells}

    type_options = run_type2type_options[options]

    job_options = PARAMS["timed_job_options"]
    job_memory = str(int(PARAMS["cellranger_memory"]) / job_threads) + "G"

    localmem = int(PARAMS["cellranger_memory"])

    statement = '''
    /usr/bin/time -v cellranger count
    --id %(tmpdir)s
    --fastqs %(fastq_dir)s
    --localcores %(job_threads)s
    --localmem=%(localmem)s
    --transcriptome=%(transcriptome_dir)s
    --nosecondary
    %(type_options)s
    > %(timing_file)s.log 2> %(timing_file)s ; checkpoint ;
    mv -f %(tmpdir)s/* %(outdir)s
    '''

    P.run()


@mkdir("cellranger")
@jobs_limit(2, "cellranger")
@follows(download10xWhitelists)
@transform(subsetRange,
           regex("raw/10X_fastqs/(\S+).fastq.1.gz"),
           add_inputs(make_cellranger_transcriptome_dir,
                      "raw/10X_fastqs/%s_whitelist.txt" % PARAMS['subset_sample'],
                      r"raw/10X_fastqs/\1.dir"),
           r"cellranger/\1-w_whitelist-%s-0/quant.time" % (PARAMS['timed_job_threads']))
def run_cellranger_subsets(infiles, timing_file):

    infile, transcriptome_sentinels, whitelist, fastq_dir = infiles

    # ugly hacky way to identify the correct index!
    species = TENX2INFO[PARAMS['subset_sample']]["species"]
    if species in ["hg", "mm"]:
        transcriptome_sentinel = [
            x for x in transcriptome_sentinels if
            os.path.basename(os.path.dirname(x))[0:2] == species][0]
    else:
        transcriptome_sentinel = [
            x for x in transcriptome_sentinels if
            os.path.basename(os.path.dirname(x)) == "hg38_and_mm10"][0]

    transcriptome_dir = os.path.dirname(transcriptome_sentinel)

    tmpdir = os.path.basename(P.getTempDir("."))
    outdir = os.path.dirname(timing_file)

    # empty the directory if it exists
    if os.path.isdir(outdir):
        shutil.rmtree(outdir)

    os.mkdir(outdir)

    n_cells = 0
    for line in open(whitelist, "r"):
        n_cells += 1

    job_threads = PARAMS['timed_job_threads']
    job_options = PARAMS["timed_job_options"]
    job_memory = str(int(PARAMS["cellranger_memory"]) / job_threads) + "G"

    localmem = int(PARAMS["cellranger_memory"])

    statement = '''
    /usr/bin/time -v cellranger count
    --id %(tmpdir)s
    --fastqs %(fastq_dir)s
    --localcores %(job_threads)s
    --localmem=%(localmem)s
    --transcriptome=%(transcriptome_dir)s
    --nosecondary
    --force-cells %(n_cells)s
    > %(timing_file)s.log 2> %(timing_file)s ; checkpoint ;
    mv -f %(tmpdir)s/* %(outdir)s
    '''

    P.run()

@mkdir("cellranger")
@jobs_limit(2, "cellranger")
@follows(download10xWhitelists)
@transform("raw/10X_fastqs/%s.fastq.1.gz" % PARAMS['subset_sample'],
           regex("raw/10X_fastqs/(\S+).fastq.1.gz"),
           add_inputs(make_basic_cellranger_transcriptome_dir,
                      "raw/10X_fastqs/%s_whitelist.txt" % PARAMS['subset_sample']),
           r"cellranger/\1_basic/quant.time")
def run_cellranger_basic(infiles, timing_file):

    infile, transcriptome_sentinel, whitelist = infiles

    transcriptome_dir = os.path.dirname(transcriptome_sentinel)

    tmpdir = os.path.basename(P.getTempDir("."))
    outdir = os.path.dirname(timing_file)

    # empty the directory if it exists
    if os.path.isdir(outdir):
        shutil.rmtree(outdir)

    os.mkdir(outdir)

    n_cells = 0
    for line in open(whitelist, "r"):
        n_cells += 1

    job_threads = PARAMS['timed_job_threads']
    job_options = PARAMS["timed_job_options"]
    job_memory = str(int(PARAMS["cellranger_memory"]) / job_threads) + "G"

    localmem = int(PARAMS["cellranger_memory"])

    sample_name = P.snip(os.path.basename(infile), ".fastq.1.gz")
    fastq_dir = sample_name + ".dir"

    statement = '''
    /usr/bin/time -v cellranger count
    --id %(tmpdir)s
    --fastqs %(fastq_dir)s
    --localcores %(job_threads)s
    --localmem=%(localmem)s
    --transcriptome=%(transcriptome_dir)s
    --nosecondary
    --force-cells %(n_cells)s
    > %(timing_file)s.log 2> %(timing_file)s ; checkpoint ;
    mv -f %(tmpdir)s/* %(outdir)s
    '''

    P.run()

################################################################################
# Run UMI-Tools
################################################################################

@jobs_limit(2, "umi")
@follows(mkdir("umi_tools"), download10xWhitelists)
@transform(download10x,
           formatter(".+/(?P<TRACK>.+).fastq.1.gz"),
           add_inputs(r"{path[0]}/{TRACK[0]}_whitelist.txt"),
           r"umi_tools/{TRACK[0]}.extracted.time")
def umi_tools_extract(infiles, timing_file):
    '''Run UMI-Tools extraction using 10X whitelist'''

    infastq1, whitelist = infiles

    infastq2 = re.sub(".fastq.1.gz", ".fastq.2.gz", infastq1)
    out1 = re.sub(".time", ".fastq.1", timing_file)
    out2 = re.sub(".time", ".fastq.2", timing_file)
    logfile = re.sub(".time", ".log", timing_file)

    job_options = PARAMS["timed_job_options"]

    statement = '''
    /usr/bin/time -v umi_tools extract
    -I %(infastq1)s
    --read2-stdout
    -S %(out2)s
    --read2-in %(infastq2)s
    --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN
    --filter-cell-barcode
    --whitelist=%(whitelist)s
    -L %(logfile)s
    2> %(timing_file)s '''

    P.run()


@jobs_limit(2, "umi")
@mkdir("umi_tools")
@transform(download10x,
           formatter(".+/(?P<TRACK>.+).fastq.1.gz"),
           r"umi_tools/{TRACK[0]}.whitelist.time")
def umi_tools_whitelist(infastq1, timing_file):
    '''Run UMI-Tools extraction using 10X whitelist'''

    infastq2 = re.sub(".fastq.1.gz", ".fastq.2.gz", infastq1)
    whitelist_out = timing_file.replace(".time", "")
    logfile = re.sub(".time", ".log", timing_file)

    job_options = PARAMS["timed_job_options"]

    statement = '''
    /usr/bin/time -v umi_tools whitelist
    -I %(infastq1)s
    -S %(whitelist_out)s
    --read2-in %(infastq2)s
    --plot-prefix=%(whitelist_out)s
    --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN
    --subset-reads=10000000
    --error-correct-threshold=0
    -L %(logfile)s
    2> %(timing_file)s '''

    P.run()


@jobs_limit(2, "umi")
@follows(mkdir("umi_tools"), download10xWhitelists)
@transform(subsetRange,
           formatter(".+/(?P<TRACK>.+).fastq.1.gz"),
           add_inputs(r"{path[0]}/%s_whitelist.txt" % PARAMS['subset_sample']),
           r"umi_tools/{TRACK[0]}.extracted.time")
def umi_tools_extract_subset(infiles, timing_file):
    '''Run UMI-Tools extraction using 10X whitelist'''

    infastq1, whitelist = infiles

    infastq2 = re.sub(".fastq.1.gz", ".fastq.2.gz", infastq1)
    out1 = re.sub(".time", ".fastq.1", timing_file)
    out2 = re.sub(".time", ".fastq.2", timing_file)
    logfile = re.sub(".time", ".log", timing_file)

    job_options = PARAMS["timed_job_options"]

    statement = '''
    /usr/bin/time -v umi_tools extract
    -I %(infastq1)s
    --read2-stdout
    -S %(out2)s
    --read2-in %(infastq2)s
    --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN
    --filter-cell-barcode
    --whitelist=%(whitelist)s
    -L %(logfile)s
    2> %(timing_file)s '''

    P.run()


################################################################################
@jobs_limit(2, "umi")
@transform((umi_tools_extract, umi_tools_extract_subset),
           regex("umi_tools/(.+).extracted.time"),
           inputs((r"umi_tools/\1.extracted.fastq.2",
                   MakeStarIndex)),
           r"umi_tools/\1.star.time", r"\1")
def map_with_star(infiles, outfile, sample_name):
    '''Map extracted umi_tools reads to merged genome'''

    fastq2, ix1, ix2, ix3 = infiles

    # very hacky way to identify the correct index!
    # note that buildSalmonIndex gets expanded into three index file names above
    if "subset" in sample_name:
        sample_name = PARAMS['subset_sample']
    species = TENX2INFO[sample_name]["species"]
    index = [x for x in (ix1, ix2, ix3) if
             os.path.basename(os.path.dirname(x)) == species][0]

    job_options = PARAMS["timed_job_options"]
    job_threads = int(PARAMS["timed_job_threads"])
    job_memory = str(max((50 / job_threads), 3))  + "G"
    outprefix = P.snip(outfile, ".star.time")
    index = os.path.dirname(index)

    statement = '''/usr/bin/time -v
                   STAR
                   --runThreadN %(job_threads)i
                   --genomeDir %(index)s
                   --readFilesIn %(fastq2)s
                   --outFilterMultimapNmax 1
                   --outSAMtype BAM SortedByCoordinate
                   --outFileNamePrefix %(outprefix)s.
                   2> %(outfile)s;'''

    P.run()
    #IOTools.zapFile(fastq2)



@jobs_limit(2, "umi")
@transform(umi_tools_extract,
           regex("umi_tools/(%s).extracted.time" % PARAMS['subset_sample']),
           inputs((r"umi_tools/\1.extracted.fastq.2",
                   MakeBasicStarIndex)),
           r"umi_tools/\1_basic.star.time")
def map_with_star_basic(infiles, outfile):
    '''Map extracted umi_tools reads to merged genome'''

    fastq2, index = infiles

    job_options = PARAMS["timed_job_options"]
    job_threads = int(PARAMS["timed_job_threads"])
    job_memory = str(max((50 / job_threads), 3))  + "G"
    outprefix = P.snip(outfile, ".star.time")
    index = os.path.dirname(index)

    statement = '''/usr/bin/time -v
                   STAR
                   --runThreadN %(job_threads)i
                   --genomeDir %(index)s
                   --readFilesIn %(fastq2)s
                   --outFilterMultimapNmax 1
                   --outSAMtype BAM SortedByCoordinate
                   --outFileNamePrefix %(outprefix)s.
                   2> %(outfile)s;'''

    P.run()
    #IOTools.zapFile(fastq2)
################################################################################
@follows(MakeMergedGTF)
@jobs_limit(2, "umi")
@transform((map_with_star, map_with_star_basic),
           regex("umi_tools/(.+).star.time"),
           inputs(r"umi_tools/\1.Aligned.sortedByCoord.out.bam"),
           r"umi_tools/\1.assign_genes.time", r"\1")
def assign_genes(inbam, outfile, sample_name):
    '''Use featureCounts to assign each read to a gene'''

    outbam = P.snip(outfile, ".time")
    job_options = PARAMS["timed_job_options"]
    job_threads = PARAMS["timed_job_threads"]
    job_memory = "5G"

    if "subset" in sample_name:
        sample_name = PARAMS['subset_sample']
    try:
        species = TENX2INFO[sample_name]['species']
        geneset = "references/%s_transcriptome.gtf" % species
    except:
        geneset = "references/hgbasic_transcriptome.gtf"

    statement = '''
                  checkpoint;
                  /usr/bin/time -v
                  featureCounts
                     -a %(geneset)s
                     -o %(outbam)s
                     -T %(job_threads)s
                     -s 1
                     -R BAM
                     %(inbam)s
                  2> %(outfile)s    ;
                  checkpoint;
                  mv %(inbam)s.featureCounts.bam %(outbam)s.bam'''

    P.run()
    #IOTools.zapFile(inbam)


################################################################################
@jobs_limit(2, "umi")
@transform(assign_genes,
           regex("umi_tools/(.+).assign_genes.time"),
           inputs(r"umi_tools/\1.assign_genes.bam"),
           r"umi_tools/\1.sorting.time")
def sort_and_index(infile, outfile):
    '''Sort and index the gene assigned bamfile'''

    outbam = P.snip(outfile, ".time") + ".bam"
    job_options = PARAMS["timed_job_options"]
    job_threads = PARAMS["timed_job_threads"]

    statement = '''/usr/bin/time -v
                   bash -c '
                   samtools sort %(infile)s -o %(outbam)s -@ %(job_threads)s;
                   samtools index %(outbam)s'
                   2> %(outfile)s '''

    P.run()
    #IOTools.zapFile(infile)

################################################################################
@jobs_limit(2, "umi")
@transform(sort_and_index,
           regex("umi_tools/(.+).sorting.time"),
           inputs(r"umi_tools/\1.sorting.bam"),
           r"umi_tools/\1.counting.time")
def run_umi_tools_count(infile, outfile):
    '''Time the umi_tools count command'''

    out_matrix = P.snip(outfile, "ing.time") + ".tsv"

    job_options = PARAMS["timed_job_options"]

    statement = '''/usr/bin/time -v
                      umi_tools count
                      -I %(infile)s
                      -S %(out_matrix)s
                      --per-cell
                      --per-gene
                      --gene-tag=XT
                      --wide-format-cell-counts
                      --log=%(outfile)s.log
                   2> %(outfile)s'''
    P.run()
    #IOTools.zapFile(infile)

################################################################################
# Run kallisto
################################################################################
@follows(mkdir("kallisto"))
@transform(MakeMergedTranscriptome,
           formatter(),
           r"kallisto/{basename[0]}.index")
def build_kallisto_index(infile, outfile):

    job_memory = "10G"

    statement = '''kallisto index -i %(outfile)s %(infile)s'''

    P.run()


################################################################################
@jobs_limit(2, "kallisto")
@follows(mkdir("kallisto"))
@transform(download10xV1,
           regex(".+/(.+).fastq.1.gz"),
           add_inputs(build_kallisto_index),
           r"kallisto/\1/config.json")
def generate_kallisto_config(infiles, outfile):

    fastq, kallisto_index = infiles

    kallisto_index = os.path.abspath(kallisto_index)
    kallisto_source = os.path.abspath(
        os.path.join(os.path.dirname(__file__),
                     "scRNA-Seq-TCC-prep/source/"))

    base_dir = os.path.abspath(P.snip(fastq, ".fastq.1.gz"))+".dir/"
    track = P.snip(outfile, "config.json")
    if not os.path.isdir(track):
        os.mkdir(track)
    save_dir = os.path.abspath(os.path.join(track, "save_dir/"))
    output_dir = os.path.abspath(os.path.join(track, "output_dir/"))
    TCC_output = os.path.abspath(os.path.join(track, "tcc_dir/"))

    import glob
    fastqs=glob.glob(base_dir + "*.fastq.gz")
    barcodes = set([re.search("-([TCGA]{8})_",x).groups()[0] for x in fastqs])
    barcodes = '["' + '","'.join(barcodes) + '"]'
    kallisto_binary = "kallisto"
    threads = PARAMS["timed_job_threads"]

    template = os.path.join(os.path.dirname(__file__),
                            "template_config.json")

    template = open(template).read()
    template = template % locals()

    open(outfile, "w").write(template)


################################################################################
@jobs_limit(2, "kallisto")
@transform(generate_kallisto_config,
           suffix("/config.json"),
           ".cell_barcodes.time")
def get_cell_barcodes(infile, outfile):

    job_options = PARAMS["timed_job_options"]
    job_threads = PARAMS["timed_job_threads"]
    job_memory = "10G" # just to be safe, let's provide plenty of mem.

    kallisto_source = os.path.abspath(
        os.path.join(os.path.dirname(__file__),
                     "scRNA-Seq-TCC-prep/source/"))

    track = P.snip(outfile, ".cell_barcodes.time")
    outfile = os.path.abspath(outfile)
    infile = os.path.abspath(infile)
    statement = '''
                    cd %(track)s;
                   /usr/bin/time -v
                   python %(kallisto_source)s/get_cell_barcodes.py
                   %(infile)s
                   2> %(outfile)s''' %locals()

    P.run()


################################################################################
@jobs_limit(2, "kallisto")
@transform(get_cell_barcodes,
           suffix(".cell_barcodes.time"),
           add_inputs(generate_kallisto_config),
           r".correct_and_split.time")
def error_correct_and_split(infiles, outfile):

    job_options = PARAMS["timed_job_options"]
    job_threads = PARAMS["timed_job_threads"]
    job_memory = str(200/job_threads) + "G" # just to be safe, let's provide plenty of mem.

    kallisto_source = os.path.abspath(
        os.path.join(os.path.dirname(__file__),
                     "scRNA-Seq-TCC-prep/source/"))

    _, config_file = infiles
    outfile = os.path.abspath(outfile)
    basedir = P.snip(outfile, ".correct_and_split.time")
    config_file = os.path.abspath(config_file)

    statement = '''cd %(basedir)s;
                   /usr/bin/time -v
                   python %(kallisto_source)s/error_correct_and_split.py
                          %(config_file)s
                   2> %(outfile)s'''

    P.run()


################################################################################
@jobs_limit(2, "kallisto")
@transform(error_correct_and_split,
           suffix(".correct_and_split.time"),
           add_inputs(generate_kallisto_config),
           r".compute_TCCs.time")
def compute_TCCs(infiles, outfile):

    _, config_file = infiles
    config_file = os.path.abspath(infile)
    job_options = PARAMS["timed_job_options"]
    job_threads = PARAMS["timed_job_threads"]
    job_memory = "50G" # just to be safe, let's provide plenty of mem.

    kallisto_source = os.path.abspath(
        os.path.join(os.path.dirname(__file__),
                     "scRNA-Seq-TCC-prep/source/"))

    outfile = os.path.abspath(outfile)
    basedir = P.snip(outfile, ".compute_TCCs.time")
    statement = '''cd %(basedir);
                   /usr/bin/time -v
                   python %(kallisto_source)s/compute_TCCs.py
                          %(config_file)s
                   2> %(outfile)s'''

    P.run()


################################################################################
# Full tasks
################################################################################
def extractFromTime(infiles, outfile, multistage=False):
    ''' extract run time / memory / cpu information from time output'''
    outf = IOTools.openFile(outfile, "w")
    outf.write("%s\n" % "\t".join(("sample", "tool", "stage", "options", "threads",
                                   "rep", "time_in_seconds", "max_mem_gb", "cpu")))

    for infile in infiles:

        max_mem, clocktime, cpu = None, None, None

        if multistage: # no replicates for kallisto/naive

            sample, stage, _ = os.path.basename(infile).split(".")
            threads = PARAMS['timed_job_threads']
            rep = 0
            tool = infile.split("/")[-2]
            options = ""


        else:
            sample, options, threads, rep = os.path.basename(
                os.path.dirname(infile)).split("-")
            stage = "complete"
            tool = infile.split("/")[-3]

        with IOTools.openFile(infile, "r") as inf:
            for line in inf:
                if "Maximum resident set size (kbytes):" in line:
                    max_mem = line.strip().split("(kbytes): ")[1]
                elif "Elapsed (wall clock) time (h:mm:ss or m:ss):" in line:
                    clocktime = line.strip().split("m:ss): ")[1]
                elif "Percent of CPU this job got:" in line:
                    cpu = line.strip().split("this job got: ")[1]
                else:
                    pass

        if max_mem and clocktime and cpu:
            h, m, s = 0.0, 0.0, 0.0
            clocktime_s = clocktime.split(':')
            if len(clocktime_s) == 2:
                m += float(clocktime_s[0])
                s += float(clocktime_s[1])
            elif len(clocktime_s) == 3:
                h += float(clocktime_s[0])
                m += float(clocktime_s[1])
                s += float(clocktime_s[2])
            else:
                raise ValueError("can't parse the elapsed time: %s" % clocktime)

            time_in_seconds = s
            time_in_seconds += (m * 60)
            time_in_seconds += (h * 3600) # convert to seconds

            time_in_seconds = round(time_in_seconds, 0)

            #For cellranger, need to use the _perf files to obtain the mem usage
            if "cellranger" in infile:
                maxrss = 0

                data = json.load(open(os.path.join(os.path.dirname(infile), "_perf")))

                for c_stage in data:
                    for x in c_stage['forks']:
                        maxrss = max(maxrss, (x['fork_stats']['maxrss']))

                #index_maxrss = data[0]['forks'][0]['fork_stats']['maxrss'] * 1024

                #max_bytes = 0
                #for b_hist in data[0]['bytehist']:
                #    max_bytes = max(max_bytes, b_hist['bytes'])

                max_mem = maxrss
                print(maxrss)

            max_mem =  round((int(max_mem) / (1024 * 1024)),4) # convert to Gb
            cpu = cpu.strip("%") # get rid of trailing "%"

            outf.write("%s\n" % "\t".join(map(str,
                (sample, tool, stage, options, threads,
                 rep, time_in_seconds, max_mem, cpu))))
        else:
            pass # TS: temp. allowing empty files
            #raise ValueError(
            #    "could not parse memory and clock_time for %s" % infile)



# currently excluding  run_alevin_V1 since we are not running V1 samples here
@mkdir("profile")
@merge((run_alevin_options, run_alevin_subsets,
        run_cellranger_options,run_cellranger_subsets),
       r"profile/profile.tsv")
def getProfile(infiles, outfile):
    ''' extract run info for single command methods = alevin and cellranger'''
    extractFromTime(infiles, outfile, multistage=False)


#currently excluding kallisto tasks since error_correct_and_split is returning errors
#@merge((umi_tools_extract,
#        map_with_star,
#        assign_genes,
#        sort_and_index,
#        run_umi_tools_count,
#        error_correct_and_split,
#        compute_TCCs,
#        get_cell_barcodes),
@mkdir("profile")
@merge((umi_tools_extract,
        map_with_star,
        assign_genes,
        sort_and_index,
        run_umi_tools_count),
       r"profile/multi_stage_profile.tsv")
def getMultiStageProfile(infiles, outfile):
    ''' extract run info for multi-stage methods = naive and kallisto'''
    extractFromTime(infiles, outfile, multistage=True)


@follows(getProfile, getMultiStageProfile)
def full():
    pass


@follows(mkdir("report"))
def build_report():
    '''build report from scratch.

    Any existing report will be overwritten.
    '''

    E.info("starting report build process from scratch")
    P.run_report(clean=True)

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
