# PABLOG: a Primer Analysis tool using a Bee-Like approach on Orthologous Genes

This guide describes in details how to use the PABLOG tool to design primers useful for qRT-PCR in non-model species.

PABLOG is a script tool developed in Python that implements an analysis pipeline for primer design.


## Installation
At the moment PABLOG is released as a Python script, pablog.py, to be run from the command-line using the interpreter (in a future update, it may be released as conda package).

The simplest way to get the tool is by downloading or cloning the GitHub repo and you are ready to go.

PABLOG requires the following Python packages to work:
  - HTSeq
  - numpy
  - pandas
  - primer3-py
  - pysam
  - tqdm

A requirements.txt file is included in the folder and can be used with pip (pip install -r requirements.txt)

A good way to have it all set up is by creating a new conda environment (Python >=3.9) and then using pip:

      conda create -n env_name python=3.9
      pip install -r requirements.txt

NOTE: PABLOG requires Python >=3.9 
 

## Usage

Run the pablog.py script using Python:

       python pablog.py FILE.bam Ref_Gene.fa output_result_file.txt [optional]
       
   where:

   - FILE.bam is the alignment file,
   - Ref_Gene.fa is the reference genome used in the alignment,
   - output_result_file.txt is the text file name where results are written in,
   - [optional] args are:
     - size [INT] (=60 by default) which is the minimun length we consider to filter out CIGAR operations in read analysis


## Guide

PABLOG analysis requires an alignment file and its reference sequence (.bam and .fa).
The pipeline consists of four main steps:

  1. Alignment analysis
  2. Coverage analysis
  3. Consensus sequence
  4. Primer3 Design

At each step a status message is printed, so the user can monitor the execution.

In the alignment analysis .bam file is scanned read-by-read and exon-exon junctions are identified (the HTSeq Python package is used to parse and read alignment).
The optional argument <em> size </em> suggests PABLOG to include only CIGAR operations of at least <em> size </em> bp.
Such a value is set to 60 by default, as other tools do (i.e. RegTools), but the user can change it as the experiment requires. Important to underline is the fact that lower values for <em> size </em> are not recommended, since they would also include regions that are not actually regions of interest (false positives). 

Once the regions are found, the coverage analysis can compute a goodness score for each of them, giving criteria for best-to/worst-to pick as primer design.

The score comes from coverage variation analysis in the region: we expect to have high coverage at the ends of the region and very low coverage inside the region, assessing a very high score. On the other hand, less variation in the coverage suggests the region is not very clear (probably for low quality sequencing data or suboptimal alignment parameters), hence low score is assigned. 

The next step in the pipeline is the consensus sequence generation for Primer3 to design the primers. 
Since we are interested in getting primers for non-model species and don't have a reference genome,
we use the alignment file to generate a sequence that contains the most frequent nucleotide base in every position from aligned reads (consensus sequence). Such operation is done with the SAMtools mpileup function and the result is parsed manually to get nucleotides. 
PABLOG saves the consensus sequence in a different text file (*consensus_sequence.fa) as a static reference to be used in other analyses. 

We want the primers to amplify exons areas (expression data) so we need to remove from the consensus sequence the substrings corresponding to introns, which are exactly the regions PABLOG found in the steps above.

After that the primer design can be submitted to Primer3. To do so we rely on the primer3-py Python package which wraps the tools to be used via Python.
In this phase we use a static file for Primer3 parameters, that contains the default settings as intended by the original authors (for more informations we suggest consulting Primer3 manual).
If the user requires a different configuration, the file can be edited as needed.

## PABLOG output

PABLOG's analysis results consist of a text file (the user can specify a particular filename in the command line when running the script) divided into 3 sections.
The first section contains the regions found and goodness scores, sorted by genomic coordinates.
In the second section, the primer designed for each exon-exon junction is printed using Primer3 output format (check manual), sorted by goodness score.
If it was not possible to design a primer in a region, a status message will be printed instead.  

The last section contains the input template sequence used by primer3 (the consensus sequence stripped of introns) as a reference for the user. 

Examples of output can be found inside the examples/ folder.  

## Running examples

Some examples are available inside examples/ folder. We added experiments for benchmarking the DAG1 and TPL genes, and each of them has all required files as well as the PABLOG's outputs. The user can view those as references or to test the tool's functionality. In the following, we provide a step-by-step guide for PABLOG tool usage:

   1) Git clone the repository or download from https://github.com/tools4plant-omics/PABLOG
   2) Navigate to PABLOG/ folder and open a command line;
   3) Run pablog on DAG with:

            python pablog.py examples/DAG/SRR22407318_2_to_DAG_Daucus_CarotaAligned.sortedByCoord.out.bam examples/DAG/DAG_genome_ref_Daucuscarota.fa outputDAG.txt
   4) Run pablog on TPL with:
            
            python pablog.py examples/TPL/SRR22407318_2_to_TPL_SesamoAligned.sortedByCoord.out.bam examples/TPL/TPL_Sesamo_ref_gene.fa outputTPL.txt 50
      Note: in TPL example a smaller <em> size </em> value is required (50) to actually detect candidate regions  

In case of doubts or report, don't hesitate to submit an issue in the dedicated section. Opinions and suggestions are welcome as well!

## Known issues / work in progress

This is the first release of PABLOG and some aspects are still work in progress. 

We are working on solving the issues related to low <em> size </em> value, as one of the effects is producing too many regions, which makes consensus-sequence's edit a little tricky (how to manage semi or full overlapping regions? maybe fuse into one region? keep one and remove others? which one?).

An other aspect that is already in progress, is managing different Primer3 settings for different regions. Right now every primer design shares the same set of parameters values. Let's say we found a region that did not produce a primer and probably a lower melting temperature (PRIMER_MIN_TM) would solve the problem, we would need to change the setting for all other primers as well.   

Moreover, a presets system may be implemented with a specific configuration (alignment + primer3) for different use-cases (about this, feedback from practical usage is crucial).

From tests and experiments performed, a good quality alignment with default settings produces very good results.
So we recommend keeping the default settings and eventually tweaking parameters on the alignment side.

## Citation
The original work is now published and available at [doi.org/10.1111/ppl.14398](https://doi.org/10.1111/ppl.14398)
