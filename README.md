# PABLOG: a Primer Analysis tool using a Bee-Like approach on Orthologous Genes

PABLOG is a script tool developed in Python that implements an analysis pipeline for primer design.


## Installation
At the moment PABLOG is released as python script pablog.py to be run from command-line using the interpreter (in a future update may be release as conda package).

The simplest way to get the tool is by downloading or cloning the github repo and you are ready to go.

PABLOG requires the following python packages to work:
  - HTSeq
  - numpy
  - pandas
  - primer3-py
  - pysam
  - tqdm

A requirements.txt file is included in the folder and can be used with pip (pip install -r requirements.txt)

A good way to have all set up is by creating a new conda environment (Python >=3.9) and then use pip:

      conda create -n env_name python=3.9
      pip install -r requirements.txt

NOTE: PABLOG requires Python >=3.9 
 

## Usage

Run pablog.py script using python:

       python pablog.py FILE.bam Ref_Gene.fa output_result_file.txt [optional]
       
   where:

   - FILE.bam is the alignment file,
   - Ref_Gene.fa is the reference genome used in the alignment,
   - output_result_file.txt text file name where results are written in,
   - [optional] args are:
     - size [INT] (=60 by default) which is the minimun length we consider to filter out CIGAR operations in read analysis


## Guide

PABLOG analysis requires an alignment file and its reference sequence (.bam and .fa).
The pipeline consists of 4 main steps:

  1. Alignment analysis
  2. Coverage analysis
  3. Consensus sequence
  4. Primer3 Design

At each step a status message is printed, so the user can monitor the execution.

In the alignment analysis .bam file is scanned read-by-read and exon-exon junctions are identified (HTSeq python package is used to parse and read alignment).
The optional argument <em> size </em> suggests PABLOG to include only CIGAR operations of at least <em> size </em> bp.
Such value is set to 60 by default, as other tools do (i.e. RegTools), but user can eventually change as experiment requires. Important to underline is the fact that lower values for <em> size </em> are not recommended, since they would include also regions that are not actually regions of interest (false positives). Moreover, at the time of the release there is a little imperfection on how these regions are used to edit the consensus sequence (more details later). So in case the experiment would require a less restrictive analysis (more false positive to be included) we suggest to tune the alignment parameters accordingly rather than <em> size </em> value. 

Once the regions are found, the coverage analysis can compute a goodness score for each of them, giving a criteria for best-to and worst-to pick as primer design.

The score comes from coverage variation analysis in the region: we expect to have high coverage at the ends of the region and very low coverage inside the region, assessing a very high score. On the other hand, less variation in the coverage suggests the region is not very clear (probably for low quality sequencing data or suboptimal alignment parameters), hence low score is assigned. 

The next step in pipeline is the consensus sequence generation for Primer3 to design the primers. 
Since we are interested in getting primers for non-model species and we don't have a reference genome,
we use the alignment file to generate a sequence that contains the most frequent nucleotide base in every position from aligned reads (consensus sequence). Such operation is done with SAMtools mpileup function and the result is parsed manually to get nucleotides frequencies and so the most frequent one. 
PABLOG saves the consensus sequence in a different text file (*consensus_sequence.fa) as static reference to be used in other analysis. 

We want the primers to amplify exons areas (expression data) so we need to remove from consensus sequence the substrings corresponding to introns, which are exactly the regions PABLOG found in the steps above.

After that the primer design can be submitted to primer3. To do so we rely on primer3-py python package which wraps the tools to be used inside python.
In this phase we use a static file for primer3 parameters, that contains the default settings as intended by original authors (for more informations we suggest consult primer3 manual).
If the user requires a different configuration, the file can be edited as needed.

## PABLOG output

PABLOG's analysis results consist of a text file (the user can specify a particular filename in command line when running the script) divided in 3 sections.
The first section contains the regions found and goodness scores, sorted by genomic coordinates.
In the second section, the primer designed in each exon-exon junction is printed using primer3 output format (check manual), sorted by goodness score.
If was not possible to design a primer in a region, a status message will be printed instead.  

The last section contains the input template sequence used by primer3 (the consensus sequence stripped of introns) as reference for user. 

Examples of output can be found inside examples/ folder.  

## Examples
Some examples are available inside examples/ folder. We added experiments for DAG and TPL genes benchmark, each of them has all required files and also the PABLOG's results. The user can see those as reference or to test the tool functionality. 

In case of doubts or report, don't hesitate to submit an issue in the dedicated section.
Opinions and suggestions are welcome as well!

## Known issues / work in progress

This is the first release of PABLOG and some aspects are still work in progress. 

Regarding some of these, we are working on solving the incosistency of low <em> size </em> values, as one of the effects is producing too many regions that makes consensus-sequence's edit a little tricky (how manage semi/full overlapping regions? maybe fuse into one region? keep one and remove others? which one?).

An other aspect that is already work in progress, is about managing different primer3 settings for different regions. Right now every primer design share the same set of parameters values: let's say we found a region that did not produce a primer and probably a lower melting temperature (PRIMER_MIN_TM) would solve the problem, we would need to change the setting for all other primers as well.   

Moreover, a presets system maybe be implemented with specific configuration (alignment + primer3) for different use-cases (about this, feedback from pratical usage is crucial).

From tests and experiments performed, a good quality alignment with default settings produces very good results.
So we recommend to keep default settings and eventually tweak parameters on the alignemnt side.
