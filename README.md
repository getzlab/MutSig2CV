# MutSig2CV

MutSig2CV ([Lawrence *et al.*, 2014](https://www.nature.com/articles/nature12912)).

## Overview

MutSig2CV analyzes somatic point mutations discovered in DNA sequencing, 
identifying genes mutated more often than expected by chance given inferred
background mutation processes.  MutSig2CV consists of three independent
statistical tests, described briefly below:

* **Abundance (CV)**: The most important step for inferring genes' mutational
significance is to properly classify whether the gene is highly mutated
relative to some background mutation rate (BMR), which varies on a macroscopic
level across patients and genes and on a microscopic level across sequence
contexts.  MutSig accounts for all three of these to renormalize BMR on a
per-gene, -patient, and -context level.

* **Clustering (CL)**: Genes often harbor mutational hotspots, specific sites that
are frequently mutated.  While abundance calculations bin mutations on the gene
level, clustering bins mutations on the local site level, which allows MutSig
to differentiate between genes with uniformly distributed mutations and genes
with localized hotspots, assigning higher significance to the latter.

* **Conservation (FN)**: MutSig uses evolutionary conservation as a proxy for
determining the functional significance of a mutated site.  It assumes that
genetic sites highly conserved across vertebrates have greater functional
significance than weakly conserved sites.  MutSig assigns a higher significance
to genes that experience frequent mutations in highly conserved sites.

For detailed descriptions of the algorithms employed in the MutSig2CV suite for
each of these tests, please visit
  https://www.broadinstitute.org/cancer/cga/mutsig


## Installing

MutSig is implemented in MATLAB. If you have a MATLAB installation and wish to run MutSig interactively on the MATLAB console, skip to the [Running](#running) section below. If you do not have MATLAB installed, or do not wish to run interactively, MutSig can be run as a standalone executable. The standalone executable is available for 64 bit Linux systems only, and requires that the MATLAB R2013a runtime be installed.
You can download and install the runtime environment from [here](https://ssd.mathworks.com/supportfiles/downloads/R2013a/deployment_files/R2013a/installers/glnxa64/MCR_R2013a_glnxa64_installer.zip). Runtime installation instructions can be found [here](http://www.mathworks.com/help/compiler/install-the-matlab-runtime.html).

Once the runtime is successfully installed, you must add it to your `LD_LIBRARY_PATH`.

```bash
MCRROOT=<path to runtime you specified when installing>
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MCRROOT/bin/glnxa64/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MCRROOT/sys/java/jre/glnxa64/jre/lib/amd64
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MCRROOT/sys/java/jre/glnxa64/jre/lib/amd64/server
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MCRROOT/sys/java/jre/glnxa64/jre/lib/amd64/native_threads
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MCRROOT/sys/os/glnxa64
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MCRROOT/bin/glnxa64
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MCRROOT/runtime/glnxa64
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MCRROOT/lib
```

MutSig requires ~3 GB of reference files. Since these files are too large to include in a GitHub repository, they are hosted elsewhere. Please download them from [here](http://software.broadinstitute.org/cancer/cga/sites/default/files/data/tools/mutsig/MutSig2CV.tar.gz), and copy the `reference` directory into this folder.

## Running <a name="running"></a>

To run on the MATLAB console, start MATLAB in this directory, and run:

```matlab
MutSig2CV(<path to mutations>, <path to output directory>, [params file])
```

To run the standalone application, `cd` to this directory, and run:

```bash
bin/MutSig2CV <path to mutations> <path to output directory> [params file]
```

MutSig looks for its reference files relative to this directory, so it is essential it is run here.

Each input is decribed below.

### Description of inputs

* **Mutations file**: Absolute path to the set of mutations to analyze.  Format is specified in the
following section, [Mutation Input Format](#mutation_inputs)

* **Output directory**: Absolute path to the directory where MutSig will save its output.  Will be
created if necessary.  NB: Any previous MutSig results in this directory will
be overwritten!  A description of the output files and formats are in the following section, [Outputs](#output_format).

* **Params file**: MutSig can take an optional parameters file, which allows for configuration of
  algorithm/run parameters.  User-configurable options are in the following section, [MutSig Configuration](#config).

### Mutation Input Format <a name="mutation_inputs"></a>

As input, MutSig takes a tab-delimited file with each line annotating a single
mutation in a single patient.  Columns can be in any order, with names and
formats as follows.  To provide maximal input flexibility, MutSig accepts
synonyms for each column name.  Column names are case sensitive.

* `chr`: Chromosome of the mutation.  MutSig only analyzes mutations on
         autosomal or sex chromosomes, and does not consider the
         mitochondrial chromosome or unplaced/alternate contigs.
  * Range: `(chr)?[1..24XY]`
  * Synonyms: Chromosome

* `pos` hg19 position of the mutation, 1-indexed.
  * Regex: `[0-9]+`
  * Synonyms: `Position`, `start`, `Start_position`

* `gene`: HUGO symbol for the gene containing this mutation, or "Unknown"
        for IGR mutations.
  * Regex: `[A-Za-z0-9]+`
  * Synonyms: `Hugo_Symbol`, `Gene_name`

* `patient`: Unique identifier for the patient.
  * Regex: `[A-Za-z0-9]+`
  * Synonyms: `Tumor_Sample_Barcode`, `Patient_name`

* `ref_allele`: hg19 reference base(s) for the position.  In the case of
                insertions, must be "-".
  * Regex: `(-|[ACGT]+)`
  * Synonyms: `Reference_Allele`

* `newbase`: Observed variant allele at the position.  In the case of
                deletions, must be "-".
  * Regex: `(-|[ACGT]+)`
  * Synonyms: `Tumor_Allele`, `Tum_allele`, `Alt_allele`, `Alternate_allele`, `Tumor_Seq_Allele2`

* `type`: Effect of mutation (e.g. missense, silent, nonsense, splice,
        UTR).  Full list of effects can be obtained as follows: 
        `awk 'NR > 1 { print $2 }' reference/mutation_type_dictionary.v4.txt | sort -u`
  * Regex: N/A
  * Synonyms: `Variant_Classification`

* classification: Mutation substitution type (e.g. SNP, INS, DEL.)
  * Regex: `(Complex_substitution|INS|DEL|[SDTO]NP)`
  * Synonyms: `Variant_Type`

If your input mutation data are missing gene/ref_allele/type/classification
fields, we recommend annotating using Oncotator, which will produce a
MutSig-ready MAF from a variety of input data.  Oncotator is available here:
  https://www.broadinstitute.org/cancer/cga/oncotator

## Outputs <a name="output_format"></a>

A MutSig run outputs several files:

### `sig_genes.txt`

A tab-delimited file containing all genes considered for analysis, sorted by
p-/q-values.  Columns are as follows:  

* `rank`: Position of the gene as sorted ascending by p-/q-value.  

* `gene`: HUGO symbol of the gene (RefSeq hg19)

* `longname`: HUGO description of the gene

* `codelen`: ORF length of the gene

* `nnei`: Number of neighboring genes in the bagel used to estimate BMR

* `nncd`: Number of noncoding mutations 

* `nsil`: Number of silent (synonymous) mutations in the gene

* `nmis`: Number of missense mutations in the gene

* `nstp`: Number of nonsense mutations in the gene

* `nspl`: Number of splice site mutations in the gene (defined as +/- 2 bases from
  the donor/acceptor site)

* `nind`: Number of insertions or deletions in the gene

* `nnon`: Number of nonsilent mutations in the gene (including all indels and splice
  site mutations, even if the codon change is synonymous in the latter case)

* `npat`: Number of patients with mutations in the gene

* `nsite`: Number of uniquely mutated sites in the gene (does not multiply count
  recurrently mutated positions)

* `pCV`: Abundance p-value

* `pCL`: Clustering p-value

* `pFN`: Functional (conservation) p-value

* `p`: Overall p-value obtained from Fisher combination of pCV, pCL, and pFN

* `q`: FDR-corrected (Benjamini-Hochberg) overall p-value

### `final_analysis_set.maf`

MutSig does not necessarily consider all mutations for analysis; for instance,
mutations in genes determined to be poorly covered, mutations determined to
belong to duplicate patients, or deep IGR mutations will be discarded.
final_analysis_set.maf contains only those mutations actually used for 
significance analysis.  Note that this file preserves the original input MAF's
columns but also contains columns used internally by MutSig.

### `patient_counts_and_rates.txt`

A tab-delimited file summarizing each patient's mutation counts, one patient
per line.  For each patient, fields are as follows:

* `name`: Patient name

* `nmut`: Total number of mutations in original input

* `n_coding`: Total number of coding mutations

* `n_coding_nonsilent`: Number of nonsilent coding mutations

* `callscheme_name`: Coverage model (e.g. whole genome, exome capture, exomeplus,
  TCGA coding only), inferred from distribution of patient's mutations

* `cov_idx`: Numerical index for this coverage model

* `N_tot`: Total number of covered bases for this coverage model

* `nsil_tot`: Number of silent mutations

* `nnon_tot`: Number of nonsilent mutations

* `n_tot`: Number of mutations considered for analysis

* `rate_sil`: Silent mutation rate

* `rate_non`: Nonsilent mutation rate

* `rate_tot`: Overall mutation rate

* `log_rate_tot`: -log10 of overall mutation rate

* `n_ind`: Number of indels

* `N_ind`: Total number of covered coding bases for this coverage model

* `rate_c`: Overall coding mutation rate

* `rate_ind`: Indel mutation rate 

### `mutcategs.txt`, `mutcateg_discovery.txt`

MutSig accounts for mutation rate heterogeneity across trinucleotide contexts
when calculating BMRs.  To increase statistical power, it clusters the 96 base
substitutions+trinucleotide contexts into k mutually exclusive categories via
entropy minimization with default k = 5.  mutcategs.txt and
mutcateg_discovery.txt contain the definitions of each category.

* `left`: Set of upstream bases: [ACGT]+

* `from`: Set of original (reference) bases: [ACGT]+

* `change`: Base substitution (transition, flip transversion, skew
  transversion): [tfs]

* `right`: Set of downstream base: [ACGT]+

* `n`: Number of mutations in this category

* `N`: Total number of genomic positions in this category times number of possible
  substitutions assigned to this category

* `rate`: Mutation rate for this category (inferred from binomial MLE for n 
  successes and N trials)

* `ci_low`: Lower bound of 95% confidence interval for rate parameter estimation

* `ci_high`: Upper bound of 95% confidence interval for rate parameter estimation

* `relrate`: N\*rate/n

* `autoname`: Used internally for parsing category names

* `name`: Human-readable name for this category

* `type`: Reserved for future use 

### `per_gene.mutation_counts.txt`

Tab-delimited file summarizing each gene's mutation counts on a per-patient
basis, one gene per line.  For each gene, fields are as follows:

* `name`: HUGO symbol

* `cov_gidx`: Numerical index for each gene

* `longname`: HUGO description

* `chr`: Chromosome

* `start`: Transcript start

* `end`: Transcript end

* `len`: Coding length

* `gc`: GC content

and then one column per patient, displaying the total number of mutations in
this gene for that patient.

### `sample_sig_gene_table.txt`

The same as per_gene.mutation_counts.txt, but only for genes with q-value <= 0.1.

### `results.mat`

MATLAB binary version of sig_genes.txt, containing much more detail.  It can be
loaded into MATLAB by typing
```matlab
  load('results.mat', 'G')
```

## MutSig Configuration <a name="config"></a>

MutSig's algorithm and run parameters are configured via a two column, tab-
delimited text file.  Here is a list of all available parameters, possible
options, and default values.  A sample parameters file (set to defaults) can
be found in test/input/params.txt.

* `number_of_categories_to_discover`: Integer.  Number of SNV categories to discover
  (see "mutcategs.txt" section above for details.) Default: 5

* `skip_permutations`: Boolean.  If true, MutSig will not perform
  clustering/functional significance tests; pCL and pFN will be set to NaN.
  Default: false

* `maxperm`: Integer.  Specifies maximum number of permutations to perform for
  pCL/pFN permutation tests.  Default: 1e5

* `remove_duplicate_patients`: Boolean.  Finds and removes duplicate patients in the
  cohort by comparing mutation overlap between patients.  This should be
  disabled when high levels of overlap are expected between samples (e.g.
  primary/met combined cohorts).  Default: true
