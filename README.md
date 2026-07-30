# Table of Contents
- [MetaPrism](#metaprism)
- [Requirements](#requirements)
- [Install](#install)
- [Tutorial](#tutorial)
- [Usages](#usages)
- [Citation](#citation)

# MetaPrism
MetaPrism: A Toolkit for Joint Analysis of Meta-genomic Sequencing Data 

MetaPrism provides joint profile (infer both taxonomical and functional profile) for shotgun metagenomic sequencing data. It also offer tools to 1) classify sequence reads and estimate the abundances for taxa-specific genes; 2) tabularize and visualize taxa-specific gene abundances; 3) build asso-ciation and prediction models for comparative analysis. 


## Requirements

1. Perl - https://www.perl.org
2. R - http://www.r-project.org
3. Perl module Statistics::R - https://metacpan.org/pod/Statistics::R
4. R library caret - https://cran.r-project.org/web/packages/caret/index.html
5. R library randomForest - https://cran.r-project.org/web/packages/randomForest/index.html
6. DIAMOND - https://github.com/bbuchfink/diamond (recommended) or USEARCH - https://www.drive5.com/usearch/
7. BWA - http://bio-bwa.sourceforge.net
8. Samtools - http://www.htslib.org
9. Centrifuge - https://ccb.jhu.edu/software/centrifuge/
10. Linux commands: sort, wget - https://www.gnu.org/software/wget/
11. MEGAHIT - https://github.com/voutcn/megahit (optional) 


## Install

If you already have Git (https://git-scm.com) installed, you can get the latest development version using Git. It will take a few seconds.
```
git clone https://github.com/jiwoongbio/MetaPrism.git
```


## Tutorial

### Current assembly-annotation and abundance workflow

The current MetaPrism workflow separates gene annotation from read-depth quantification. `MetaPrism.pl` identifies and annotates genes in an assembled genome or metagenome, and `MetaPrism.abundance.pl` estimates the abundance of the annotated regions from sequencing reads mapped back to the same assembly. The optional reference-assembly preparation adds taxonomic assignments and creates a ranked single-copy marker catalog for genome-normalized abundance estimation.

A complete paired-end short-read example is shown below:

```bash
git clone https://github.com/jiwoongbio/MetaPrism.git
cd MetaPrism

threads=16

perl MetaPrism_data.gene.pl -p $threads
perl MetaPrism_data.taxon.pl
perl MetaPrism_data.gene.base.pl -p $threads

perl MetaPrism.pl -p $threads genome_assembly.fasta > MetaPrism.txt

bwa index genome_assembly.fasta
bwa mem -t $threads genome_assembly.fasta sample.1.fastq.gz sample.2.fastq.gz \
    | samtools sort - > sample.sorted.bam
samtools index sample.sorted.bam

perl MetaPrism.abundance.pl -p $threads \
    MetaPrism.txt sample.sorted.bam \
    > MetaPrism.abundance.txt
```

`MetaPrism_data.taxon.pl` and `MetaPrism_data.gene.base.pl` are not required for gene annotation alone. They are included in the complete workflow because they enable automatic taxonomic annotation and generate a marker catalog tailored to the reference assemblies used by MetaPrism. If the marker catalog is not available, `MetaPrism.abundance.pl` can still normalize abundance using its built-in set of ten base-abundance genes.

#### Prepare the gene database

```bash
perl MetaPrism_data.gene.pl -p $threads
```

With no database-building option and no `-s` option, `MetaPrism_data.gene.pl` downloads the prebuilt gene database currently designated as the default by the MetaPrism data server. It then creates the local DIAMOND database used by `MetaPrism.pl`. In this command, `-p $threads` specifies the number of threads used by `diamond makedb`; it does not control the downloads.

The optional `-s` argument selects an available orthology database by the beginning of its build timestamp. For example, the following command requests the database version whose timestamp begins with `2024`:

```bash
perl MetaPrism_data.gene.pl -p $threads -s 2024
```

The normal workflow omits `-s` and installs the current default database.

#### Prepare the taxonomy reference and single-copy marker catalog

```bash
perl MetaPrism_data.taxon.pl
perl MetaPrism_data.gene.base.pl -p $threads
```

`MetaPrism_data.taxon.pl` selects prioritized RefSeq and GenBank assemblies from the configured taxonomic groups. By default, the selected groups are Bacteria, Archaea, and Fungi. It combines the selected assemblies in `data/species_assembly.fasta` and builds the minimap2 index `data/species_assembly.mni`.

When `data/species_assembly.mni` is present, `MetaPrism.pl` detects it automatically and uses it to add taxonomic assignments to gene annotations. The same reference assembly collection is analyzed by `MetaPrism_data.gene.base.pl` to construct a marker catalog for abundance normalization.

`MetaPrism_data.gene.base.pl` annotates each reference assembly independently with the active MetaPrism gene database. A gene contributes one count for an assembly only when it is detected exactly once in that assembly. The script then ranks genes by the number of supporting assemblies and writes the result to:

```text
data/single_copy_gene.count.txt
```

This ranked catalog provides broadly supported single-copy markers for estimating a sample-level abundance baseline. For `MetaPrism_data.gene.base.pl`, `-p $threads` is the maximum number of reference assemblies analyzed concurrently. Each assembly-level job runs `MetaPrism.pl` with one thread. Completed assembly annotations are retained under `data/species_assembly.MetaPrism/` and reused when the command is run again.

#### Annotate genes in an assembly

```bash
perl MetaPrism.pl -p $threads genome_assembly.fasta > MetaPrism.txt
```

`genome_assembly.fasta` may contain one or more contigs and may be uncompressed or gzip-compressed. For nucleotide input, `MetaPrism.pl`:

1. examines both strands in all three reading frames and identifies candidate coding regions from the configured start and termination codons;
2. translates the candidate regions using NCBI genetic code 11 by default;
3. searches the translated sequences against the active gene database with DIAMOND `blastp`;
4. filters mappings by E-value, alignment coverage, and any database-specific identity threshold; and
5. reports the coding coordinates, gene assignment, matched reference protein, and protein-sequence differences.

For this command, `-p $threads` controls parallel sequence translation and the DIAMOND search. When taxonomic matching is performed, the same thread count is passed to minimap2.

With the standard downloaded database, the annotation table contains these columns:

```text
chromosome  start  end  strand  gene  definition  protein  variant
```

`chromosome`, `start`, `end`, and `strand` identify the annotated region in the input assembly. `gene` and `definition` describe the assignment, `protein` is the matched reference-protein identifier, and `variant` reports comma-separated differences in `reference_position|reference|query` format.

When `data/species_assembly.mni` is available, `MetaPrism.pl` also appends:

```text
taxonId  taxonName  taxonRank
```

The default taxonomic rank is `species`; another rank can be requested with `-r`. A protein FASTA file can be analyzed directly with `-P`, in which case translation is skipped and an `input` column replaces the genomic-coordinate columns.

#### Map sequencing reads to the annotated assembly

`MetaPrism.abundance.pl` requires one or more coordinate-sorted BAM files whose reference names and coordinates correspond to the assembly used by `MetaPrism.pl`. The script creates a BAM index when one is not already available. Any suitable aligner may be used. In the example below, `samtools sort -` reads the BWA alignment stream from standard input, and shell redirection writes the coordinate-sorted BAM to `sample.sorted.bam`:

```bash
bwa index genome_assembly.fasta
bwa mem -t $threads genome_assembly.fasta sample.1.fastq.gz sample.2.fastq.gz \
    | samtools sort - > sample.sorted.bam
samtools index sample.sorted.bam
```

When data from multiple lanes or read sets have been aligned to the same assembly, their coordinate-sorted BAM files may be supplied together to `MetaPrism.abundance.pl`. Their depth contributions are added during abundance calculation.

#### Calculate and normalize gene abundance

```bash
perl MetaPrism.abundance.pl -p $threads \
    MetaPrism.txt sample.sorted.bam \
    > MetaPrism.abundance.txt
```

For every row in `MetaPrism.txt`, `MetaPrism.abundance.pl` calculates the mean mapped-read depth across the annotated interval. It preserves all annotation columns and appends an `abundance` column.

When `-b` is omitted, the script chooses the default base-abundance genes as follows:

1. If `data/single_copy_gene.count.txt` is readable and nonempty, it selects the first 10 genes from that ranked catalog.
2. Otherwise, it uses the built-in gene set `K02950,K02874,K02946,K02948,K02867,K02952,K02886,K02988,K02992,K02965`.

The normalization then proceeds as follows:

1. Mean depth is calculated separately for every annotated interval.
2. Interval depths assigned to the same gene are summed to obtain one abundance value per selected marker gene. A selected marker that is not detected contributes zero.
3. The marker abundances are sorted, and the two lowest and two highest values are removed by default.
4. The remaining values are averaged to obtain the sample's base abundance.
5. The mean depth of every annotated interval is divided by the base abundance.

When the ranked marker catalog has been prepared, this procedure estimates the abundance baseline from genes that are repeatedly observed as single copy across the reference assembly collection. It reduces the direct effect of sequencing depth and total microbial DNA on the reported values and makes gene abundance more interpretable relative to the microbial genome content represented in the sample. A value near 1 indicates depth close to the estimated marker baseline; it should be interpreted as a relative genome-normalized abundance rather than an exact integer copy number.

The marker set can be changed with `-b`. An integer selects that many genes from the top of `data/single_copy_gene.count.txt`; a comma-separated string specifies gene identifiers directly; and a readable file supplies one gene identifier per line. The number of marker values removed from each end before calculating the mean is controlled by `-e`.

For example, this command uses the first 20 genes in the ranked catalog and removes the four lowest and four highest marker abundances:

```bash
perl MetaPrism.abundance.pl -p $threads -b 20 -e 4 \
    MetaPrism.txt sample.sorted.bam \
    > MetaPrism.abundance.txt
```

### Original MetaPrism workflow

We present a short tutorial to help users quickly get started on their own analysis. The example datasets are based on [https://www.ncbi.nlm.nih.gov/bioproject/PRJNA397906](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA397906) and the full command list is  available at [example/example.sh](example/example.sh).

1. Prepare database files

This step will download necessary databases for MetaPrism in your current directory.

```
perl MetaPrism_gene_prepare.pl
```

2. De novo metagenome assembly (per sample)

For each sample, users need to perform de novo assembly. Suppose the input paired end sequence files are `sample1.1.fastq.gz` and `sample1.2.fastq.gz`. The commnad line is:

```
megahit -1 sample1.1.fastq.gz -2 sample1.2.fastq.gz -o sample1.megahit
```

3. Gene annotation and abundance quantification (per sample)

For each sample, MetaPrism first quantify all gene abundances. The file `sample1.megahit/final.contigs.fa` is the output from the previous step.

```
perl MetaPrism_gene.pl sample1.gene sample1.megahit/final.contigs.fa sample1.1.fastq.gz,sample1.2.fastq.gz
```

4. Taxon annotation (per sample)

MetaPrism will next infer the taxonomy for each contig. The file `sample1.gene.region.abundance.txt` is the output from the previous step. The result will be outputted to the console. We redirected it to the result file `sample1.gene_taxon.region.abundance.txt`.

```
perl MetaPrism_taxon_centrifuge.pl sample1.gene.region.abundance.txt sample1.megahit/final.contigs.fa centrifuge/data/p_compressed > sample1.gene_taxon.region.abundance.txt
```

5. Compare sample groups and identify differentially-abundant genes

Suppose that you repeat step 2 to step 4, and you get a list of joint features (`sample1.gene_taxon.region.abundance.txt`, `sample2.gene_taxon.region.abundance.txt`, ..., `sample6.gene_taxon.region.abundance.txt`). MetaPrism can perform comparative analysis using the following command:

```
perl MetaPrism_comparison.pl -F gene sample.group.txt \
	sample1=sample1.gene_taxon.region.abundance.txt \
	sample2=sample2.gene_taxon.region.abundance.txt \
	sample3=sample3.gene_taxon.region.abundance.txt \
	sample4=sample4.gene_taxon.region.abundance.txt \
	sample5=sample5.gene_taxon.region.abundance.txt \
	sample6=sample6.gene_taxon.region.abundance.txt \
	> gene.comparison.txt

awk -F'\t' '(NR == 1 || ($4 >= 1 && $5 <= 0.01))' gene.comparison.txt > gene.comparison.filtered.txt
```

* sample.group.txt is a text file containing lines of tab-delimited sample and group like following:

````
  sample1 group1
  sample2 group1
  sample3 group1
  sample4 group2
  sample5 group2
  sample6 group2
````

Here `sample1`, `sample2`, and `sample3` are from `group1`, and the rest are from `group2`. For another example group file, see [example/sample.group.txt](example/sample.group.txt).

6. Generate a heatmap webpage

You can also generate a heatmap webpage using the `MetaPrism_heatmap.pl` command.

```
perl MetaPrism_heatmap.pl -F gene -s -g gene.comparison.filtered.txt -r both \
	sample1=sample1.gene_taxon.region.abundance.txt \
	sample2=sample2.gene_taxon.region.abundance.txt \
	sample3=sample3.gene_taxon.region.abundance.txt \
	sample4=sample4.gene_taxon.region.abundance.txt \
	sample5=sample5.gene_taxon.region.abundance.txt \
	sample6=sample6.gene_taxon.region.abundance.txt \
	> gene.heatmap.html
```

7. Generate a tabular result file

Users may also want to have a tabular file for their own analysis. This command will produce such tabular text file: 

```
perl MetaPrism_table.pl -F taxon_average -s \
	sample1=sample1.gene_taxon.region.abundance.txt \
	sample2=sample2.gene_taxon.region.abundance.txt \
	sample3=sample3.gene_taxon.region.abundance.txt \
	sample4=sample4.gene_taxon.region.abundance.txt \
	sample5=sample5.gene_taxon.region.abundance.txt \
	sample6=sample6.gene_taxon.region.abundance.txt \
	> taxon.table.txt
```

8. Build a prediction model

Users can build a prediction model using `MetaPrism_prediction.pl`. The option `-t xgbTree` uses the xgboost algorithm with leave-one-out cross validation. The result file `prediction.feature.txt` lists the feature importances, and the file `prediction.txt` lists the prediction accuracies.

```
perl MetaPrism_prediction.pl -t xgbTree -f prediction.feature.txt sample.group.txt \
	sample1=sample1.gene_taxon.region.abundance.txt \
	sample2=sample2.gene_taxon.region.abundance.txt \
	sample3=sample3.gene_taxon.region.abundance.txt \
	sample4=sample4.gene_taxon.region.abundance.txt \
	sample5=sample5.gene_taxon.region.abundance.txt \
	sample6=sample6.gene_taxon.region.abundance.txt \
	> prediction.txt
```


## Usages

The current assembly-annotation and abundance commands are listed first, followed by the commands from the original MetaPrism workflow.

### Current assembly-annotation and abundance commands

* `MetaPrism_data.gene.pl`

Downloads a prebuilt gene database or builds a new orthology or antimicrobial-resistance gene database. With neither `-b` nor `-a`, it installs a prebuilt database. Omitting `-s` installs the database currently designated as the default.

```text
Usage:   perl MetaPrism_data.gene.pl [options]

Options: -h       display this help message
         -p INT   number of threads used by diamond makedb [1]
         -r       redownload data
         -b       build orthology database
         -a       build AMR gene database
         -t STR   comma-separated NCBI taxonomy IDs [2,2157,4751]
         -i INT   UniRef identity: 50, 90, or 100 [100]
         -s STR   select a prebuilt orthology database by build-time prefix
```

The `-t` and `-i` values define the orthology-database family when a database is built with `-b` or an archived version is selected with `-s`. The `-b` and `-a` modes are database-maintenance functions and are not required for the normal prebuilt-database workflow.

* `MetaPrism_data.taxon.pl`

Builds the reference assembly collection used for taxonomic matching and single-copy marker discovery. It selects prioritized RefSeq and GenBank assemblies within the taxonomic groups specified by `-t`, writes the combined sequences to `data/species_assembly.fasta`, and builds `data/species_assembly.mni` for minimap2.

```text
Usage:   perl MetaPrism_data.taxon.pl [options]

Options: -h       display this help message
         -r       redownload data
         -t STR   comma-separated NCBI taxonomy IDs [2,2157,4751]
```

The default taxonomy IDs cover Bacteria, Archaea, and Fungi. Run this command before `MetaPrism_data.gene.base.pl`. Once `data/species_assembly.mni` is available, `MetaPrism.pl` detects it automatically and adds taxonomic fields to its output.

* `MetaPrism_data.gene.base.pl`

Builds the ranked single-copy marker catalog from the assemblies prepared by `MetaPrism_data.taxon.pl`. Each reference assembly is annotated independently. A gene receives one count for an assembly only when it occurs exactly once, and the result is written to `data/single_copy_gene.count.txt` as the gene identifier followed by the number of supporting assemblies. The file is sorted with the most broadly supported genes first. When this file is readable and nonempty, `MetaPrism.abundance.pl` uses its first 10 genes as the default normalization markers instead of the built-in marker set.

```text
Usage:   perl MetaPrism_data.gene.base.pl [options]

Options: -h       display this help message
         -t DIR   directory for temporary files [$TMPDIR or /tmp]
         -p INT   maximum number of assemblies processed concurrently [1]
```

Each assembly-level child process runs `MetaPrism.pl` with one thread. Completed results are cached under `data/species_assembly.MetaPrism/` and reused on subsequent runs.

* `MetaPrism.pl`

Annotates genes and protein variants from a nucleotide assembly FASTA file or, with `-P`, from a protein FASTA file. Results are written to standard output.

```text
Usage:   perl MetaPrism.pl [options] genome.fasta > MetaPrism.txt

Options: -h       display this help message
         -t DIR   directory for temporary files [$TMPDIR or /tmp]
         -p INT   number of threads [1]
         -P       input is protein fasta file
         -C STR   codon and translation e.g. ATG=M [NCBI genetic code 11 (Bacterial, Archaeal and Plant Plastid)]
         -S STR   comma-separated start codons [GTG,ATG,CTG,TTG,ATA,ATC,ATT]
         -T STR   comma-separated termination codons [TAG,TAA,TGA]
         -l INT   minimum translation length [10]
         -e FLOAT maximum e-value to report alignments [10]
         -k INT   maximum target [100]
         -c FLOAT minimum coverage [0.8]
         -d FILE  diamond path [diamond]
         -D FILE  diamond database file
         -G FILE  gene definition file
         -g FILE  protein-to-gene file
         -m FILE  minimap2 file
         -r STR   taxonomic rank [species]
```

Unless custom files are supplied with `-D`, `-G`, or `-g`, `MetaPrism.pl` uses the active database prepared by `MetaPrism_data.gene.pl`. For nucleotide input, `-p` controls parallel translation and the DIAMOND search. It is also passed to minimap2 when taxonomic matching is performed. The resulting annotation table is the first positional input to `MetaPrism.abundance.pl`.

* `MetaPrism.abundance.pl`

Calculates the mean mapped-read depth of each genomic region in a nucleotide-mode `MetaPrism.pl` annotation table and normalizes that depth to a base abundance. The input BAM files must be coordinate-sorted and aligned against the same assembly used to create `MetaPrism.txt`. Multiple BAM files may be supplied; their depth contributions are summed. All input annotation columns are retained, and an `abundance` column is appended.

```text
Usage:   perl MetaPrism.abundance.pl [options] MetaPrism.txt sample.sorted.bam [...] > MetaPrism.abundance.txt

Options: -h       display this help message
         -t DIR   directory for temporary files [$TMPDIR or /tmp]
         -p INT   number of parallel worker processes [1]
         -q INT   minimum mapping quality [0]
         -f INT   require these SAM flag bits [0]
         -F INT   exclude reads with these SAM flag bits [0]
         -S STR   stranded library orientation, "f" or "r"
         -B FLOAT use a specified base abundance
         -b STR   base-abundance genes or number of catalog genes [built-in 10-gene set, or first 10 catalog genes when available]
         -e INT   number of low and high marker values to remove [2]
         --use_merged_bam_file  use a temporary merged BAM file
         --use_samtools_depth   use samtools depth instead of CIGAR parsing
```

`-p` controls the number of annotation batches processed concurrently and is also passed to `samtools merge` when temporary BAM merging is requested.

If `-b` is omitted, the default marker set is resolved automatically:

- when `data/single_copy_gene.count.txt` is readable and nonempty, the first 10 genes in that ranked catalog are used;
- otherwise, the built-in genes `K02950,K02874,K02946,K02948,K02867,K02952,K02886,K02988,K02992,K02965` are used.

For each selected marker, the script sums the mean depths of all annotation intervals assigned to that gene. Missing markers contribute zero. It then removes the lowest and highest marker values according to `-e` (two from each side by default), averages the remaining values, and divides every interval's mean depth by this base abundance.

The `-b` argument accepts three forms:

- an integer, which selects that many genes from the top of `data/single_copy_gene.count.txt`;
- a comma-separated list of gene identifiers; or
- a readable text file containing one gene identifier per line.

An integer `-b` value requires `data/single_copy_gene.count.txt`. Supplying an empty `-b` value does not disable normalization in this version; it activates the same automatic default-marker selection used when `-b` is omitted. When marker genes are selected, their calculated mean replaces any value supplied with `-B`. The mapping-quality, SAM-flag, and strandedness filters are applied before depth is calculated.

If the optional `Bio::DB::Sam` Perl module is available, MetaPrism can calculate interval coverage through that module; otherwise, it falls back to samtools. `--use_samtools_depth` explicitly selects the temporary-BAM and `samtools depth` implementation.

### Original MetaPrism commands


In this section, we list the command line option for all available `MetaPrism` functions.

* MetaPrism_gene_prepare.pl

```
Usage:   perl MetaPrism_gene_preapare.pl [options]

Options: -h       display this help message
         -r       redownload data
         -m FILE  executable file path of mapping program, "diamond" or "usearch" [diamond]
         -k       download prebuilt KEGG files
         -a       download ARDB database
         -b       download beta-lactamase database
```

* MetaPrism_gene.pl

```
Usage:   perl MetaPrism_gene.pl [options] output.prefix genome.fasta [input.fastq|input.R1.fastq,input.R2.fastq [...]]

Options: -h       display this help message
         -A STR   prepared genome prefix
         -B       input indexed sorted BAM file instead of FASTQ file
         -m FILE  executable file path of mapping program, "diamond" or "usearch" [diamond]
         -p INT   number of threads [1]
         -e FLOAT maximum e-value to report alignments [10]
         -t DIR   directory for temporary files [$TMPDIR or /tmp]
         -a FLOAT search acceleration for ublast [0.5]
         -C STR   codon and translation e.g. ATG=M [NCBI genetic code 11 (Bacterial, Archaeal and Plant Plastid)]
         -S STR   comma-separated start codons [GTG,ATG,CTG,TTG,ATA,ATC,ATT]
         -T STR   comma-separated termination codons [TAG,TAA,TGA]
         -l INT   minimum translation length [10]
         -c FLOAT minimum coverage [0.8]
         -q INT   minimum mapping quality [0]
         -s STR   strand specificity, "f" or "r"
         -P STR   contig prefix used for abundance estimation
```

* MetaPrism_taxon_centrifuge.pl

```
Usage:   perl MetaPrism_taxon_centrifuge.pl [options] MetaPrism_gene.region.txt genome.fasta centrifuge.index > MetaPrism_gene.region.taxon.txt

Options: -h       display this help message
         -p INT   number of threads [1]
```

* MetaPrism_comparison.pl

```
Usage:   perl MetaPrism_comparison.pl [options] sample.group.txt [sample=]abundance.txt [...] > MetaPrism_comparison.txt

Options: -h       display this help message
         -A STR   abundance column [meanDepth/genome]
         -R STR   taxon rank [genus]
         -F STR   feature type, "gene_taxon", "gene", "gene_average", "taxon", "taxon_average" [gene_taxon]
         -t STR   statistical test for comparing sample groups, "kruskal", "anova", "poisson", "quasipoisson", "metagenomeSeq" [kruskal]
         -o FLOAT offset [1]
```

* MetaPrism_prediction.pl

```
Usage:   perl MetaPrism_prediction.pl [options] sample.group.txt [sample=]abundance.txt [...]

Options: -h       display this help message
         -A STR   abundance column [meanDepth/genome]
         -R STR   taxon rank [genus]
         -F STR   feature type, "gene_taxon", "gene", "gene_average", "taxon", "taxon_average" [gene_taxon]
         -t STR   train method [rf]
         -c STR   train control method [LOOCV]
         -m FILE  model file
         -f FILE  important feature file
         -s INT   seed [1]
```

* MetaPrism_table.pl

```
Usage:   perl MetaPrism_table.pl [options] [sample=]abundance.txt [...] > table.txt

Options: -h       display this help message
         -A STR   abundance column [meanDepth/genome]
         -R STR   taxon rank [genus]
         -F STR   feature type, "gene_taxon", "gene", "gene_average", "taxon", "taxon_average" [gene_taxon]
         -s       scale
```

* MetaPrism_heatmap.pl

```
Usage:   perl MetaPrism_heatmap.pl [options] [sample=]abundance.txt [...] > heatmap.html

Options: -h       display this help message
         -A STR   abundance column [meanDepth/genome]
         -R STR   taxon rank [genus]
         -F STR   feature type, "gene_taxon", "gene", "gene_average", "taxon", "taxon_average" [gene_taxon]
         -s       scale
         -g FILE  feature file
         -t INT   taxon abbreviation length [4]
         -f INT   HTML font size [15]
         -w INT   HTML table cell width [60]
```


## Citation

Kim J, Jiang S, Wang Y, Xiao G, Xie Y, Liu DJ, Li Q, Koh A, Zhan X. MetaPrism: A versatile toolkit for joint taxa/gene analysis of metagenomic sequencing data. G3 (Bethesda). 2021 Apr 15;11(4):jkab046. doi: 10.1093/g3journal/jkab046. PMID: 33713107; PMCID: PMC8049424.
