# Routine Sequence QC

![Tests](https://github.com/BCCDC-PHL/routine-sequence-qc/actions/workflows/push_main.yml/badge.svg)

A generic pipeline that can be run routinely on all Illumina sequence runs, regardless of the project or organism of interest.

* Sequence quality information
* Possible contamination

## Analyses

* Parse run-level QC statistics from the 'InterOp' directory and write to `.csv` and `.json` format.
* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/): sample-level sequence quality metrics
* [`mash`](https://github.com/marbl/Mash): Estimate genome size and depth of coverage
* [`seqtk fqchk`](https://github.com/lh3/seqtk): Measure average sequence quality score, percent of bases above Q30 and GC content.
* [`kraken2`](https://github.com/DerrickWood/kraken2) + [`bracken`](https://github.com/jenniferlu717/Bracken): Taxonomic classification
of reads. Estimation of relative abundances of taxonomic groups (genus, species) in each sample.
* [MultiQC](https://github.com/ewels/MultiQC): Collect several QC metrics into a single interactive HTML report.

## Usage

```
nextflow run BCCDC-PHL/routine-sequence-qc \
  [--instrument_type nextseq] \
  [--kraken2_db /path/to/kraken2_db] \
  [--bracken_db /path/to/bracken_db] \
  [--seqtk_fqchk_threshold <threshold>] \
  [--mash_sketch_kmer_size <kmer_size>] \
  [--mash_sketch_minimum_copies <copies>] \
  --run_dir <your illumina run directory> \
  --outdir <output directory>
```

## Outputs

```
<outdir>
├── abundance_top_n
│   ├── top_3_abundances_genus.csv
│   └── top_5_abundances_species.csv
├── basic_qc_stats
│   └── basic_qc_stats.csv
├── bracken
│   ├── <sample_id>_Genus_bracken_abundances.tsv
│   ├── <sample_id>_Genus_bracken_abundances_adjusted.tsv
│   ├── <sample_id>_Genus_bracken.txt
│   ├── <sample_id>_Genus_bracken_adjusted.txt
│   ├── <sample_id>_Species_bracken_abundances.tsv
│   ├── <sample_id>_Species_bracken_abundances_adjusted.tsv
│   ├── <sample_id>_Species_bracken.txt
│   ├── <sample_id>_Species_bracken_adjusted.txt
│   ├── ...
├── fastqc
│   ├── <sample_id>_R1_fastqc
│   ├── <sample_id>_R2_fastqc
│   ├── ...
├── interop_summary
│   ├── interop_index-summary.csv
│   ├── interop_summary.csv
│   └── interop_summary.json
├── kraken2
│   ├── <sample_id>_kraken2.txt
│   ├── ...
├── mash_sketch
│   ├── <sample_id>_R1_mash_sketch.txt
├── mash_sketch_summary
│   └── mash_sketch_summary.csv
├── multiqc
│   ├── multiqc_data
│   └── multiqc_report.html
└── parse_sample_sheet
│   └── sample_sheet.json
├── pipeline_complete.json
├── seqtk_fqchk
│   ├── <sample_id>_seqtk_fqchk_all_positions.csv
│   ├── <sample_id>_seqtk_fqchk_by_position.csv
└── seqtk_fqchk_summary
    └── seqtk_fqchk_summary.csv
```
