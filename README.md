# Routine Sequence QC
A generic pipeline that can be run routinely on all Illumina sequence runs, regardless of the project or organism of interest.

* Sequence quality information
* Possible contamination

## Analyses

* Parse run-level QC statistics from the 'InterOp' directory and write to `.csv` and `.json` format.
* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/): sample-level sequence quality metrics
* [Kraken2](https://github.com/DerrickWood/kraken2) + [Bracken](https://github.com/jenniferlu717/Bracken): Taxonomic classification
of reads. Estimation of relative abundances of taxonomic groups (genus, species) in each sample.
* [MultiQC](https://github.com/ewels/MultiQC): Collect several QC metrics into a single interactive HTML report.

## Usage

```
nextflow run BCCDC-PHL/routine-sequence-qc-nf --run_dir <your illumina run directory> --outdir <output directory>
```

## Outputs

```
<outdir>
├── abundance_top_n
│   ├── top_3_abundances_genus.csv
│   └── top_5_abundances_species.csv
├── bracken
│   ├── <sample_id>_Genus_bracken_abundances.tsv
│   ├── <sample_id>_Genus_bracken.txt
│   ├── <sample_id>_Species_bracken_abundances.tsv
│   ├── <sample_id>_Species_bracken.txt
│   ├── ...
├── fastqc
│   ├── <sample_id>_R1_fastqc
│   ├── <sample_id>_R2_fastqc
│   ├── ...
├── interop_summary
│   ├── interop_index-summary.csv
│   └── interop_summary.csv
├── kraken2
│   ├── <sample_id>_kraken2.txt
│   ├── ...
├── multiqc
│   ├── multiqc_data
│   └── multiqc_report.html
└── parse_sample_sheet
    └── sample_sheet.json
```
