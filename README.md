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