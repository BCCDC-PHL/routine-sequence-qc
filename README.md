# Routine Sequence QC
A generic pipeline that can be run routinely on all Illumina sequence runs, regardless of the project or organism of interest.

* Sequence quality information
* Possible contamination

## Usage

```
nextflow run BCCDC-PHL/routine-sequence-qc-nf --run_dir <your illumina run directory> --outdir <output directory>
```