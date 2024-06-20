# deepvariant_nf

This repository contains a Nextflow variant calling pipeline for analyzing Next-Generation Sequencing (NGS) data using [deepvariant](https://github.com/google/deepvariant).

```mermaid
%%{init: {'theme':'base'}}%%
flowchart TB
    subgraph "input"
    v0["reads"]
    v1["reference"]
    end
    v3([QCONTROL])
    subgraph "output"
    v18["vcf"]
    v25["html"]
    end
    v5([TRIM])
    v8([ALIGN])
    v9([FLAGSTAT])
    v10([QUALIMAP])
    v12([FAINDEX])
    v13([BAMINDEX])
    v14([VARCALL])
    v17([ANNOTATE])
    v24([REPORT])
    v2([REFINDEX])
    v19(( ))
    v20(( ))
    v21(( ))
    v22(( ))
    v23(( ))
    v0 --> v3
    v0 --> v5
    v1 --> v2
    v3 --> v20
    v5 --> v8
    v5 --> v19
    v1 --> v8
    v2 --> v8
    v8 --> v9
    v8 --> v10
    v8 --> v13
    v8 --> v14
    v9 --> v21
    v10 --> v22
    v1 --> v12
    v12 --> v14
    v13 --> v14
    v1 --> v14
    v14 --> v17
    v17 --> v18
    v17 --> v23
    v19 --> v24
    v20 --> v24
    v21 --> v24
    v22 --> v24
    v23 --> v24
    v24 --> v25
```

## Description

The pipeline is implemented in Nextflow and includes several stages for NGS data analysis:

1. **REFINDEX:** Index creation using BWA (Burrows-Wheeler Aligner).
2. **QCONTROL:** Quality control using FastQC.
3. **TRIM:** Data preprocessing using Fastp.
4. **ALIGN:** Sequence alignment using BWA mem.
5. **FLAGSTAT:** Alignment quality control using Samtools flagstat.
6. **QUALIMAP:** Alignment quality control using Qualimap.
7. **FAINDEX:** Fai index creation using Samtools faidx.
8. **BAMINDEX:** Bai index creation using Samtools index.
9. **PREPARE:** File processing and preparation using Samtools.
10. **VARCALL:** Variant calling using deepvariant.
11. **ANNOTATE:** Annotation using VEP (Variant Effect Predictor).
12. **REPORT:** Compiling report using MultiQC.

## Usage

### Quick Start

To quickly run the pipeline, use the following command:

```bash
nextflow run glebus-sasha/deepvariant_nf \
-profile <docker/singularity> \
--reference <path-to-reference> \
--reads "<path-to-reads-folder>/*[rR]{1,2}*.fastq*" \ # quotes are important
--outdir results
```

### Requirements

- Nextflow (https://www.nextflow.io/docs/latest/install.html)
- Docker (https://docs.docker.com/engine/install/) or
- Singularity (https://github.com/sylabs/singularity/blob/main/INSTALL.md)

### Running the Pipeline

1. Install all the necessary dependencies such as Nextflow, Singularity.
3. Clone this repository: `git clone https://github.com/glebus-sasha/deepvariant_nf.git`
4. Navigate to the pipeline directory: `cd deepvariant_nf`
5. Edit the `nextflow.config` file to set the required parameters, if necessary.
6. Run the pipeline, setting the required parameters, for example:

```bash
nextflow run main.nf
```

## License

This project is licensed under the [MIT License](LICENSE).
