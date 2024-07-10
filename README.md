# deepvariant_nf

This repository contains a Nextflow variant calling pipeline for analyzing Next-Generation Sequencing (NGS) data using [deepvariant](https://github.com/google/deepvariant).

```mermaid
%%{init: {'theme':'base'}}%%
flowchart TB
    subgraph "input"
    v0["trimmed reads"]
    v1["reference"]
    v30["indices"]
    v31["vep_DB"]
    end
    subgraph "output"
    v18["vcf"]
    v25["html"]
    end
    v8([ALIGN])
    v9([FLAGSTAT])
    v10([QUALIMAP])
    v13([BAMINDEX])
    v14([VARCALL])
    v17([ANNOTATE])
    v24([REPORT])
    v21(( ))
    v22(( ))
    v23(( ))
    v0 --> v8
    v1 --> v8
    v30 --> v8
    v8 --> v9
    v8 --> v10
    v8 --> v13
    v8 --> v14
    v9 --> v21
    v10 --> v22
    v13 --> v14
    v1 --> v14
    v14 --> v17
    v17 --> v18
    v17 --> v23
    v31 --> v14
    v30 --> v14
    v21 --> v24
    v22 --> v24
    v23 --> v24
    v24 --> v25
```


## Description

The pipeline is implemented in Nextflow and includes several stages for NGS data analysis:

1. **ALIGN:** Sequence alignment using BWA mem.
2. **FLAGSTAT:** Alignment quality control using Samtools flagstat.
3. **QUALIMAP:** Alignment quality control using Qualimap.
4. **FAINDEX:** Fai index creation using Samtools faidx.
5. **BAMINDEX:** Bai index creation using Samtools index.
6. **PREPARE:** File processing and preparation using Samtools.
7. **VARCALL:** Variant calling using deepvariant.
8. **ANNOTATE:** Annotation using VEP (Variant Effect Predictor).
9. **REPORT:** Compiling report using MultiQC.

## Usage

### Quick Start

To quickly run the pipeline, use the following command:

```bash
nextflow run glebus-sasha/deepvariant_nf \
-profile <docker/singularity> \
--reference <path-to-reference> \
--reads <path-to-reads-folder> \
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
