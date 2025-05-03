
# alndv

A nextflow (DSL 2) Whole‑Genome Short‑Read Pipeline (BWA‑MEM2 → DeepVariant → GLnexus) for small‑variant discovery from paired‑end FASTQ files.


---
## 1. Overview

```

FASTQ ─► FASTQC ─► BWA‑MEM2 ──► MERGEB ─► QUALIMAP 
                                   └─► DEPTH (mosdepth)
                                   └─► DEEPVARIANT_AUTOSOMES ─┬─► per‑sample VCF
                                   |         |                  └─► per‑sample gVCF
                                   |         │
                                   |         ▼
                                   |       GLnexus (cohort merge) ─► joint autosomes BCF
                                   | 
                                   └─► DEEPVARIANT_SEX ─┬─► per‑sample VCF
                                            │            └─► per‑sample gVCF
                                            │
                                            ▼
                                          GLnexus (cohort merge) ─► joint SEX BCF

```

| Step | Tool | Purpose |
|------|------|---------|
| **FASTQC** | FastQC 0.12 | Raw‑read quality control |
| **BWAMEM** | `bwa‑mem2` 2.3 / `bwa` 0.7 | Alignment, alt‑aware mapping, optional HLA typing |
| **MERGEB** | samtools 1.19 | Merge technical replicates into a single CRAM |
| **QUALIMAP** | Qualimap 2.3 | Alignment QC (streaming from CRAM via FIFO) |
| **DEPTH** | mosdepth 0.3 | Depth & breadth of coverage |
| **DEEPVARIANT** | DeepVariant 1.8 | SNP/indel calling (WGS model) |
| **GLNEXUS_DEEPVARIANT** | GLnexus 1.4 | Joint genotyping of all gVCFs |

---

## 2. File‑tree after a successful run (default `--outdir results`)

```

results/
├── QC/
│   └── FASTQC/         # \*.html and \*.zip reports
├── BWA/
│   ├── \*.log.bwamem
│   └── HLA/            # .hla.all and logs (if --hla)
├── CRAM/
│   ├── \*.cram
│   └── \*.cram.crai    
├── DEPTH/
│   ├── \*.mosdepth.dist.txt
│   └── \*.mosdepth.summary.txt
├── qualimap/
│   └── SAMPLE/         # Qualimap HTML + tables
├── deepvariant\_persample/
│   ├── \*.vcf.gz
│   └── \*.g.vcf.gz
└── glx/
└── all.glx.bcf

```

---

## 3. Quick start

```bash
nextflow run main.nf --csv reads-test.cvs --ref test.fa --debug true
```

### Required inputs

| Flag      | Description                                                  |
| --------- | ------------------------------------------------------------ |
| `--reads` | Glob pattern for paired FASTQ (`sample_*{1,2}.fq.gz`) **or** |
| `--csv`   | CSV with columns `sampleId,part,read1,read2`                 |
| `--ref`   | Reference FASTA indexed with `.fai` and `.dict`              |

---

## 4. Key parameters

| Parameter   | Default    | Notes                                                  |
| ----------- | ---------- | ------------------------------------------------------                       |
| `--alt`     | `false`    | Activate *alt‑aware* alignment (`*.alt` file required) |
| `--hla`     | `false`    | Run **xHLA / k8 Alt-aware** typing (human data only)   |
| `--debug`   | `false`    | Dry‑run style outputs; generates dummy files           |

---

## 5. Process outputs

| Process               | Channel             | Data type | Path / pattern                         |
| --------------------- | ------------------- | --------- | -------------------------------------- |
| `FASTQC`              | `fqc`               | directory | `<sample>-<part>.fastqc/`              |
| `BWAMEM`              | `bams`              | tuple     | `<sample>-<part>.mkdup.cram` + `.crai` |
| `MERGEB`              | `mbams`             | tuple     | `<sample>.merged.cram` + `.crai`       |
| `QUALIMAP`            | `qualimap_results`  | dir       | `qualimap/<sample>/`                   |
| `DEPTH`               | `depth`             | files     | `<sample>.depth.*`                     |
| `DEEPVARIANT`         | `vcf` / `gvcf_file` | files     | `<sample>.vcf.gz`, `<sample>.g.vcf.gz` |
| `GLNEXUS_DEEPVARIANT` | `bcf`               | file      | `all.glx.bcf`                          |

---

## 6. Joint‑calling logic

All per‑sample gVCFs are collected into **one** tuple:

```nextflow
allgvcf = DEEPVARIANT.out.gvcf_file.collect()
          .map{ gvcfs -> tuple('all', gvcfs) }
```

The tuple `('all', [list_of_gvcfs])` feeds `GLNEXUS_DEEPVARIANT`, which emits
`all.glx.bcf` (Bgzipped BCF, ready for analysis or left‑alignment/annotation).

---

## 7. Containers & reproducibility

Every process runs inside an OCI/Singularity container pinned by digest:

* `bwa-mem2_instrain_multiqc_qualimap_samtools:850f96...`
* `google/deepvariant:1.8.0`
* `biocontainers/glnexus:1.4.1--h40d77a6_0`

Pull happens automatically via Nextflow Wave/Seqera Registry when `-with-docker` or `-with-singularity` is enabled.

---

## 8. Troubleshooting

| Issue                       | Remedy                                                                             |
| --------------------------- | ---------------------------------------------------------------------------------- |
| `qualimap` hangs            | Ensure JVM has enough heap (`--java-mem-size`) and CRAM reference matches `--ref`. |
| DeepVariant “invalid range” | Reference/CRAM mismatch or truncated CRAM; regenerate index.                       |
| GLnexus OOM                 | Lower `--threads` or set `NX_MEM_GB` via `params.memory`.                          |

---

## 9. Developer

1. Alex Di Genova

## 10. License

Licensed under the MIT License.  
Pull requests welcome!
