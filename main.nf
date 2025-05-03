
process FASTQC {
    tag "$sampleId-mem"
    //label 'process_medium'
    publishDir "$params.outdir/QC/FASTQC", mode: "copy"



    input:
    tuple val(sampleId), val(part), file(read1), file(read2)

    output:
    path("${sampleId}-${part}.fastqc"), emit: fqc

    script:
    if(params.debug == true){
    """
    echo fastqc -o ${sampleId}-${part}.fastqc $read1 $read2
    mkdir -p ${sampleId}-${part}.fastqc
    touch ${sampleId}-${part}.fastqc/report.fastqc

    """
    } else{
    """
    mkdir -p ${sampleId}-${part}.fastqc
    fastqc -t $task.cpus -o ${sampleId}-${part}.fastqc $read1 $read2
    """
    }

}

//we do run bwa-mem2 or bwa mem
process BWAMEM {

    tag "$sampleId-mem"
    //label 'process_high'
    publishDir "$params.outdir/BWA", mode: "copy", pattern: '*.log.*'
    publishDir "$params.outdir/BWA/HLA", mode: "copy", pattern: '*.hla.all'

    input:
    tuple val(sampleId), val(part), file(read1), file(read2)

    output:
    tuple val("${sampleId}"), val("${part}"), file("${sampleId}-${part}.mkdup.cram"),file("${sampleId}-${part}.mkdup.cram.crai"), emit: bams
    path("${sampleId}-${part}.log.bwamem")
    path("${sampleId}-${part}.hla.all") , optional: true
    path("${sampleId}-${part}.log.hla") , optional: true

   script:
    def aln="bwa-mem2"
    //we define the aln tool
    if(params.aligner=="bwa"){
	aln="bwa"
    }
    if(params.debug == true){
    """
    echo "seqtk mergepe $read1 $read2 | ${aln} mem -p -t $task.cpus -R'@RG\\tID:${sampleId}-${part}\\tSM:${sampleId}\\tPL:ill' ${params.ref} - 2> ${sampleId}-${part}.log.bwamem | k8 ${params.alt_js} -p ${sampleId}-${part}.hla ${params.ref}.alt | samtools view -1 - > ${sampleId}-${part}.aln.bam"
    echo "run-HLA ${sampleId}-${part}.hla > ${sampleId}-${part}.hla.top 2> ${sampleId}-${part}.log.hla;"
    echo "touch ${sampleId}-${part}.hla.HLA-dummy.gt; cat ${sampleId}-${part}.hla.HLA*.gt | grep ^GT | cut -f2- > ${sampleId}-${part}.hla.all"
    echo "rm -f ${sampleId}-${part}.hla.HLA*;"
    touch ${sampleId}-${part}.mkdup.cram
    touch ${sampleId}-${part}.mkdup.cram.crai
    touch ${sampleId}-${part}.log.bwamem
    touch ${sampleId}-${part}.hla.all
    """
    }else{
    if(params.hla == true){
    """
	seqtk mergepe $read1 $read2 \\
        | ${aln} mem -p -t $task.cpus -R'@RG\\tID:${sampleId}-${part}\\tSM:${sampleId}\\tPL:ill' ${params.ref} - 2> ${sampleId}-${part}.log.bwamem \\
        | k8 ${params.alt_js} -p ${sampleId}-${part}.hla ${params.ref}.alt \
        | samtools view -Sb -  \
        | samtools fixmate -m - -  \
        | samtools sort -m 1G -@ ${task.cpus - 8 } -  \
        | samtools markdup -O CRAM  --write-index --reference ${params.ref} -@ ${task.cpus - 8} - ${sampleId}-${part}.mkdup.cram

   run-HLA ${sampleId}-${part}.hla > ${sampleId}-${part}.hla.top 2> ${sampleId}-${part}.log.hla;
	 touch ${sampleId}-${part}.hla.HLA-dummy.gt; cat ${sampleId}-${part}.hla.HLA*.gt | grep ^GT | cut -f2- > ${sampleId}-${part}.hla.all;
	 rm -f ${sampleId}-${part}.hla.HLA*;
    """
    }
    else if (params.alt == true){
     """
	seqtk mergepe $read1 $read2  \\
  	| ${aln} mem -p -t $task.cpus  -R'@RG\\tID:${sampleId}-${part}\\tSM:${sampleId}\\tPL:ill' ${params.ref} - 2> ${sampleId}-${part}.log.bwamem \\
  	| k8 ${params.alt_js} -p ${sampleId}-${part}.hla hs38DH.fa.alt \
    | samtools view -Sb -  \
    | samtools fixmate -m - -  \
    | samtools sort -m 1G -@ ${task.cpus - 8 } -  \
    | samtools markdup -O CRAM  --write-index --reference ${params.ref} -@ ${task.cpus - 8} - ${sampleId}-${part}.mkdup.cram

     """
    }else{
	//normal mapping mode
     """
	  seqtk mergepe $read1 $read2 \\
  	| ${aln} mem -p -t $task.cpus  -R'@RG\\tID:${sampleId}-${part}\\tSM:${sampleId}\\tPL:ill' ${params.ref} - 2> ${sampleId}-${part}.log.bwamem \\
    | samtools view -Sb -  \
    | samtools fixmate -m - -  \
    | samtools sort -@ ${task.cpus} -  \
    | samtools markdup -O CRAM  --write-index --reference ${params.ref} -@ ${task.cpus} - ${sampleId}-${part}.mkdup.cram

     """

    }
  }

}

//merge bams by sample
process MERGEB{

  tag "$sampleId-merge"
  publishDir "$params.outdir/CRAM", mode: "copy"
  container "oras://community.wave.seqera.io/library/bwa-mem2_instrain_multiqc_qualimap_samtools:850f96dbd001458f"

  input:
  tuple val(sampleId), val(parts), file(cramFiles), file(craiFiles)
  path reference

  output:
  tuple val(sampleId), file("${sampleId}.merged.cram"), file("${sampleId}.merged.cram.crai") ,emit: mbams

  script:
  def filesb = cramFiles instanceof List ? cramFiles : [cramFiles]
  if ( filesb.size() == 1 ) {
	if ( params.debug == true ) {
                """
                echo ln -s ${cramFiles[0]} ${sampleId}.merged.cram
                touch ${sampleId}.merged.cram.crai
                touch ${sampleId}.merged.cram
                """
            } else {
                """
                ln -s ${cramFiles[0]} ${sampleId}.merged.cram
                ln -s ${craiFiles[0]} ${sampleId}.merged.cram.crai
                """
            }
  }else{
  if(params.debug == true){
  """
    echo samtools merge --write-index --reference  ${reference} -O CRAM -@ $task.cpus -f ${sampleId}.merged.cram ${cramFiles}
    touch ${sampleId}.merged.cram
    touch ${sampleId}.merged.cram.crai
  """
  }else{
  """
  samtools merge --write-index --reference ${reference} -O CRAM -@ $task.cpus -f ${sampleId}.merged.cram ${cramFiles}
  """
  }
 }
}

process QUALIMAP {
    tag "${sampleId}"
    publishDir "${params.outdir}/qualimap", mode: 'copy'


    container "oras://community.wave.seqera.io/library/bwa-mem2_instrain_multiqc_qualimap_samtools:850f96dbd001458f"

    input:
    tuple val(sampleId), path(cram), path(crai)
    path reference
    path fai

    output:
    path "${sampleId}", emit: qualimap_results

    script:
    if(params.debug){
    """
    echo qualimap bamqc \
        -bam ${cram} \
        -outdir qualimap_results/${sampleId} \
        -nt ${task.cpus} \
        --java-mem-size=${task.memory.toGiga()}G
        mkdir ${sampleId}
    """
    }else{
    """
    # reference FASTA (necesario para descomprimir CRAM)
    REF=${reference}

     # ---------- 1)  Crear FIFO ----------
    fifo="${sampleId}.bam.pipe"
     mkfifo "\$fifo"

    # ---------- 2)  Convertir CRAM → BAM en segundo plano ----------
     samtools view -@ ${task.cpus - 2 } -T ${reference} -b ${cram} > "\$fifo" &
     sam_pid=\$!

     # ---------- 3)  Lanzar Qualimap leyendo del FIFO ----------
      qualimap bamqc \
        -bam    "\$fifo" \
        -outdir ${sampleId} \
        -nt     ${task.cpus} \
        --java-mem-size=${task.memory.toGiga()}G
       rc=\$?

     # ---------- 4)  Limpiar ----------
     kill "\$sam_pid" 2>/dev/null || true   # por si qualimap terminó antes
     rm -f "\$fifo"

    exit "\$rc"
    """
    }
}


// Run deepVariant for a single sample
process DEEPVARIANT_AUTOSOMES{

	tag "$sampleId-deepVariant"
	publishDir "$params.outdir/deepvariant_persample", mode: "copy"

    container "docker.io/google/deepvariant:1.8.0"

	input:
	tuple val(sampleId), file(cram), file(crai)
  path reference
  path fai
	output:
  tuple val(sampleId), file("${sampleId}.autosomes.vcf.gz")  , emit: vcf
	tuple val(sampleId), file("${sampleId}.autosomes.g.vcf.gz") , emit: gvcf
  path "${sampleId}.autosomes.g.vcf.gz", emit: gvcf_file

	script:
	if (params.debug == true){
  """
  REGION_AUTOSOMES="${params.autosomes}"
        echo  /opt/deepvariant/bin/run_deepvariant --model_type WGS \\
                --ref=${reference} \\
                --reads=${cram} \\
                --output_vcf=${sampleId}.autosomes.vcf.gz \\
                --output_gvcf=${sampleId}.autosomes.g.vcf.gz \\
                --intermediate_results_dir=tmp \\
                --num_shards=${task.cpus} \\
                --vcf_stats_report=true \\
                --regions \"\$REGION_AUTOSOMES\"

       touch ${sampleId}.autosomes.vcf.gz ${sampleId}.autosomes.g.vcf.gz
	"""
	}
	else {
    	"""
      # we call in autosomes
      REGION_AUTOSOMES="${params.autosomes}"
      /opt/deepvariant/bin/run_deepvariant --model_type WGS  \\
      --ref=${params.ref} \\
      --reads=${cram} \\
      --output_vcf=${sampleId}.autosomes.vcf.gz \\
      --output_gvcf=${sampleId}.autosomes.g.vcf.gz \\
      --intermediate_results_dir=tmp \\
      --num_shards=${task.cpus} \\
      --vcf_stats_report=true \\
      --regions \"\$REGION_AUTOSOMES\"
   	"""
	}
}

// Run GLnexus to merge deepvariants g.vcfs

process GLNEXUS_DEEPVARIANT_AUTOSOMES{

	tag "$sampleId-glx"
	publishDir "$params.outdir/glx", mode: "copy"

  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/glnexus:1.4.1--h40d77a6_0' :
        'biocontainers/glnexus:1.4.1--h40d77a6_0' }"

	input:
  tuple val(sampleId), path(gvcfs)

	output:
   tuple val(sampleId), path("${sampleId}.autosomes.glx.bcf"), emit: bcf

	script:
	def files = gvcfs.join(' ')

  if (params.debug == true){
    """
    echo   glnexus_cli  --threads $task.cpus \\
          --mem-gbytes $task.memory.giga \\
          --config DeepVariant_unfiltered ${files}
    touch ${sampleId}.autosomes.glx.bcf

    """
  }else{
    """
    glnexus_cli  --threads $task.cpus \\
        --mem-gbytes $task.memory.giga \\
        --config DeepVariant_unfiltered ${files} >  ${sampleId}.autosomes.glx.bcf
    """
  }
}


// Run deepVariant for a single sample
process DEEPVARIANT_SEX{

	tag "$sampleId-deepVariant"
	publishDir "$params.outdir/deepvariant_persample", mode: "copy"

    container "docker.io/google/deepvariant:1.8.0"

	input:
	tuple val(sampleId), file(cram), file(crai)
  path reference
  path fai
	output:
  tuple val(sampleId), file("${sampleId}.sex.vcf.gz")  , emit: vcf
	tuple val(sampleId), file("${sampleId}.sex.g.vcf.gz") , emit: gvcf
  path "${sampleId}.sex.g.vcf.gz", emit: gvcf_file

	script:
	if (params.debug == true){
  """
  REGION_SEX="${params.sex}"
  HAPLOID_CONTIGS="chrX,chrY"
  PAR_BED="GRCh38_PAR.bed"
        echo  /opt/deepvariant/bin/run_deepvariant --model_type WGS \\
                --ref=${reference} \\
                --reads=${cram} \\
                --output_vcf=${sampleId}.autosomes.vcf.gz \\
                --output_gvcf=${sampleId}.autosomes.g.vcf.gz \\
                --intermediate_results_dir=tmp \\
                --num_shards=${task.cpus} \\
                --vcf_stats_report=true \\
                --haploid_contigs "\$HAPLOID_CONTIGS" \\
                --par_regions_bed "${baseDir}/aux/\$PAR_BED" \\
                --regions \"\$REGION_SEX\"

       touch ${sampleId}.sex.vcf.gz ${sampleId}.sex.g.vcf.gz
	"""
	}
	else {
    	"""
      #we call in sex chromosomes in hg38
      REGION_SEX="${params.sex}"
      HAPLOID_CONTIGS="chrX,chrY"
      PAR_BED="GRCh38_PAR.bed"
        /opt/deepvariant/bin/run_deepvariant --model_type WGS \\
        --ref=${reference} \\
        --reads=${cram} \\
        --output_vcf=${sampleId}.sex.vcf.gz \\
        --output_gvcf=${sampleId}.sex.g.vcf.gz \\
        --intermediate_results_dir=tmp \\
        --num_shards=${task.cpus} \\
        --vcf_stats_report=true \\
        --haploid_contigs \"\$HAPLOID_CONTIGS\" \\
        --par_regions_bed "${baseDir}/aux/\$PAR_BED" \\
        --regions \"\$REGION_SEX\"
   	"""
	}
}

// Run GLnexus to merge deepvariants g.vcfs

process GLNEXUS_DEEPVARIANT_SEX{

	tag "$sampleId-glx"
	publishDir "$params.outdir/glx", mode: "copy"

  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/glnexus:1.4.1--h40d77a6_0' :
        'biocontainers/glnexus:1.4.1--h40d77a6_0' }"

	input:
  tuple val(sampleId), path(gvcfs)

	output:
   tuple val(sampleId), path("${sampleId}.sex.glx.bcf"), emit: bcf

	script:
	def files = gvcfs.join(' ')

  if (params.debug == true){
    """
    echo   glnexus_cli  --threads $task.cpus \\
          --mem-gbytes $task.memory.giga \\
          --config DeepVariant_unfiltered ${files}
    touch ${sampleId}.sex.glx.bcf

    """
  }else{
    """
    glnexus_cli  --threads $task.cpus \\
        --mem-gbytes $task.memory.giga \\
        --config DeepVariant_unfiltered ${files} >  ${sampleId}.sex.glx.bcf
    """
  }
}





process DEPTH{
    tag "$sampleId-depth"
    //label 'process_medium'

    publishDir "$params.outdir/DEPTH", mode: "copy"

    input:
    tuple val(sampleId), file(cram), file(crai)


    output:
    path("${sampleId}.depth.*"), emit: depth

    script:
    if(params.debug == true){
    """
    echo mosdepth -f ${params.ref} -t $task.cpus ${sampleId}.depth $cram
    touch ${sampleId}.depth.mosdepth.dist.txt
    touch ${sampleId}.depth.mosdepth.summary.txt
    """
    }else{
    """
    mosdepth -f ${params.ref} -t $task.cpus ${sampleId}.depth $cram
    """
    }
}



//we declare the workflow for index

workflow {
    // TODO do a parameter check
    //we read pairs from regex
    if(params.reads != null){
    // --reads "./reads/B087*.merge.{1,2}.fq.gz"
    read_pairs_ch = channel.fromFilePairs(params.reads)
    }else if(params.csv != null){
    //we reads pairs from csv
    read_pairs_ch=Channel.fromPath(params.csv) \
        | splitCsv(header:true) \
        | map { row-> tuple(row.sampleId, row.part,  file(row.read1), file(row.read2)) }
    }else{
        println "Error: reads regex or path"
    }

   // read_pairs_ch.view()
    ref = file(params.ref)
    ref_fai = file(params.ref+".fai")
    //fastqc read quality
    FASTQC(read_pairs_ch)
    //read aligment alt/hla
    BWAMEM(read_pairs_ch)
    //we do merge the bams by sample ID
    groups=BWAMEM.out.bams.groupTuple(by: 0)
    //groups.view()
   // groups.view()
    MERGEB(groups,ref)
    //Quality of alignments
    //MERGEB.out.mbams.view()
    QUALIMAP(MERGEB.out.mbams,ref)
    //Coverage Stats from cram files
    DEPTH(MERGEB.out.mbams)
    //variant calling witn DeepVariant
    DEEPVARIANT_AUTOSOMES(MERGEB.out.mbams,ref,ref_fai)
    allgvcf=DEEPVARIANT_AUTOSOMES.out.gvcf_file.collect().map { gvcf_list ->
        tuple('all', gvcf_list)
    }
    //allgvcf.view()
    GLNEXUS_DEEPVARIANT_AUTOSOMES(allgvcf)

    DEEPVARIANT_SEX(MERGEB.out.mbams,ref,ref_fai)
    allgvcfsex=DEEPVARIANT_SEX.out.gvcf_file.collect().map { gvcf_list ->
        tuple('all', gvcf_list)
    }
    //allgvcf.view()
    GLNEXUS_DEEPVARIANT_SEX(allgvcfsex) 
}
