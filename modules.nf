

/*
 * Process 1. 
 */
process TRIM_GALORE {

    label 'low'

    publishDir "$params.results.directory/trim_galore/fastq", mode : 'copy', pattern : "*_trimmed.fq.gz"
    publishDir "$params.results.directory/trim_galore/logs", mode : 'copy', pattern : "*_trimming_report.txt"
    publishDir "$params.results.directory/trim_galore/collapsed", mode : 'copy', pattern : "*.trimmed.uniq.fa"

    input : 
        val min_length
        val max_length
        tuple val(sampleID), val(fastq)
    
    output : 
        path("${sampleID}*_trimming_report.txt")
        tuple val(sampleID), path("${sampleID}_trimmed.fq.gz")
        tuple val(sampleID), path("${sampleID}.trimmed.uniq.fa"), emit : collapsed_fa

    script : 
    """
    #!/bin/bash 

    source activate smrnaseq
    
    trim_galore -j ${task.cpus} -q 30 -e 0.1 --gzip --length ${min_length} --max_length ${max_length} --fastqc --basename ${sampleID} ${fastq}

    zcat ${sampleID}_trimmed.fq.gz | awk '(NR%4==2)' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\\t"; print ">" \$2 ":" \$1 "\\n" \$2 }' > ${sampleID}.trimmed.uniq.fa
    
    """
}

/*
 * Process 2. 
 */
process BOWTIE_INDEX {

    label 'low'
    publishDir "$params.bowtie.directory", mode : 'copy'
    
    input :
        val genome
        val juncs

    output :
        path("*.ebwt")
        path("*chrom_sizes")
        path("${name}*"), emit : bowtie_index

    script : 
    f = file("${genome}")
    name = "${f.baseName}"
    """
    #!/bin/bash

    source activate smrnaseq

    # bowtie
    bowtie-build --quiet ${genome} ${name} --threads ${task.cpus}
    bowtie-build --quiet ${juncs} ${name}_juncs --threads ${task.cpus}

    samtools faidx ${genome}
    cat ${genome}.fai | cut -f1,2 > ${name}_chrom_sizes

    """

}

/*
 * Process 3. Bowtie align genome 
 */
process BOWTIE_ALIGN_GENOME {

    label 'low'

    publishDir "$params.results.directory/alignment", mode : 'copy', pattern : "*.ntm"
    publishDir "$params.results.directory/alignment", mode : 'copy', pattern : "*.bai"
    publishDir "$params.results.directory/alignment", mode : 'copy', pattern : "*.bam"
    publishDir "$params.results.directory/logs", mode : 'copy', pattern : "*.log"

    input : 
        val idx
        tuple val(sampleID), path(fasta)

    output : 
        tuple val(sampleID), path("${sampleID}.unmapped.v0.uni.fq"), emit : tailor_input
        tuple val(sampleID), path("${sampleID}.Aligned.out.ntm"), emit : bowtie_alignment
        path("${sampleID}.depth.tsv"), emit : normalization_constants
        path("${sampleID}*.alignment.log")
        path("*.bam")
        path("*.bai")

    script : 
    """
    #!/bin/bash

    source activate smrnaseq

    #############################################################################
    # Do Alignment with -a --best --strata, get all aligned reads in best stratum
    #############################################################################
    bowtie \\
        -x ${idx} \\
        -f ${fasta} \\
        -p ${task.cpus} \\
        -a \\
        --un ${sampleID}.tmp \\
        --best \\
        --strata ${params.bowtie.v} ${params.bowtie.m} \\
        -S > Aligned.out.genome.sam 2> genome.log
    
    bowtie \\
        -x ${idx}_juncs \\
        -f ${sampleID}.tmp \\
        -p ${task.cpus} \\
        -a \\
        --un ${sampleID}.unmapped.uni.fa \\
        --best \\
        --strata ${params.bowtie.v} ${params.bowtie.m} \\
        -S > Aligned.out.junc.sam 2> junc.log

    ##############################################################################
    # Do Alignment with -v 0, get reads that do not map perfectly, send to Tailor
    ##############################################################################
    bowtie \\
        -x ${idx} \\
        -f ${fasta} \\
        -p ${task.cpus} \\
        -v 0 \\
        -a \\
        --un unmapped.genome.v0.tmp \\
        --best \\
        --strata \\
        -S > mapped.v0.sam

    bowtie \\
        -x ${idx}_juncs \\
        -f unmapped.genome.v0.tmp \\
        -p ${task.cpus} \\
        -v 0 \\
        -a \\
        --un unmapped.genome.junc.v0.tmp \\
        --best \\
        --strata \\
        -S > mapped.junc.v0.sam
    
    ${params.python.version} ${params.bin}/uniq_fasta_to_uniq_fastq.py unmapped.genome.junc.v0.tmp > ${sampleID}.unmapped.v0.uni.fq

    ########################################################################
    # Process alignment
    ########################################################################
    samtools sort -m 1G -@ ${task.cpus} -o ${sampleID}.Aligned.out.genome.sorted.bam Aligned.out.genome.sam

    samtools index -@ ${task.cpus} ${sampleID}.Aligned.out.genome.sorted.bam
    
    bam2bed < ${sampleID}.Aligned.out.genome.sorted.bam | awk -F'\\t' -v OFS='\\t' '{split(\$4,a,":"); split(\$16,b,":"); print \$1,\$2,\$3,a[1],a[2],\$6,b[3] }' > Aligned.out.genome.sorted.bed

    sam2bed < Aligned.out.junc.sam  | awk -F'\\t' -v OFS='\\t' '{split(\$4,a,":"); split(\$16,b,":"); print \$1,\$2,\$3,a[1],a[2],\$6,b[3] }' > tmp

    ${params.python.version} ${params.bin}/junc_bed2bed.py -i tmp -o Aligned.out.junc.bed

    cat Aligned.out.genome.sorted.bed Aligned.out.junc.bed > ${sampleID}.Aligned.out.bed

    ${params.python.version} ${params.bin}/bed_to_ntm.py ${sampleID}.Aligned.out.bed ${sampleID}.Aligned.out.ntm

    depth=\$(${params.python.version} ${params.bin}/add_columns.py ${sampleID}.Aligned.out.ntm 4)

    total_reads=\$(cat ${fasta} | grep "^>" | sed -e 's/.*://g' | ${params.bin}/addCols stdin)
    
    percent_aligned=\$(${params.python.version} -c "print(\$depth/\$total_reads)")

    echo \$depth > depth

    ## Alignment summary
    echo -e "depth\ttotal_reads\tpercent_aligned" > ${sampleID}.alignment.log
    echo -e \$depth'\t'\$total_reads'\t'\$percent_aligned >> ${sampleID}.alignment.log

    ## Normalization constants 
    echo -e 'total\t\$depth' > ${sampleID}.depth.tsv

    """

}

/*
 * Process 4
 */
process COUNT_FEATURES {

    label 'low'

    publishDir "$params.results.directory/bed", mode : 'copy', pattern : "*.bed.tsv"
    publishDir "$params.results.directory/logs", mode : 'copy', pattern : "*.log"
    publishDir "$params.results.directory/normalization", mode : 'copy', pattern : "*-constants.tsv"

    input :
        val features
        val normalize_to
        val reference_annotation
        tuple val(sampleID), val(alignment)
    
    output : 
        tuple val(sampleID), path("${sampleID}*normalization-constants.tsv"), emit : normalization_constants
        tuple val(sampleID), path("*.counts.tsv"), emit : counts
        path("*.counts.tsv"), emit : master_table_input
        path("*.bed.tsv"), emit : bed_counts
        path("*.log")
       
    script : 
    """
    #!/bin/bash

    source activate smrnaseq

    time ${params.python.version} ${params.bin}/count.py \\
        -f ${features} \\
        -a ${reference_annotation} \\
        -i ${alignment} \\
        -o ${sampleID} \\
        ${normalize_to}
        
    """

}

/*
 * Process 5. 
 */
process TRANSCRIPTS {

    label 'low' 

    publishDir "$params.results.directory/transcripts", mode : 'copy', pattern : "*.tsv"

    input : 
        val transcripts
        tuple val(sampleID), val(fasta), val(normalization)
    
    output : 
        tuple val(sampleID), path("${sampleID}.transcripts.counts.tsv"), emit : transcripts_output
        path("*.counts.tsv"), emit : master_table_input

    script : 
    """
    #!/bin/bash

    source activate smrnaseq
    ${params.python.version} ${params.bin}/align_transcripts.py -f ${fasta} \\
        -t ${params.transcripts.file} \\
        -o ${sampleID}.transcripts.counts.tsv \\
        -n ${normalization}
    """
}

/*
 * Process 6. 
 */
process RBIND_COUNTS {

    label 'local'

    publishDir "$params.results.directory/counts", mode : 'copy', pattern : "*.tsv"

    input : 
        tuple val(sampleID), val(files)
    
    output : 
        path("*.tsv"), emit : tables

    script : 
    """
    #!/bin/bash

    source activate smrnaseq
    
    ${params.python.version} ${params.bin}/rbind_tables.py -f "${files}" -o ${sampleID}.counts.tsv

    """
}

/*
 * Process 7. 
 */
process MASTER_TABLE {

    label 'local'

    publishDir "$params.results.directory/master_tables", mode : 'copy', pattern : "*.tsv"

    input :
        val project_name
        val counts

    output : 
        path("*.tsv"), emit : tables
    
    script :
    """
    #!/bin/bash

    source activate smrnaseq

    time ${params.python.version} ${params.bin}/make_master.py "${counts}" ${project_name}

    """
}

/*
 * Process 8. 
 */
process DGE {

    label 'local'

    publishDir "$params.results.directory/DGE", mode : 'copy', pattern : "*.tsv"

    input : 
        val comparisons
        val counts
    
    output : 
        path("*.tsv")
    
    script : 
    """
    #!/bin/bash

    source activate smrnaseq

    time ${params.python.version} ${params.bin}/compare/compare_samples.py -c ${comparisons} -f ${counts}

    """
}

/*
 * Process 9. 
 */
process TAILOR_INDEX {

    label 'low'
    publishDir "$params.tailor.directory", mode : 'copy'

    input :
        val fasta
    
    output :
        path("*")
        path("${name}*"), emit : tailor_index

    script : 
    f = file("${fasta}")
    name = "${f.baseName}"
    """
    #!/bin/bash

    source activate smrnaseq

    ${params.bin}/tailor_v11 build -i ${fasta} -p ${name}
    """

}

/*
 * Process 10. 
 */
process TAILOR_MAP {

    label 'low'

    publishDir "$params.results.directory/tailor/alignment", mode : 'copy', pattern : "*.tailed.bed"
    publishDir "$params.results.directory/tailor/counts", mode : 'copy', pattern : "*.tsv"

    input :
        val genome
        val talor_index
        val features
        val reference_annotation
        tuple val(sampleID), val(fastq), val(normalization_constants)

    output :
        path("${sampleID}.Aligned.out.tailed.bed")
        path("*.tsv")        

    script :
    """
    #!/bin/bash

    source activate smrnaseq

    # Bed file format
    # 1-3 -> chrom, start, end
    # 4 -> seq
    # 5 -> ntm
    # 6 -> strand (all should be +)
    # 7 -> tail
    # 8 -> tail length
    # 9 -> number of locations mapped

    # MUST map sense and be unique mapper
    
    ${params.bin}/tailor_v11 map \\
        -i ${fastq} \\
        -p ${talor_index} \\
        -n ${task.cpus} \\
        2> tailor.log | \\
    tee aligned.sam | \\
    ${params.bin}/tailor_sam_to_bed | \\
    awk -v num=\$nTag -F'\\t' -v OFS="\\t" '{
        if (\$8!="*")
        {
            split(\$4,a,":")

            print \$1,\$2,\$3,a[1],a[2]/\$5,\$6,\$8,\$5

        }
    }' > ${sampleID}.Aligned.out.tailed.bed

    # Count reads
    ${params.python.version} ${params.bin}/tailor_count.py \\
        -i ${sampleID}.Aligned.out.tailed.bed \\
        -a ${reference_annotation} \\
        -f ${features} \\
        -n ${normalization_constants} \\
        -o ${sampleID}-tailed \\
        -g ${genome}

    """
}