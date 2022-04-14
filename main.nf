#!/usr/bin/env nextflow


////////////////////////////////////////////////////
/* --           Gather Parameters              -- */
////////////////////////////////////////////////////
/*
 * Set index directory path
 */
index_dir =  "${baseDir}/index"

/*
 * Bowtie
 */
if (params.bowtie.run != false){
    // Check if index is build 
    genome = file("${params.bowtie.genome}")
    genome_name = "${genome.baseName}"
    genome_path = "${index_dir}/bowtie/${genome_name}/${genome_name}"
    params.bowtie.directory = "${index_dir}/bowtie/${genome_name}/${genome_name}"
    genome_test = file("${index_dir}/bowtie/${genome_name}/${genome_name}.1.ebwt")
    genome_exists = "${genome_test.exists()}"

    if (genome_exists == true){
        params.bowtie.build = false
    } else {
        params.bowtie.build = true
    }

} else {
    params.bowtie.build = false
}

/*
 * Tailor
 */
if (params.tailor.run != false) {
    // Check if tailor index is build 
    tailor = file("${params.bowtie.genome}")
    tailor_name = "${genome.baseName}"
    tailor_path = "${index_dir}/tailor/${genome_name}/${genome_name}"
    params.tailor.directory = "${index_dir}/tailor/${genome_name}/${genome_name}"
    tailor_test = file("${index_dir}/tailor/${genome_name}/${genome_name}.t_bwt.bwt")
    tailor_exists = "${tailor_test.exists()}"

    if (tailor_exists == true){
        params.tailor.build = false
    } else {
        params.tailor.build = true
    }

} else {
    params.tailor.build = false
}

/*
 * Features
 */
if (params.counter.normalization == false){
    params.counter.normalizeCommand = ""
} else {
    params.counter.normalizeCommand = "-n ${params.counter.normalization}"
}

/*
 * Bowtie alignment params
 */
if (params.bowtie.mismatch == false){
    params.bowtie.v = "-v 0"
} else {
    params.bowtie.v = "-v ${params.bowtie.mismatch}"
}

if (params.bowtie.multimap == false){
    params.bowtie.m = ""
} else if (params.bowtie.multimap == 0){
    params.bowtie.m = ""
} else {
    params.bowtie.m = "-m ${params.bowtie.multimap}"
}

/*
 * Bin
 */
params.bin = "${baseDir}/bin"

/*
 * Reads 
 */
reads = Channel
        .fromPath( params.reads.directory )
        .ifEmpty { error "No reads in reads directory" }
        .map { file -> tuple(file.simpleName, file) }

////////////////////////////////////////////////////
/* --         PRINT PARAMETER SUMMARY          -- */
////////////////////////////////////////////////////
def summary = [:]
if(workflow.revision)          summary['Pipeline Release'] = workflow.revision

summary['Project Name']        = params.project.name
summary['Project Directory']   = params.project.directory
summary['Results']             = params.results.directory
summary['Index Directory']     = "${baseDir}/index"
summary['Reads']               = params.reads.directory

summary['Bowtie build']        = params.bowtie.build
summary['Tailor build']        = params.tailor.build

summary['Genome']              = params.bowtie.genome
summary['Junctions']           = params.bowtie.juncs
summary['Genome Generate']     = params.bowtie.build
summary['Mismatches']          = params.bowtie.v
summary['Multimapper']         = params.bowtie.m

summary['TrimGalore']          = params.trimgalore.run
summary['Min Length']          = params.trimgalore.minL
summary['Max Length']          = params.trimgalore.maxL

summary['Counter']             = params.counter.run
summary['Bed Ref.']            = params.counter.reference
summary['Features']            = params.counter.features
summary['Normalization']       = params.counter.normalization

summary['DGE']                 = params.dge.run
summary['Comparisons']         = params.dge.comparisons

summary['python version']      = params.python.version

summary['Current home']        = "$HOME"
summary['Current user']        = "$USER"
summary['Current path']        = "$PWD"
summary['Script dir']          = workflow.projectDir

logfile = file("${params.project.directory}/NF-run-${params.project.name}.log")
logfile.write("Nextflow parameter summary")
logfile.write("")
summary.each {  logfile.append("\n${it}\n") }

/*
 * Enable dsl2 syntax and import modules
 * d-s-l = design specific language
 */
nextflow.enable.dsl=2

/*
 * Include modules 
 */
include {
    TRIM_GALORE;
    BOWTIE_INDEX;
    BOWTIE_ALIGN_GENOME;
    COUNT_FEATURES;
    TRANSCRIPTS;
    DGE;
    MASTER_TABLE;
    RBIND_COUNTS;
    TAILOR_INDEX;
    TAILOR_MAP
} from './modules.nf' 

////////////////////////////////////////////////////
/* --               Pipeline Logic             -- */
////////////////////////////////////////////////////
workflow {

    // STEP 1: Trimming...
    if (params.trimgalore.run == true) {
        TRIM_GALORE(params.trimgalore.minL, params.trimgalore.maxL, reads)
        TRIM_GALORE.out.collapsed_fa.tap{ collapsed_fa }
    } else {
        collapsed_fa = Channel.fromPath( params.reads.directory ).ifEmpty { error "No reads in reads directory" }.map { file -> tuple(file.simpleName, file) }
    }

    // STEP 2: Bowtie alignment
    if (params.bowtie.build == true){
        BOWTIE_INDEX(params.bowtie.genome, params.bowtie.juncs)
        BOWTIE_INDEX.out.bowtie_index.tap{ bowtie_index }
    } else {
        bowtie_index = genome_path
    }

    BOWTIE_ALIGN_GENOME(bowtie_index, collapsed_fa)

    // STEP 3: Feature count --optional--
    if (params.counter.run == true){
        COUNT_FEATURES(params.counter.features, params.counter.normalizeCommand, params.counter.reference, BOWTIE_ALIGN_GENOME.out.bowtie_alignment)
    }

    // STEP 4: Transcript alignment --optional--
    if (params.transcripts.run == true){

        // STEP 4.1: IF counter is true one output channel will be normalization constants, 
        // however, if not true normalization constants will be set to the depth output by bowtie 
        // alignment.
        if (params.counter.run == true) {
            transcripts_input = collapsed_fa.join(COUNT_FEATURES.out.normalization_constants)
        } else {
            transcripts_input = collapsed_fa.join(BOWTIE_ALIGN_GENOME.out.normalization_constants)
        }

        TRANSCRIPTS(params.transcripts.file, transcripts_input)
    }

    // STEP 5: Make master table input, if transcripts is true, or if counter is true or if both are true
    if (params.transcripts.run == true && params.counter.run == true) {
        counts = COUNT_FEATURES.out.counts.mix(TRANSCRIPTS.out.transcripts_output).groupTuple()
        RBIND_COUNTS( counts )
        RBIND_COUNTS.out.tables.tap{ master_table_input }
        make_master_table = true
    } else if (params.transcripts.run == true && params.counter.run == false) {
        counts = TRANSCRIPTS.out.transcripts_output
        RBIND_COUNTS( counts )
        RBIND_COUNTS.out.tables.tap{ master_table_input }
        make_master_table = true
    } else if (params.transcripts.run == false && params.counter.run == true) {
        counts = COUNT_FEATURES.out.counts
        RBIND_COUNTS( counts )
        RBIND_COUNTS.out.tables.tap{ master_table_input }
        make_master_table = true
    } else {
        make_master_table = false
    }

    // STEP 6: make master table
    if (make_master_table == true) {
        MASTER_TABLE( params.project.name, master_table_input.collect() )
    }

    // STEP 7. differential gene expression analysis
    if (params.dge.run == true) {
        DGE( params.dge.comparisons, MASTER_TABLE.out.tables.flatten() )
    }

    // STEP 8. Tailor --optional--
    if (params.tailor.run == true) {
        
        if (params.tailor.build == true){
            TAILOR_INDEX(params.bowtie.genome)
            TAILOR_INDEX.out.tailor_index.tap{ tailor_index }
        } else {
            tailor_index = tailor_path
        }

        tailor_input = BOWTIE_ALIGN_GENOME.out.tailor_input.join(COUNT_FEATURES.out.normalization_constants)    
        TAILOR_MAP(params.bowtie.genome, tailor_index, params.counter.features, params.counter.reference, tailor_input)
    }

}