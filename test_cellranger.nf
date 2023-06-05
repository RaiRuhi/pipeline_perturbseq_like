// Declare syntax version
nextflow.enable.dsl=2
// Script parameters
//params.GTF_GZ_LINK = 'http://ftp.ensembl.org/pub/release-106/gtf/homo_sapiens/Homo_sapiens.GRCh38.106.gtf.gz'
//params.CELLRANGER_REF = '/data/gersbachlab/Ruhi/pipeline_perturbseq_like/refdata-gex-GRCh38-2020-A'
//params.TRANSCRIPTOME_REFERENCE = "human"
//params.KALLISTO_BIN = '/data/reddylab/Ruhi/miniconda3/envs/perturbseq_like_pipeline/bin/kallisto'
//params.GENOME = 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz'
//params.GUIDE_FEATURES = '/data/gersbachlab/Ruhi/Nextflow/pipeline_perturbseq_like/df_from_gasperini_tss.xlsx'
//params.CHEMISTRY = '10xv2'
//params.THREADS = 15
//params.DISTANCE_NEIGHBORS = 1000000
//params.IN_TRANS = "FALSE"
//params.FASTQ_FILES_TRANSCRIPTS = ['/data/reddylab/Alex/IGVF/Jamborees/230403/data/bam_pilot_scrna_1/K1000_CRISPRi_cells_r1_SI-GA-G1_MissingLibrary_1_HMMKLBGX3/bamtofastq_S1_L001_R1_001.fastq.gz /data/reddylab/Alex/IGVF/Jamborees/230403/data/bam_pilot_scrna_1/K1000_CRISPRi_cells_r1_SI-GA-G1_MissingLibrary_1_HMMKLBGX3/bamtofastq_S1_L001_R2_001.fastq.gz ']
//params.FASTQ_NAMES_TRANSCRIPTS = ['S1_L1']
//params.FASTQ_FILES_GUIDES = ['/data/reddylab/Alex/IGVF/Jamborees/230403/data/bam_pilot_guide_1/K1000_CRISPRi_gRNA_r1_CGTTACCG_MissingLibrary_1_HMMKLBGX3/bamtofastq_S1_L001_R1_001.fastq.gz /data/reddylab/Alex/IGVF/Jamborees/230403/data/bam_pilot_guide_1/K1000_CRISPRi_gRNA_r1_CGTTACCG_MissingLibrary_1_HMMKLBGX3/bamtofastq_S1_L001_R2_001.fastq.gz ']
//params.FASTQ_NAMES_GUIDES = ['S1_L1']
//params.WHITELIST= '/data/gersbachlab/Ruhi/Nextflow/737K-august-2016.txt'
//params.ADDGENENAMES= '/data/gersbachlab/Ruhi/pipeline_perturbseq_like/test_genename.txt'
//params.CREATE_REF = false
//params.DIRECTION  = "both"





 

workflow {
   dir_images_composition_scrna =  compositionREADSscRNA(Channel.from(params.FASTQ_NAMES_TRANSCRIPTS), Channel.from(params.FASTQ_FILES_TRANSCRIPTS))
   dir_images_composition_guides = compositionREADSGuides(Channel.from(params.FASTQ_NAMES_GUIDES), Channel.from(params.FASTQ_FILES_GUIDES)) 
   gtf_out = downloadGTF(params.GTF_GZ_LINK)
   downloadReference(Channel.of(params.TRANSCRIPTOME_REFERENCE), Channel.of(params.KALLISTO_BIN) )
   downloadGenome (Channel.of(params.GENOME))
   guide_feature_preprocessed = guidePreprocessing(params.GUIDE_FEATURES)
   creatingGuideRef( downloadGenome.out.genome, Channel.of(params.KALLISTO_BIN), guide_feature_preprocessed.guide_features, params.CREATE_REF )
    
    
    cellranger_out = cellranger (
                 Channel.from(params.FASTQ_NAMES_TRANSCRIPTS),
                 Channel.from(params.FASTQ_FILES_TRANSCRIPTS),
                 Channel.from(params.FASTQ_NAMES_GUIDES),
                 Channel.from(params.FASTQ_FILES_GUIDES),
                 Channel.from(params.GUIDE_FEATURES),
                 Channel.from(params.CELLRANGER_REF),
                 Channel.of(params.CHEMISTRY).collect(),
                 Channel.of(params.THREADS).collect(),
                 Channel.of(params.WHITELIST).collect(),
                )
           
    
     //dir_count_files = cellranger_guide_out.ks_guide_out.join(cellranger_rna_out.ks_transcripts_out, remainder: true).flatten().toList().view()
    dir_count_files = cellranger_out.mtx_file.flatten().toList().view()
    //dir_count_files.view()
    
    df_initialized = preprocessing(dir_count_files)
    

    
    //out_processed = filtering(df_initialized.df_initial_files, dir_count_files )
    
    
    concact_prefiltering_out = concact_prefiltering(df_initialized.df_initial_files,
                                                    dir_count_files,
                                                params.EXPECTED_CELL_NUMBER,
                                                params.MITO_SPECIE,
                                                params.MITO_EXPECTED_PERCENTAGE,
                                                params.PERCENTAGE_OF_CELLS_INCLUDING_TRANSCRIPTS,
                                                params.TRANSCRIPTS_UMI_TRHESHOLD)
                                                
    //concact_prefiltering_out.guide_ann
    //concact_prefiltering_out.transcripts_ann
    
    
    moun_raw_creation = moun_raw_creation(  concact_prefiltering_out.guide_ann,
                                            concact_prefiltering_out.transcripts_ann,
                                            downloadReference.out.t2t_transcriptome_index,     
                                            gtf_out.gtf )
                                            
     //moun_raw_creation = moun_raw_creation( cellranger_out.mtx_file )                                        
    
    merge_bin_and_muon_out =  merge_bin_and_muon( moun_raw_creation.raw_moun_data, params.GUIDE_UMI_LIMIT )    
    
    pert_loader_out = PerturbLoaderGeneration(merge_bin_and_muon_out.muon_processed, gtf_out.gtf , params.DISTANCE_NEIGHBORS, params.IN_TRANS, params.ADDGENENAMES )
    runSceptre_out = runSceptre(pert_loader_out.perturb_piclke, params.DIRECTION)
    create_anndata_from_sceptre_out = create_anndata_from_sceptre(runSceptre_out.sceptre_out_dir,df_initialized)
    
    
    
}


process downloadGTF {
    input:
    val gtf_gz_path
    output:
    path "transcripts.gtf" , emit: gtf
    script:
    """    
    wget -O - $gtf_gz_path | gunzip -c > transcripts.gtf
    """
}

process downloadReference {
    input:
    val ref_name 
    path k_bin 
    
    output:
    path "transcriptome_index.idx" , emit: transcriptome_idx
    path "transcriptome_t2g.txt"   , emit: t2t_transcriptome_index

    script:
    """
        kb ref -d $ref_name -i transcriptome_index.idx -g transcriptome_t2g.txt --kallisto ${k_bin}
    """
}

process downloadGenome {
    input:
    val genome_path
    output:
    path "genome.fa.gz" , emit: genome
    script:
    """
        wget $genome_path -O genome.fa.gz 
    """


}

process guidePreprocessing {
    cache 'lenient'
    debug true
    input:
    path (guide_input_table)
    output:
    path "guide_features.txt" , emit: guide_features
    script:
    """    
    guide_table_processing.py  $guide_input_table
    """
}





process creatingGuideRef {
    cache 'lenient'
    input:
    val genome_path
    path k_bin
    path guide_features
    val create_ref
    output:
    path "guide_index.idx" ,  emit: guide_index
    path "t2guide.txt" , emit: t2tguide_index
    script:
    
    """
        kb ref -i guide_index.idx -f1 $genome_path -g t2guide.txt --kallisto $k_bin  --workflow kite $guide_features 

    """


    
}

process compositionREADSscRNA {
    debug true
    input:
    tuple val(out_name_dir)
    tuple val(string_fastqz)
    output:
    path ("transcript_${out_name_dir}_composition"),  emit: composition_plot_dir

    script:
    
        """
        mkdir transcript_"${out_name_dir}_composition"
        fq_composition.py $string_fastqz transcript_${out_name_dir}                                                                
        """
} 


process compositionREADSGuides {
    debug true
    input:
    tuple val(out_name_dir)
    tuple val(string_fastqz)
    output:
    path ("guide_${out_name_dir}_composition"),  emit: composition_plot_dir

    script:
    
        """
        mkdir "guide_${out_name_dir}_composition"
        fq_composition.py $string_fastqz guide_${out_name_dir}                                                                
        """
} 

process cellranger {
        input:
        tuple val(rna_sample_name)
        tuple val(rna_string_fastqz)
        tuple val(guide_sample_name)
        tuple val(guide_string_fastqz)
        tuple val(guide_features)
        tuple val(cellranger_ref)
        tuple val(chemistry)
        tuple val(threads)
        tuple val(whitelist)
        output:
        path ('outs/filtered_feature_bc_matrix'),  emit: mtx_file

        script:
    
        """
        
        chmod 700 /data/gersbachlab/Ruhi/pipeline_perturbseq_like/generating_cell_range_inputs.py
        /data/gersbachlab/Ruhi/pipeline_perturbseq_like/generating_cell_range_inputs.py --rna_sample_name $rna_sample_name --rna_string_fastqz "$rna_string_fastqz" --guide_sample_name $guide_sample_name --guide_string_fastqz "$guide_string_fastqz" --guide_features $guide_features --cellranger_ref $cellranger_ref --chemistry $chemistry --threads 15 --whitelist $whitelist
        #The script generating_cell_range_inputs.py should create the library.csv and feature_ref.csv 
        #remove the comment to run cellranger
        /data/gersbachlab/Ruhi/pipeline_perturbseq_like/cellranger-7.1.0/bin/cellranger count --id=gasperini_01 --libraries=/data/gersbachlab/Ruhi/pipeline_perturbseq_like/library.csv --transcriptome=$cellranger_ref --feature-ref=/data/gersbachlab/Ruhi/pipeline_perturbseq_like/feature_ref.csv 
        
        #Artificially creating the output. Remove it when cellranger starts running.
        
        mkdir  outs
        mkdir outs/filtered_feature_bc_matrix/
        mkdir 'outs/filtered_feature_bc_matrix/matrix.mtx.gz'
        
      
        
        """
} 

process capture_variables_and_save_list{
    debug true
    input:
    val received
    val out_name
    output:
    path "${out_name}.txt",  emit: out_file
    script:
    """
    echo  '${received}' > ${out_name}.txt 
    """
}


process preprocessing {
    debug true
    input:
    path (count_list)
    output:
    path 'initial_preprocessing_file_names.txt', emit: df_initial_files
    script:
    """    
    preprocessing.py  ${count_list} 

    """   
    
}


process filtering{
    debug true
    input:
    path (path_df)
    path (all_files_context)
    output:
    path 'results_per_lane/processed_anndata_guides_data.h5ad', emit: guide_ann
    path 'results_per_lane/processed_anndata_transcripts_data.h5ad',  emit: transcripts_ann
    script:
    """
    #use -merge to merge the guides
    # I need to add these parameters to the pipeline config  mito...cellnumber...merge...guide_limit
    filtering_and_lane_merging.py --path ${path_df} --expected_cell_number 8000 --mito_specie hsapiens --mito_expected_percentage 0.2 --percentage_of_cells_to_include_transcript 0.2  --guide_umi_limit 5

    """
}


process concact_prefiltering{
    debug true
    input:
    path (path_df)
    path (dir_count_context)
    val  (expected_cell_number)
    val  (mito_specie)
    val  (mito_expected_percentage)
    val (percentage_of_cells_to_include_transcript)
    val  (transcripts_umi_treshold)
    output:
    path 'results_per_lane/full_raw_guide_ann_data.h5ad', emit: guide_ann
    path 'results_per_lane/full_raw_scrna_ann_data.h5ad',  emit: transcripts_ann
    
    
    
    
    script:
    """
    #chmod 700 /n/scratch3/users/l/lf114/pipeline_perturbseq_like/bin/concact_and_pre_filtering.py; 
    concact_and_pre_filtering.py --path ${path_df} --expected_cell_number $expected_cell_number --mito_specie $mito_specie --mito_expected_percentage $mito_expected_percentage --percentage_of_cells_to_include_transcript  $percentage_of_cells_to_include_transcript --transcripts_umi_treshold $transcripts_umi_treshold 
    """
}

process moun_raw_creation{
    debug true
    input:
    path (ann_guide)
    path (ann_exp)
    path (transcript_file)
    path (gtf_in)
    output:
    path 'raw_mudata_guide_and_transcripts.h5mu', emit: raw_moun_data
      
    script:
    """
    #chmod 700 /n/scratch3/users/l/lf114/pipeline_perturbseq_like/bin/muon_creation.py; 
    muon_creation.py --ann_guide $ann_guide --ann_exp $ann_exp --gtf_in $gtf_in --transcript_file $transcript_file
    
    """
}


process merge_bin_and_muon {
    debug true
    input:
    path (muon_data)
    val  (guide_umi_limit)
    output:
    path 'processed_mudata_guide_and_transcripts.h5mu', emit: muon_processed
    script:
    """
    echo "merge option will  work in  future versions";
    merge_bin_and_muon.py --muon_data  $muon_data  --guide_umi_limit $guide_umi_limit
    """
}

process PerturbLoaderGeneration {
    debug true
    input:
    path (muon_data)
    path (gtf_in)
    val (distance_from_guide)
    val (in_trans)
    val (addgenes)
    output:
    path 'perturbdata.pkl', emit: perturb_piclke
    
   """ 
   PerturbLoader_generation.py --muon_data $muon_data --gtf_in $gtf_in  --distance_from_guide $distance_from_guide --in_trans $in_trans --add_gene_names $addgenes 

    """   
}



process runSceptre {
    debug true
    input:
    path (perturbloader_pickle)
    val (direction)
    output:
    path 'sceptre_out' , emit: sceptre_out_dir
    
   """ 
   mkdir sceptre_out
   mv $perturbloader_pickle sceptre_out
   cd sceptre_out
   runSceptre.py $perturbloader_pickle  $direction
   """   
}

process create_anndata_from_sceptre {
    debug true
    input:
    path (sceptre_results_dir)
    path (mudata_processed)
    output:
    path 'mudata_results.h5mu' , emit:  mudata_perturbation_results
    
   """ 
   sceptre_anndata_creation.py $sceptre_results_dir $mudata_processed
   """   
}

