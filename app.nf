nextflow.enable.dsl = 2

params.outdir = "./genome_annotation_results"

params.genome1 = ""
params.genome2 = ""
params.config = ""
params.mcool_intragenomic = ""
params.mcool_intergenomic = ""
params.genes1 = ""
params.genes2 = ""

process Repeater2Genome1 {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path genome

    output:
    path "repeats1.bed"

    script:
    """
    source /opt/conda/etc/profile.d/conda.sh
    conda activate nextflow_env

    java -jar -Xms16g -Xmx64g /tools/Repeater2.jar \\
        ${genome} kmer=20 sln=200 -seqshow

    # Convert all GFFs to BED: chrom, start-1, end, repeat_class
    cat *.gff | \\
        awk -F'\\t' '/^[^#]/ && \$4 ~ /^[0-9]+\$/ {
            split(\$1, a, " ")
            print a[1], \$4-1, \$5, \$3
        }' OFS='\\t' > repeats1.bed
    """
}

// Same for genome2
process Repeater2Genome2 {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path genome

    output:
    path "repeats2.bed"

    script:
    """
    source /opt/conda/etc/profile.d/conda.sh
    conda activate nextflow_env

    java -jar -Xms16g -Xmx64g /tools/Repeater2.jar \\
        ${genome} kmer=20 sln=200 -seqshow

    cat *.gff | \\
        awk -F'\\t' '/^[^#]/ && \$4 ~ /^[0-9]+\$/ {
            split(\$1, a, " ")
            print a[1], \$4-1, \$5, \$3
        }' OFS='\\t' > repeats2.bed
    """
}

process ExtractGenes1 {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path extract_genes_script
    path genes

    output:
    path "genes1.bed"

    script:
    """
    source /opt/conda/etc/profile.d/conda.sh
    conda activate syri

    python ${extract_genes_script} --genes ${genes}
    cp genes.bed genes1.bed
    """
}

process ExtractGenes2 {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path extract_genes_script
    path genes

    output:
    path "genes2.bed"

    script:
    """
    source /opt/conda/etc/profile.d/conda.sh
    conda activate syri

    python ${extract_genes_script} --genes ${genes}
    cp genes.bed genes2.bed
    """
}

process CreateSyntheticGenomes {
    publishDir "${params.outdir}/synthetic_genomes", mode: 'copy'

    input:
    path genome1
    path genome2
    path config

    output:
    path "new_genome1.fa"
    path "new_genome2.fa"
    path "chromosome_mapping.tsv"

    script:
    """
    source /opt/conda/etc/profile.d/conda.sh
    conda activate syri

    > new_genome1.fa
    > new_genome2.fa
    > chromosome_mapping.tsv

    echo "synthetic_chrom\tgenome1_original_chroms\tgenome2_original_chroms\tgenome1_lengths\tgenome2_lengths" > chromosome_mapping.tsv
    chrom_num=1

    while IFS= read -r line; do
        # Skip empty lines and comments
        [[ -z "\$line" ]] && continue
        [[ "\$line" =~ ^# ]] && continue

        ref_chromosomes="\${line%:*}"
        qry_chromosomes="\${line#*:}"

        # Build synthetic chromosome for genome1
        echo ">chrom\${chrom_num}" >> new_genome1.fa
        genome1_lengths=""
        genome1_total=0
        for chr in \$(echo "\$ref_chromosomes" | tr ',' ' '); do
            awk -v target="\$chr" '\$1 == ">"target {flag=1; next} /^>/{flag=0} flag{print}' ${genome1} >> new_genome1.fa
            chr_length=\$(seqkit fx2tab -n -l ${genome1} | awk -v target="\$chr" '\$1 == target {print \$NF}')
            genome1_lengths="\$genome1_lengths,\$chr_length"
            genome1_total=\$((genome1_total + chr_length))
        done
        genome1_lengths="\${genome1_lengths:1}"

        # Build synthetic chromosome for genome2
        echo ">chrom\${chrom_num}" >> new_genome2.fa
        genome2_lengths=""
        genome2_total=0
        for chr in \$(echo "\$qry_chromosomes" | tr ',' ' '); do
            awk -v target="\$chr" '\$1 == ">"target {flag=1; next} /^>/{flag=0} flag{print}' ${genome2} >> new_genome2.fa
            chr_length=\$(seqkit fx2tab -n -l ${genome2} | awk -v target="\$chr" '\$1 == target {print \$NF}')
            genome2_lengths="\$genome2_lengths,\$chr_length"
            genome2_total=\$((genome2_total + chr_length))
        done
        genome2_lengths="\${genome2_lengths:1}"

        # Store mapping for later coordinate back-conversion
        echo "chrom\$chrom_num\t\$ref_chromosomes\t\$qry_chromosomes\t\$genome1_lengths\t\$genome2_lengths" >> chromosome_mapping.tsv

        ((chrom_num++))
    done < ${config}
    """
}

process SyRI {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path genome1
    path genome2

    output:
    path "syri.out"

    script:
    """
    source /opt/conda/etc/profile.d/conda.sh
    conda activate syri

    mkdir -p alignments
    minimap2 -ax asm5 --eqx -t 8 \\
      ${genome1} ${genome2} > alignments/aln.sam
    samtools view -b -h -@ 8 alignments/aln.sam > alignments/aln.bam
    samtools sort -@ 8 -o alignments/aln.sorted.bam alignments/aln.bam
    syri -c alignments/aln.sorted.bam -r ${genome1} -q ${genome2} -k -F B --nc 8
    """
}

process MappingSyriResults {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path mapping_syri_results_script
    path syri
    path chromosome_mapping
    path genome1
    path genome2

    output:
    path "mapped_svs.tsv"

    script:
    """
    source /opt/conda/etc/profile.d/conda.sh
    conda activate syri

    python3 ${mapping_syri_results_script} ${syri} ${chromosome_mapping} ${genome1} ${genome1}
    """
}

process EagleC {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path hic

    output:
    path "eaglec_results/output.SV_calls.txt"

    script:
    """
    # Get resolutions and convert to comma-separated list
    RES_LIST=\$(cooler ls ${hic} | awk -F'/' '{printf "%s,", \$NF}' | sed 's/,\$//')
    echo "Resolutions list: \$RES_LIST"

    # Get first resolution for chromosome filtering
    FIRST_RES=\$(echo "\$RES_LIST" | cut -d',' -f1)
    echo "First Resolution: \$FIRST_RES"

    # Get chromosomes as space-separated list
    CHROMOSOMES=\$(cooler dump -t chroms "${hic}::/resolutions/\$FIRST_RES" | awk '\$2 >= 2000000 {printf "%s ", \$1}' | sed 's/ \$//')

    echo "Chromosomes: \$CHROMOSOMES"

    # Calculate intra and inter extend sizes based on resolutions
    # Higher resolution = more confidence = larger extend size
    INTRA_EX=""
    INTER_EX=""

    # Convert RES_LIST to array for processing
    RES_ARRAY=(\$(echo "\$RES_LIST" | tr ',' ' '))
    for res in "\${RES_ARRAY[@]}"; do
        if [ \$res -le 25000 ]; then
            INTRA_EX+="3,"
            INTER_EX+="2,"
        elif [ \$res -le 100000 ]; then
            INTRA_EX+="2,"
            INTER_EX+="1,"
        else
            INTRA_EX+="1,"
            INTER_EX+="1,"
        fi
    done

    INTRA_EX=\${INTRA_EX%,}
    INTER_EX=\${INTER_EX%,}

    echo "Intra extend sizes: \$INTRA_EX"
    echo "Inter extend sizes: \$INTER_EX"

    mkdir -p eaglec_results
    predictSV --mcool ${hic} \\
      --resolutions "\$RES_LIST" \\
      --intra-extend-size "\$INTRA_EX" \\
      --inter-extend-size "\$INTER_EX" \\
      --model-path /tools/EagleC2-models \\
      -g other \\
      -C \$CHROMOSOMES \\
      -p 16 --prob-cutoff-1 0.05 \\
      -O eaglec_results/output
    """
}

process RunComparison {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path compare_script
    path syri
    path eaglec2

    output:
    path "syri_matched.tsv"
    path "matches.tsv"
    path "comparison_stats.txt"

    script:
    """
    source /opt/conda/etc/profile.d/conda.sh
    conda activate python_sv

    python ${compare_script} \
        --syri ${syri} \
        --eaglec2 ${eaglec2} > comparison_stats.txt
    """
}

process AnnotateBreakpointsWithRepeats {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path annotate_breakpoints_script
    path syri_mapped
    path repeats1
    path repeats2
    path genome1
    path genome2
    path genes1
    path genes2

    output:
    path "breakpoints_with_repeats.tsv"

    script:
    """
    source /opt/conda/etc/profile.d/conda.sh
    conda activate syri

    python ${annotate_breakpoints_script} \
        --syri ${syri_mapped} \
        --repeats1 ${repeats1} \
        --repeats2 ${repeats2} \
        --genes1 ${genes1} \
        --genes2 ${genes2}
    """
}

process PrepareIGVSessions {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path split_breakpoints_script
    path breakpoints_tsv
    path genes1
    path genes2
    path repeats1
    path repeats2
    path genome1_fa
    path genome2_fa
    path tad

    output:
    path "session_genome1.xml"
    path "session_genome2.xml"
    path "breakpoints_genome1.bed"
    path "breakpoints_genome2.bed"

    script:
    """
    source /opt/conda/etc/profile.d/conda.sh
    conda activate syri

    # Create FASTA indices for IGV
    samtools faidx ${genome1_fa}
    samtools faidx ${genome2_fa}

    # Split breakpoints by genome
    python ${split_breakpoints_script} ${breakpoints_tsv}

    # IGV session for genome1
    cat > session_genome1.xml << EOF
    <?xml version="1.0" encoding="UTF-8"?>
    <Session genome="${genome1_fa}" version="8">
        <Resources>
            <Resource path="./${genome1_fa}"/>
            <Resource path="./${genes1}"/>
            <Resource path="./${repeats1}"/>
            <Resource path="./${tad}"/>
            <Resource path="./breakpoints_genome1.bed"/>
        </Resources>
    </Session>
    EOF

    # IGV session for genome2
    cat > session_genome2.xml << EOF
    <?xml version="1.0" encoding="UTF-8"?>
    <Session genome="${genome2_fa}" version="8">
        <Resources>
            <Resource path="./${genome2_fa}"/>
            <Resource path="./${genes2}"/>
            <Resource path="./${repeats2}"/>
            <Resource path="./breakpoints_genome2.bed"/>
        </Resources>
    </Session>
    EOF
    """
}

process Visualize {
    publishDir "${params.outdir}/visualizer", mode: 'copy'

    input:
    path visualizer_script
    path structural_variants
    path genome1
    path genome2
    path repeats1
    path repeats2
    path genes1
    path genes2
    path config

    output:
    path "visualizer"

    script:
    """
    source /opt/conda/etc/profile.d/conda.sh
    conda activate python_sv

    mkdir visualizer

    echo -e "${genome1}\tgenome1" > genomes.txt
    echo -e "${genome2}\tgenome2" >> genomes.txt

    python ${visualizer_script} \
        --sr ${structural_variants} \
        --genomes genomes.txt \
        --ref-repeats ${repeats1} \
        --qry-repeats ${repeats2} \
        --ref-genes ${genes1} \
        --qry-genes ${genes2} \
        --mapping ${config}
    """
}

process SpectralTADAnalysis {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path spectraltad_script
    path hic

    output:
    path "all_tads.bed"

    shell:
    """
    source /opt/conda/etc/profile.d/conda.sh
    conda activate spectraltad

    mkdir -p spectraltad_results

    RES_LIST=\$(cooler ls ${hic} | awk -F'/' '{print \$NF}' | tr '\n' ' ')
    echo "Resolutions: \$RES_LIST"

    # Pick best resolution: ≤200kb, closest to 50kb
    MAIN_RES=50000
    BEST_RES=0
    for res in \$RES_LIST; do
        if [ \$res -le 200000 ]; then
            if [ \$res -ge \$MAIN_RES ]; then
                if [ \$((MAIN_RES - BEST_RES)) -gt \$((res - MAIN_RES)) ]; then
                    BEST_RES=\$res
                fi
                break
            fi
            BEST_RES=\$res
        fi
    done

    echo "BEST_RES: \$BEST_RES"

    cooler dump --join ${hic}::/resolutions/\$BEST_RES > spectraltad_matrix.txt

    Rscript ${spectraltad_script} \\
        --matrix spectraltad_matrix.txt \\
        --output spectraltad_results/

    # Combine all chromosomes into single BED
    tail -n +2 -q spectraltad_results/*.bed >> all_tads.bed
    """
}

process RunRepeatOBserver {
    publishDir "${params.outdir}/centromeres/histograms", pattern: "*/Summary_output/histograms/*", mode: 'copy'
    publishDir "${params.outdir}/centromeres/Shannon_div", pattern: "*/Summary_output/Shannon_div/*", mode: 'copy'
    publishDir "${params.outdir}/centromeres/spectra", pattern: "*/Summary_output/spectra/*", mode: 'copy'

    input:
    path genome

    output:
    path "output_chromosomes/SpeciesName_H0-AT/Summary_output/histograms/"
    path "output_chromosomes/SpeciesName_H0-AT/Summary_output/Shannon_div/"
    path "output_chromosomes/SpeciesName_H0-AT/Summary_output/spectra/"

    script:
    """
    source /opt/conda/etc/profile.d/conda.sh
    conda activate repeatobserver

    bash Setup_Run_Repeats.sh \\
        -i SpeciesName \\
        -f ${genome} \\
        -h H0 \\
        -c 15 \\
        -m 100000 \\
        -g FALSE
    """
}

process Compartments {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path hic
    path genes
    path genome

    output:
    path "compartments_tandem_genes.txt"

    script:
    """
    source /opt/conda/etc/profile.d/conda.sh
    conda activate hicsca

    # Detect tandem repeats
    /tools/trf-mod ${genome} > tandem_repeats.bed

    # Find best resolution (≤200kb, closest to 100kb)
    RES_LIST=\$(cooler ls ${hic} | awk -F'/' '{print \$NF}' | tr '\n' ' ')
    BEST_RES=0
    for res in \$RES_LIST; do
        if [ \$res -le 200000 ]; then
            if [ \$res -ge 100000 ]; then
                if [ \$((100000 - BEST_RES)) -gt \$((res - 100000)) ]; then
                    BEST_RES=\$res
                fi
                break
            fi
            BEST_RES=\$res
        fi
    done

    echo "Using resolution: \$BEST_RES"

    FIRST_RES=\$(echo \$RES_LIST | awk '{print \$1}')
    CHROMOSOMES=\$(cooler dump -t chroms "${hic}::/resolutions/\$FIRST_RES" | awk '\$2 >= 2000000 {printf "%s ", \$1}' | sed 's/ \$//')

    # Convert cooler to .hic for hic-sca compatibility
    hictk convert ${hic} hic.hic -r \$BEST_RES

    # Run compartment analysis
    hic-sca -f hic.hic -r \$BEST_RES -p hicsca_compartments -o compartmets_results -c \$CHROMOSOMES --bed
    tail -n +2 compartmets_results/hicsca_compartments_\${BEST_RES}bp.bed > compartments.bed
    awk '{print \$1, \$2, \$3, \$4}' OFS='\\t' compartments.bed > compartments_raw.bed

    # Calculate tandem repeat density per compartment
    bedtools coverage -a compartments_raw.bed -b tandem_repeats.bed | awk '{print \$1, \$2, \$3, \$4, \$8*100}' OFS='\\t' > compartments_tandems.txt

    # Calculate gene density per compartment
    bedtools coverage -a compartments_tandems.txt -b ${genes} | awk '{print \$1, \$2, \$3, \$4, \$5, \$6, \$9*100}' OFS='\\t' > compartments_tandem_genes.txt
    """
}

workflow {
    // Input validation
    genome1 = params.genome1 ? file(params.genome1) : null
    genome2 = params.genome2 ? file(params.genome2) : null
    config = params.config ? file(params.config) : null
    mcool_intragenomic = params.mcool_intragenomic ? file(params.mcool_intragenomic) : null
    mcool_intergenomic = params.mcool_intergenomic ? file(params.mcool_intergenomic) : null
    genes1_file = params.genes1 ? file(params.genes1) : null
    genes2_file = params.genes2 ? file(params.genes2) : null

    // Scripts (must be in workflow directory)
    compare_script = file("${projectDir}/compare_syri_eaglec2.py")
    extract_genes_script = file("${projectDir}/extract_genes.py")
    mapping_syri_results_script = file("${projectDir}/mapping_syri_results.py")
    annotate_breakpoints_script = file("${projectDir}/annotate_breakpoints.py")
    split_breakpoints_script = file("${projectDir}/split_breakpoints.py")
    visualizer_script = file("${projectDir}/visualizer.py")
    spectraltad_script = file("${projectDir}/spectraltad.R")

    repeats1 = genome1 ? Repeater2Genome1(genome1) : null
    repeats2 = genome2 ? Repeater2Genome2(genome2) : null

    synthetic_genomes = (genome1 && genome2 && config) ? CreateSyntheticGenomes(genome1, genome2, config) : null
    syri_out = synthetic_genomes ? SyRI(synthetic_genomes[0], synthetic_genomes[1]) : null

    syri_mapped = (syri_out && synthetic_genomes && genome1 && genome2) ? MappingSyriResults(mapping_syri_results_script, syri_out, synthetic_genomes[2], genome1, genome2) : null

    eaglec = mcool_intergenomic ? EagleC(mcool_intergenomic) : null

    comparison = (syri_mapped && eaglec) ? RunComparison(compare_script, syri_mapped, eaglec) : null
    syri_filtered = comparison ? comparison[0] : (syri_mapped ? syri_mapped : null)

    genes1 = genes1_file ? ExtractGenes1(extract_genes_script, genes1_file) : null
    genes2 = genes2_file ? ExtractGenes2(extract_genes_script, genes2_file) : null

    breakpoint_annotation = (syri_filtered && repeats1 && repeats2 && genome1 && genome2 && genes1 && genes2) ? AnnotateBreakpointsWithRepeats(annotate_breakpoints_script, syri_filtered, repeats1, repeats2, genome1, genome2, genes1, genes2) : null

    tad = mcool_intragenomic ? SpectralTADAnalysis(spectraltad_script, mcool_intragenomic) : null

    igv = (breakpoint_annotation && genes1 && genes2 && repeats1 && repeats2 && genome1 && genome2 && tad) ? PrepareIGVSessions(split_breakpoints_script, breakpoint_annotation, genes1, genes2, repeats1, repeats2, genome1, genome2, tad) : null

    visualize = (breakpoint_annotation && genome1 && genome2 && repeats1 && repeats2 && genes1 && genes2 && config) ? Visualize(visualizer_script, breakpoint_annotation, genome1, genome2, repeats1, repeats2, genes1, genes2, config) : null

    centromeres = genome1 ? RunRepeatOBserver(genome1) : null

    compartments = (mcool_intragenomic && genes1 && genome1) ? Compartments(mcool_intragenomic, genes1, genome1) : null
}

