nextflow.enable.dsl = 2

params.genome1 = ""
params.genome2 = ""
params.config = ""
params.hic = ""
params.genes1 = ""
params.genes2 = ""

process Repeater2Genome1 {
    input:
    path genome

    output:
    path "repeats1.bed"

    script:
    """
    source /nfs/home/tklimentiev/miniconda3/etc/profile.d/conda.sh
    conda activate nextflow_env

    java -jar -Xms16g -Xmx64g ${projectDir}/../repeater2/dist/Repeater2.jar \\
        ${genome} kmer=20 sln=200 -seqshow

    cat *.gff | \\
        awk -F'\\t' '/^[^#]/ && \$4 ~ /^[0-9]+\$/ {
            split(\$1, a, " ")
            print a[1], \$4-1, \$5, \$3
        }' OFS='\\t' > repeats1.bed
    """
}

process Repeater2Genome2 {
    input:
    path genome

    output:
    path "repeats2.bed"

    script:
    """
    source /nfs/home/tklimentiev/miniconda3/etc/profile.d/conda.sh
    conda activate nextflow_env

    java -jar -Xms16g -Xmx64g ${projectDir}/../repeater2/dist/Repeater2.jar \\
        ${genome} kmer=20 sln=200 -seqshow

    cat *.gff | \\
        awk -F'\\t' '/^[^#]/ && \$4 ~ /^[0-9]+\$/ {
            split(\$1, a, " ")
            print a[1], \$4-1, \$5, \$3
        }' OFS='\\t' > repeats2.bed
    """
}

process CreateSyntheticGenomes {
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
    source /nfs/home/tklimentiev/miniconda3/etc/profile.d/conda.sh
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

        source_chromosomes="\${line%:*}"
        target_genomes="\${line#*:}"

        echo ">chrom\${chrom_num}" >> new_genome1.fa
        genome1_lengths=""
        genome1_total=0
        for chr in \$(echo "\$source_chromosomes" | tr ',' ' '); do
            awk -v target="\$chr" '\$1 == ">"target {flag=1; next} /^>/{flag=0} flag{print}' ${genome1} >> new_genome1.fa
            chr_length=\$(seqkit fx2tab -n -l ${genome1} | awk -v target="\$chr" '\$1 == target {print \$NF}')
            genome1_lengths="\$genome1_lengths,\$chr_length"
            genome1_total=\$((genome1_total + chr_length))
        done
        genome1_lengths="\${genome1_lengths:1}"

        echo ">chrom\${chrom_num}" >> new_genome2.fa
        genome2_lengths=""
        genome2_total=0
        for chr in \$(echo "\$target_genomes" | tr ',' ' '); do
            awk -v target="\$chr" '\$1 == ">"target {flag=1; next} /^>/{flag=0} flag{print}' ${genome2} >> new_genome2.fa
            chr_length=\$(seqkit fx2tab -n -l ${genome2} | awk -v target="\$chr" '\$1 == target {print \$NF}')
            genome2_lengths="\$genome2_lengths,\$chr_length"
            genome2_total=\$((genome2_total + chr_length))
        done
        genome2_lengths="\${genome2_lengths:1}"

        # Write mapping info
        echo "chrom\$chrom_num\t\$source_chromosomes\t\$target_genomes\t\$genome1_lengths\t\$genome2_lengths" >> chromosome_mapping.tsv

        ((chrom_num++))
    done < ${config}
    """
}

process SyRI {
    input:
    path genome1
    path genome2

    output:
    path "syri.out"

    script:
    """
    source /nfs/home/tklimentiev/miniconda3/etc/profile.d/conda.sh
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
    input:
    path syri
    path chromosome_mapping
    path genome1
    path genome2

    output:
    path "mapped_svs.tsv"

    script:
    """
    source /nfs/home/tklimentiev/miniconda3/etc/profile.d/conda.sh
    conda activate syri

    python3 ${projectDir}/mapping_syri_results.py ${syri} ${chromosome_mapping} ${genome1} ${genome1}
    """
}

process EagleC {
    input:
    path hic
    path genome

    output:
    path "eaglec_results/output.SV_calls.txt"

    shell:
    '''
    source /nfs/home/tklimentiev/miniconda3/etc/profile.d/conda.sh
    conda activate EagleC2

    echo !{projectDir}

    # Get resolutions and convert to comma-separated list
    RES_LIST=$(cooler ls !{hic} | grep -oP '\\d+' | sort -n | tr '\n' ',' | sed 's/,$//')
    echo "Resolutions list: $RES_LIST"

    # Get first resolution for chromosome filtering
    FIRST_RES=$(echo "$RES_LIST" | cut -d',' -f1)

    # Get chromosomes as space-separated list
    CHROMOSOMES=$(cooler dump -t chroms "!{hic}::/resolutions/$FIRST_RES" 2>/dev/null | \\
      awk '$2 >= 2000000 && $1 !~ /(ChrUn|random|scaffold|mito|Unknown)/ {printf "%s ", $1}' | sed 's/ $//')

    echo "Chromosomes: $CHROMOSOMES"

    # Calculate intra and inter extend sizes based on resolutions
    INTRA_EX=""
    INTER_EX=""

    # Convert RES_LIST to array for processing
    IFS=',' read -ra RES_ARRAY <<< "$RES_LIST"
    for res in "${RES_ARRAY[@]}"; do
        if [ $res -le 25000 ]; then
            INTRA_EX+="3,"
            INTER_EX+="2,"
        elif [ $res -le 100000 ]; then
            INTRA_EX+="2,"
            INTER_EX+="1,"
        else
            INTRA_EX+="1,"
            INTER_EX+="1,"
        fi
    done

    INTRA_EX=${INTRA_EX%,}
    INTER_EX=${INTER_EX%,}

    echo "Intra extend sizes: $INTRA_EX"
    echo "Inter extend sizes: $INTER_EX"

    mkdir -p eaglec_results
    predictSV --mcool !{hic} \\
      --resolutions "$RES_LIST" \\
      --intra-extend-size "$INTRA_EX" \\
      --inter-extend-size "$INTER_EX" \\
      --model-path !{projectDir}/../EagleC2-models \\
      -g other \\
      -C $CHROMOSOMES \\
      -p 16 --prob-cutoff-1 0.05 \\
      -O eaglec_results/output
    '''
}

process RunComparison {
    input:
    path syri
    path eaglec2

    output:
    path "matches.tsv"
    path "comparison_stats.txt"

    script:
    """
    source /nfs/home/tklimentiev/miniconda3/etc/profile.d/conda.sh
    conda activate python_sv

    echo "SyRI file: $syri"
    echo "EagleC2 file: $eaglec2"

    python ${projectDir}/compare_syri_eaglec2.py \
        --syri "$syri" \
        --eaglec2 "$eaglec2" \
        --output "matches.tsv" \
        --sample "rice" \
        --pad 100000 > comparison_stats.txt
    """
}

process AnnotateBreakpointsWithRepeats {
    input:
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
    source /nfs/home/tklimentiev/miniconda3/etc/profile.d/conda.sh
    conda activate syri
    
    python ${projectDir}/annotate_breakpoints.py \
        --syri ${syri_mapped} \
        --repeats1 ${repeats1} \
        --repeats2 ${repeats2} \
        --genes1 ${genes1} \
        --genes2 ${genes2}
    """
}

workflow {
    genome1 = file(params.genome1)
    genome2 = file(params.genome2)
    config = file(params.config)
    hic_file = file(params.hic)

    repeats1 = Repeater2Genome1(genome1)
    repeats2 = Repeater2Genome2(genome2)

    synthetic_genomes = CreateSyntheticGenomes(genome1, genome2, config)
    syri_out = SyRI(synthetic_genomes[0], synthetic_genomes[1])

    syri_mapped = MappingSyriResults(syri_out, synthetic_genomes[2], genome1, genome2)

    eaglec = EagleC(hic_file, genome1)

    RunComparison(syri_mapped, eaglec)

    genes1 = file(params.genes1)
    genes2 = file(params.genes2)

    breakpoint_annotation = AnnotateBreakpointsWithRepeats(
        syri_mapped,
        repeats1,
        repeats2,
        genome1,
        genome2,
        genes1,
        genes2
    )
}
