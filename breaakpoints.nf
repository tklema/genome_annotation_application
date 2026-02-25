nextflow.enable.dsl = 2

params.genome1 = ""
params.genome2 = ""
params.repeats1 = ""
params.repeats2 = ""
params.syri = ""
params.genes1 = ""
params.genes2 = ""

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

process InferDirectionality {
    input:
    path breakpoint_annotation
    path syri_mapped

    output:
    path "directional_svs.tsv"

    script:
    """
    source /nfs/home/tklimentiev/miniconda3/etc/profile.d/conda.sh
    conda activate syri

    python ${projectDir}/infer_direction.py \
        --annotation ${breakpoint_annotation} \
        --syri ${syri_mapped}
    """
}

workflow {
    genome1 = file(params.genome1)
    genome2 = file(params.genome2)

    repeats1 = file(params.repeats1)
    repeats2 = file(params.repeats2)

    syri_mapped = file(params.syri)

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

    InferDirectionality(breakpoint_annotation, syri_mapped)
}
