import pandas as pd
import random


samples = pd.read_table("config.tsv").set_index("sample_name", drop=False).to_dict(orient='index')
tissue_names = list(samples.keys())

rule all:
    input:
        "plots/dimplots_tissues_1.png",
        "plots/dimplots_tissues_2.png",
        "plots/sil.pdf",
        "plots/knn.pdf",
        "plots/umap_remo_10x_B_lymphoma.png",
        "plots/umap_remo_brain_cancer.png",
        "plots/disease_metrics.pdf",
        "plots/whole_body_runtime.png",
        "plots/wholebody_dimplot.png",
        "plots/term_enrichment_umap.png",
        "plots/term_enrichment_metrics.pdf"

rule install_dependencies:
    input: "code/packages.txt"
    output: touch("install.done")
    shell:
        """
        # R packages
        Rscript -e 'setRepositories(ind=1:3); packages <- readLines({input}); install.packages(packages)'

        # Rust packages
        cargo install fragtk
        cargo install spars
        """

rule get_annotations:
    output: "data/annotations.rds"
    shell:
        """
        Rscript code/build_annotations.R
        """

rule get_remo:
    output: "data/REMOv1_GRCh38.bed.gz"
    shell: "wget -O {output} https://github.com/stuart-lab/REMO.v1.GRCh38/blob/main/inst/extdata/REMOv1_GRCh38.bed.gz"

rule download:
    input: "data/{tissue}/download.sh"
    output:
        "data/{tissue}/fragments.tsv.gz",
        "data/{tissue}/gex.mtx",
        "data/{tissue}/genes.tsv",
        "data/{tissue}/rna_cells.txt"
    shell:
        """
        cd data/{wildcards.tissue}
        sh download.sh
        """

rule call_peaks:
    input: "data/{tissue}/fragments.tsv.gz"
    output: 
        peaks="data/{tissue}/peaks.bed"
    shell:
        """
        macs2 callpeak \
            -f BED --nomodel --shift -100 --extsize 200 --name {wildcards.tissue} \
            -t {input} \
        	--outdir data/{wildcards.tissue}/

        # cut narrowPeak to  bed file
        cut -f1-3 data/{wildcards.tissue}/{wildcards.tissue}_peaks.narrowPeak > {output} 
        """

rule count_barcodes:
    input: "data/{tissue}/fragments.tsv.gz"
    output:
        barcode_counts="data/{tissue}/barcode_counts.tsv",
        barcodes="data/{tissue}/barcodes_atac.txt"
    params:
        ncells=10000
    shell:
        """
        fragtk count \
            -f {input} \
            -o {output.barcode_counts} \
            -n {params.ncells} \
            > {output.barcodes}
        """

rule remo_matrix:
    input:
        frags="data/{tissue}/fragments.tsv.gz",
        regions="data/REMOv1_GRCh38.bed.gz",
        barcodes="data/{tissue}/barcodes_atac.txt"
    output: directory("data/{tissue}/remo_multiome/")
    shell:
        """
        fragtk matrix \
        	-f {input.frags} \
        	-c {input.barcodes} \
            -o {output} \
        	-b {input.regions} \
        	--group \
            --pic
         """

rule ccre_matrix:
    input:
        frags="data/{tissue}/fragments.tsv.gz",
        regions="data/REMOv1_GRCh38.bed.gz",
        barcodes="data/{tissue}/barcodes_atac.txt"
    output: directory("data/{tissue}/ccre_multiome/")
    shell:
        """
        fragtk matrix \
        	-f {input.frags} \
        	-c {input.barcodes} \
            -o {output} \
        	-b {input.regions} \
            --pic
         """

rule peak_matrix:
    input:
        frags="data/{tissue}/fragments.tsv.gz",
        regions="data/{tissue}/peaks.bed",
        barcodes="data/{tissue}/barcodes_atac.txt"
    output: directory("data/{tissue}/peaks_multiome/")
    shell:
        """
        fragtk matrix \
        	-f {input.frags} \
        	-c {input.barcodes} \
            -o {output} \
        	-b {input.regions} \
            --pic
         """

rule build_linkpeaks_object:
    input:
        ccre=directory("data/{tissue}/ccre_multiome/"),
        gex="objects/{tissue}_multiome_gex.rds",
        annotation="data/annotations.rds",
        frags="data/{tissue}/fragments.tsv.gz"
    output: obj="objects/{tissue}_multiome_linkpeaks.rds"
    script: "code/build_linkpeaks_obj.R"

rule aggregate_linkpeaks:
    input: expand("objects/{tissue}_multiome_linkpeaks.rds", tissue=tissue_names)
    output: touch("linkpeaks.done")

rule QC_object:
    input:
        frags="data/{tissue}/fragments.tsv.gz",
        gex="data/{tissue}/gex.mtx",
        genes="data/{tissue}/genes.tsv",
        rna_cells="data/{tissue}/rna_cells.txt",
        peak_counts=directory("data/{tissue}/peaks_multiome/"),
        annotations="data/annotations.rds"
    output:
        qc_cells="data/{tissue}/cells_{tissue}_multiome.txt"
    params:
        nCount_ATAC_above=lambda wildcards: samples[wildcards.tissue]["nCount_ATAC_above"],
        nCount_ATAC_below=lambda wildcards: samples[wildcards.tissue]["nCount_ATAC_below"],
        TSS_above=lambda wildcards: samples[wildcards.tissue]["TSS_above"],
        nCount_RNA_above=lambda wildcards: samples[wildcards.tissue]["nCount_RNA_above"],
        nCount_RNA_below=lambda wildcards: samples[wildcards.tissue]["nCount_RNA_below"],
        percent_mt_below=lambda wildcards: samples[wildcards.tissue]["percent_mt_below"]
    script:
        "code/QC.R"

rule build_objects:
    input:
        gex="data/{tissue}/gex.mtx",
        genes="data/{tissue}/genes.tsv",
        rna_cells="data/{tissue}/rna_cells.txt",
        qc_cells="data/{tissue}/cells_{tissue}_multiome.txt",
        peak_counts=directory("data/{tissue}/peaks_multiome/"),
        remo_counts=directory("data/{tissue}/remo_multiome/"),
        annotations="data/annotations.rds",
        remo="data/REMOv1_GRCh38.bed.gz"
    output:
        remo="objects/{tissue}_multiome_remo.rds",
        peaks="objects/{tissue}_multiome_peaks.rds",
        gex="objects/{tissue}_multiome_gex.rds"
    params:
        dims=lambda wildcards: samples[wildcards.tissue]["dims"]
    script:
        "code/analyze.R"

rule download_pbmc_annotated:
    output: "objects/pbmc_10k_v3.rds"
    shell: "curl https://signac-objects.s3.amazonaws.com/pbmc_10k_v3.rds -o {output}"

rule predict_celltype_pbmc:
    input:
        ref="objects/pbmc_10k_v3.rds",
        obj="objects/pbmc_multiome_gex.rds"
    output: "data/pbmc/multiome_labels.rds"
    script:
        "code/pbmc/label_transfer.R"

rule calculate_metrics:
    input:
        remo="objects/{tissue}_multiome_remo.rds",
        peaks="objects/{tissue}_multiome_peaks.rds"
    output:
        metrics="data/metrics/{tissue}_metrics.rds"
    params:
        dims=lambda wildcards: samples[wildcards.tissue]["dims"],
        tissue_name=lambda wildcards: wildcards.tissue
    script:
        "code/metrics/calculate_metrics.R"

rule plot_metric:
    input:
        pbmc="data/metrics/pbmc_metrics.rds",
        brain="data/metrics/brain_metrics.rds",
        jejunum="data/metrics/jejunum_metrics.rds",
        heart="data/metrics/heart_metrics.rds",
        left_colon="data/metrics/left_colon_metrics.rds",
        liver="data/metrics/liver_metrics.rds",
        lung="data/metrics/lung_metrics.rds",
        muscle="data/metrics/muscle_metrics.rds",
        pancreas="data/metrics/pancreas_metrics.rds",
        bile_duct="data/metrics/bile_duct_metrics.rds",
        fallopian_tube="data/metrics/fallopian_tube_metrics.rds",
        grey_matter="data/metrics/grey_matter_metrics.rds",
        kidney_cortex="data/metrics/kidney_cortex_metrics.rds",
        pituitary="data/metrics/pituitary_metrics.rds",
        placenta="data/metrics/placenta_metrics.rds",
        retina="data/metrics/retina_metrics.rds"
    output:
        all_metrics="data/metrics/tissue_metrics.rds",
        all_sil="data/metrics/tissue_metrics_sil.rds",
        all_knn="data/metrics/tissue_metrics_knn.rds",
        plots_sil="plots/sil.pdf",
        plots_ch="plots/ch.pdf",
        plots_knn="plots/knn.pdf",
        plots_combined="plots/metrics_boxplot.pdf"
    script:
        "code/metrics/plot_metrics_all.R"

rule plot_umap:
    input:
        objects=expand("objects/{tissue}_multiome_remo.rds", tissue=tissue_names),
        pbmc_labels="data/pbmc/multiome_labels.rds"
    output:
        umap_1="plots/dimplots_tissues_1.png",
        umap_2="plots/dimplots_tissues_2.png",
        pbmc_umap_peaks="plots/umap_peaks_pbmc.png",
        pbmc_umap_remo="plots/umap_remo_pbmc.png",
        pbmc_umap_gex="plots/umap_gex_pbmc.png"
    params:
        tissues=tissue_names
    script:
        "code/metrics/plot_umap.R"

## Module stats -- fig 1B, 1C, 1D, 2A, 2B

rule module_stats:
    input:
        directory("data/pbmc/remo_multiome/"),
        directory("data/pbmc/peaks_multiome/"),
        directory("data/pbmc/ccre_multiome/"),
        "data/pbmc/gex.mtx",
        "data/pbmc/cells_pbmc_multiome.txt"
    output:
        "plots/fig1/module_size_distance.pdf",
        "plots/fig1/module_distance.pdf",
        "plots/fig1/module_coaccess.pdf",
        "plots/fig2/ncount_density.pdf",
        "plots/fig2/matrix_count_var.pdf",
        "plots/fig2/nonzero_mean.png",
        "plots/fig2/mean_variance.png"
    shell:
        """
        Rscript code/stats/remo_stats.R
        """

## Whole body dataset downsampling
N_SAMPLES = range(50000, 600001, 50000)
TISSUES = ['fetal', 'zhang']
SEED = 42

def sample_cells(input_file, output_file, N, seed=42):
    random.seed(seed)
    with open(input_file) as f:
        lines = f.readlines()
    
    sampled_lines = random.sample(lines, N)  # Random sampling without replacement
    
    with open(output_file, "w") as f:
        f.writelines(sampled_lines)


# Fetal and adult chromatin atlas data aggregation: https://github.com/stuart-lab/human-chromatin-atlas

rule get_zhang:
    output: "data/zhang/fragments.tsv.gz", "data/zhang/cells.txt", "data/zhang/peaks.bed", "data/zhang/counts.tsv.gz"
    shell:
        """
        aws s3 cp s3://stuartlab/zhang/fragments.tsv.gz data/zhang/ --request-payer
        aws s3 cp s3://stuartlab/zhang/fragments.tsv.gz.tbi data/zhang/ --request-payer
        aws s3 cp s3://stuartlab/zhang/peaks.bed data/zhang/ --request-payer
        aws s3 cp s3://stuartlab/zhang/cells.txt.gz data/zhang/ --request-payer
        gzip -d cells.txt.gz
        """

rule get_fetal:
    output: "data/fetal/fragments.tsv.gz", "data/fetal/cells.txt", "data/fetal/peaks.bed"
    shell:
        """
        aws s3 cp s3://stuartlab/fetal/fragments.tsv.gz data/fetal/ --request-payer
        aws s3 cp s3://stuartlab/fetal/fragments.tsv.gz.tbi data/fetal/ --request-payer
        aws s3 cp s3://stuartlab/fetal/peaks.bed.gz data/fetal/ --request-payer
        aws s3 cp s3://stuartlab/fetal/cells.txt.gz data/fetal/ --request-payer
        gzip -d cells.txt.gz
        gzip -d peaks.bed.gz
        """

rule downsample_cells:
    input: "data/{tissue}/cells.txt"
    output: "data/{tissue}/sampled_{N}.txt"
    params:
        N=lambda wildcards: int(wildcards.N),
        seed=SEED
    run:
        sample_cells(input[0], output[0], params.N, params.seed)

rule combine_peaks:
    input: "data/fetal/peaks.bed", "data/zhang/peaks.bed"
    output: "data/atlas_peaks.bed"
    script: "code/zhang/combine_peaks.R"

rule quantify_wholebody_peaks:
    input:
        frags="data/{tissue}/fragments.tsv.gz",
        cells="data/{tissue}/cells.txt",
        peaks="data/atlas_peaks.bed"
    output: directory("data/{tissue}/fullcounts_peaks")
    shell:
        """
        fragtk matrix -f {input.frags} -c {input.cells} -b {input.peaks} -o {output} --pic
        """

rule quantify_wholebody_remo:
    input:
        frags="data/{tissue}/fragments.tsv.gz",
        cells="data/{tissue}/cells.txt",
        regions="data/REMOv1_GRCh38.bed.gz"
    output: directory("data/{tissue}/fullcounts_remo")
    shell:
        """
        fragtk matrix -f {input.frags} -c {input.cells} -b {input.regions} -o {output} --pic --group
        """

rule remo_stats:
    input: directory("data/{tissue}/fullcounts_remo")
    output: "data/{tissue}/stats/remo_row.tsv", "data/{tissue}/stats/remo_col.tsv"
    shell:
        """
        spars stats -i {input}/matrix.mtx.gz -o data/{wildcards.tissue}/stats/remo
        """

rule select_features:
    input: expand("data/{tissue}/stats/remo_row.tsv", tissue=TISSUES)
    output: features="data/atlas/remo_features.txt"
    script: "code/zhang/select_features.R"

rule subset_counts_peaks:
    input:
        counts=directory("data/{tissue}/fullcounts_peaks"),
        cells="data/{tissue}/sampled_{N}.txt"
    output: directory("data/{tissue}/peaks_{N}")
    shell:
        """
        spars subset -i {input.counts} --cols {input.cells} -o {output} --colnames
        """

rule subset_counts_remo:
    input:
        counts=directory("data/{tissue}/fullcounts_remo"),
        cells="data/{tissue}/sampled_{N}.txt",
        features="data/atlas/remo_features.txt"
    output: directory("data/{tissue}/remo_{N}")
    shell:
        """
        spars subset -i {input.counts} --cols {input.cells} --rows {input.features} -o {output} --colnames
        """

# full dataset for each
rule subset_features_remo:
    input:
        counts=directory("data/{tissue}/fullcounts_remo"),
        features="data/atlas/remo_features.txt"
    output: directory("data/{tissue}/top_features_remo")
    shell:
        """
        spars subset -i {input.counts} --rows {input.features} -o {output}
        """

rule process_downsamplings_peaks:
    input:
        lambda wildcards: expand("data/{tissue}/peaks_{N}", tissue=TISSUES, N=wildcards.N)
    output: object="data/atlas/objects/ds_peaks_{N}.rds"
    params: dims=50
    benchmark: "benchmarks/atlas/peaks_{N}.txt"
    script:
        "code/zhang/process_peaks.R"

rule process_downsamplings_remo:
    input:
        lambda wildcards: expand("data/{tissue}/remo_{N}", tissue=TISSUES, N=wildcards.N)
    output: object="data/atlas/objects/ds_remo_{N}.rds"
    params: dims=50
    benchmark: "benchmarks/atlas/remo_{N}.txt"
    script:
        "code/zhang/process_remo.R"

rule plot_wholebody_runtime:
    input: runtimes=expand("benchmarks/atlas/{method}_{N}.txt", method=['peaks', 'remo'], N=N_SAMPLES)
    output: plot="plots/whole_body_runtime.png"
    script:
        "code/zhang/plot_runtime.R"

rule plot_wholebody:
    input: emb="data/atlas/objects/ds_remo_600000.rds"
    output: plot="plots/wholebody_dimplot.png"
    script:
        "code/zhang/plot_zhang.R"

## disease
brain_cancer_sample = ["I2", "I3", "I4", "I5", "M8", "M7"]

rule download_brain_cancer:
    input: "data/brain_cancer/download.sh"
    output:
        expand("data/brain_cancer/{sample}/fragments.tsv.gz", sample=brain_cancer_sample),
        "data/brain_cancer/GSE206579_snATAC_PF_EPN_Seurat.rds"
    shell:
        """
        cd data/brain_cancer
        sh download.sh
        """

rule barcodes_brain_cancer:
    input:
        source_obj="data/brain_cancer/GSE206579_snATAC_PF_EPN_Seurat.rds"
    output:
        peaks_bed="data/brain_cancer/peaks.bed",
        barcodes=expand("data/brain_cancer/{sample}/{sample}_barcodes.txt", sample=brain_cancer_sample)
    script:
        "code/disease/get_barcodes_brain_cancer.R"

rule matrix_brain_cancer:
    input:
        frags="data/brain_cancer/{sample}/fragments.tsv.gz",
        peaks_bed="data/brain_cancer/peaks.bed",
        remo_bed="data/REMOv1_GRCh38.bed.gz",
        barcodes="data/brain_cancer/{sample}/{sample}_barcodes.txt"
    output:
        remo=directory("data/brain_cancer/{sample}/remo"),
        peaks=directory("data/brain_cancer/{sample}/peaks")
    shell:
        """
        fragtk matrix \
        	-f {input.frags} \
        	-c {input.barcodes} \
            -o {output.remo} \
        	-b {input.remo_bed} \
        	--group \
            --pic

        fragtk matrix \
        	-f {input.frags} \
        	-c {input.barcodes} \
            -o {output.peaks} \
        	-b {input.peaks_bed} \
            --pic
         """

rule analyze_brain_cancer:
    input:
        source_obj="data/brain_cancer/GSE206579_snATAC_PF_EPN_Seurat.rds",
        remo_matrices=directory(expand("data/brain_cancer/{sample}/remo", sample=brain_cancer_sample)),
        peak_matrices=directory(expand("data/brain_cancer/{sample}/peaks", sample=brain_cancer_sample))
    output:
        remo_obj="objects/brain_cancer_remo.rds",
        peaks_obj="objects/brain_cancer_peaks.rds",
        umap_remo="plots/umap_remo_brain_cancer.png",
        umap_peaks="plots/umap_peaks_brain_cancer.png"
    script:
        "code/disease/analyze_brain_cancer.R"

rule get_B_lymphoma_annotations:
    output: "data/10x_B_lymphoma/10x_B_Lymphoma_Cell_Types.csv"
    shell:
        """
        # celltype annotations from loupe file 
        # https://cf.10xgenomics.com/samples/cell-arc/2.0.0/lymph_node_lymphoma_14k/lymph_node_lymphoma_14k_cloupe.cloupe
        aws s3 cp s3://stuartlab/REMO/disease/10x_B_Lymphoma_Cell_Types.csv data/10x_B_lymphoma/10x_B_Lymphoma_Cell_Types.csv
        """

rule annotate_B_lymphoma:
    input:
        gex="objects/10x_B_lymphoma_multiome_gex.rds",
        remo="objects/10x_B_lymphoma_multiome_remo.rds",
        peaks="objects/10x_B_lymphoma_multiome_peaks.rds",
        celltype_annot="data/10x_B_lymphoma/10x_B_Lymphoma_Cell_Types.csv"
    output:
        umap_gex="plots/umap_gex_10x_B_lymphoma.png",
        umap_remo="plots/umap_remo_10x_B_lymphoma.png",
        umap_peaks="plots/umap_peaks_10x_B_lymphoma.png"
    script:
        "code/disease/annotate_B_lymphoma.R" 

rule disease_metrics:
    input:
        B_lymphoma_remo="objects/10x_B_lymphoma_multiome_remo.rds",
        B_lymphoma_peaks="objects/10x_B_lymphoma_multiome_peaks.rds",
        brain_cancer_remo="objects/brain_cancer_remo.rds",
        brain_cancer_peaks="objects/brain_cancer_peaks.rds",
        celltype_annot="data/10x_B_lymphoma/10x_B_Lymphoma_Cell_Types.csv"
    output:
        B_lymphoma_metrics="data/metrics/10x_B_lymphoma_metrics.rds",
        brain_cancer_metrics="data/metrics/brain_cancer_metrics.rds",
        disease_metrics_plot="plots/disease_metrics.pdf"
    script:
        "code/disease/disease_metrics.R"

## Term enrichment
rule annotate_brain:
    input: obj="objects/brain_multiome_gex.rds"
    output: obj="data/term_enrichment/brain/brain_gex.rds"
    script: "code/term_enrichment/brain/annotate_celltypes.R"

rule brain_term_enrichment:
    input:
        gex="data/term_enrichment/brain/brain_gex.rds",
        remo="objects/brain_multiome_remo.rds"
    output:
        remo="data/term_enrichment/brain/brain_remo.rds",
        predictions="data/term_enrichment/brain/brain_predictions.rds",
        organ="data/term_enrichment/brain/brain_predictions_organ.rds"
    script: "code/term_enrichment/brain/make_remo_obj.R"

rule download_pbmc_term_enrichment:
    output:
        "data/term_enrichment/pbmc/fragments.tsv.gz",
        "data/term_enrichment/pbmc/gex.mtx",
        "data/term_enrichment/pbmc/genes.tsv",
        "data/term_enrichment/pbmc/rna_cells.txt"
    shell:
        """
        cd data/term_enrichment/pbmc/
        sh download.sh
        """

rule qc_pbmc:
    input:
        frags="data/term_enrichment/pbmc/fragments.tsv.gz",
        gex="data/term_enrichment/pbmc/gex.mtx",
        genes="data/term_enrichment/pbmc/genes.tsv",
        cells="data/term_enrichment/pbmc/rna_cells.txt",
        annotations="data/annotations.rds"
    output: cells="data/term_enrichment/pbmc/cells_pbmc_multiome.txt"
    script: "code/term_enrichment/pbmc/QC.R"

rule quant_pbmc_remo:
    input:
        frags="data/term_enrichment/pbmc/fragments.tsv.gz",
        cells="data/term_enrichment/pbmc/cells_pbmc_multiome.txt",
        remo="data/REMOv1_GRCh38.bed.gz"
    output: directory("data/term_enrichment/pbmc/remo")
    shell:
        """
        fragtk matrix \
            -f {input.frags} \
            -c {input.cells} \
            -o data/term_enrichment/pbmc/remo/ \
            -b {input.remo} \
            --pic \
            --group
        """

rule build_pbmc_multiome_term_enrichment:
    input:
        remo=directory("data/term_enrichment/pbmc/remo"),
        allcells="data/term_enrichment/pbmc/rna_cells.txt",
        cells="data/term_enrichment/pbmc/cells_pbmc_multiome.txt",
        gex="data/term_enrichment/pbmc/gex.mtx",
        genes="data/term_enrichment/pbmc/genes.tsv",
        ref="objects/pbmc_10k_v3.rds"
    output: obj="data/term_enrichment/pbmc/pbmc_multiome.rds"
    script: "code/term_enrichment/pbmc/make_multiome_obj.R"

rule term_enrichment_pbmc:
    input:
        multiome="data/term_enrichment/pbmc/pbmc_multiome.rds",
        remo=directory("data/term_enrichment/pbmc/remo")
    output:
        remo="data/term_enrichment/pbmc/pbmc_remo.rds",
        predictions="data/term_enrichment/pbmc/pbmc_predictions.rds",
        organ="data/term_enrichment/pbmc/pbmc_predictions_organ.rds"
    script:
        "code/term_enrichment/pbmc/make_remo_obj.R"

rule download_islet:
    output:
        "data/term_enrichment/islet/fragments.tsv.gz",
        "data/term_enrichment/islet/gex.mtx",
        "data/term_enrichment/islet/genes.tsv",
        "data/term_enrichment/islet/rna_cells.txt"
    shell:
        """
        cd data/term_enrichment/islet
        aws s3 cp s3://stuartlab/REMO/annotate/islet/ ./ --recursive
        mv atac_fragments.tsv.gz fragments.tsv.gz
        mv atac_fragments.tsv.gz.tbi fragments.tsv.gz.tbi
        Rscript ../../../code/extract_matrix.R filtered_feature_bc_matrix.h5 ./
        """

rule qc_islet:
    input:
        frags="data/term_enrichment/islet/fragments.tsv.gz",
        gex="data/term_enrichment/islet/gex.mtx",
        genes="data/term_enrichment/islet/genes.tsv",
        cells="data/term_enrichment/islet/rna_cells.txt",
        annotations="data/annotations.rds"
    output: cells="data/term_enrichment/islet/cells_islet_multiome.txt"
    script: "code/term_enrichment/islet/QC.R"

rule quant_islet_remo:
    input:
        frags="data/term_enrichment/islet/fragments.tsv.gz",
        cells="data/term_enrichment/islet/cells_islet_multiome.txt",
        remo="data/REMOv1_GRCh38.bed.gz"
    output: directory("data/term_enrichment/islet/remo")
    shell:
        """
        fragtk matrix \
            -f {input.frags} \
            -c {input.cells} \
            -o data/term_enrichment/islet/remo/ \
            -b {input.remo} \
            --pic \
            --group
        """

rule build_islet_multiome:
    input:
        allcells="data/term_enrichment/islet/rna_cells.txt",
        cells="data/term_enrichment/islet/cells_islet_multiome.txt",
        gex="data/term_enrichment/islet/gex.mtx",
        genes="data/term_enrichment/islet/genes.tsv"
    output: obj="data/term_enrichment/islet/islet_multiome.rds"
    script: "code/term_enrichment/islet/make_multiome_obj.R"

rule term_enrichment_islet:
    input:
        multiome="data/term_enrichment/islet/islet_multiome.rds",
        remo=directory("data/term_enrichment/islet/remo")
    output:
        remo="data/term_enrichment/islet/islet_remo.rds",
        predictions="data/term_enrichment/islet/islet_predictions.rds",
        organ="data/term_enrichment/islet/islet_predictions_organ.rds"
    script:
        "code/term_enrichment/islet/make_remo_obj.R"

rule term_enrichment_metrics:
    input:
        brain_remo="data/term_enrichment/brain/brain_remo.rds",
        brain_predictions="data/term_enrichment/brain/brain_predictions.rds",
        pbmc_remo="data/term_enrichment/pbmc/pbmc_remo.rds",
        pbmc_predictions="data/term_enrichment/pbmc/pbmc_predictions.rds",
        islet_remo="data/term_enrichment/islet/islet_remo.rds",
        islet_predictions="data/term_enrichment/islet/islet_predictions.rds",
        label_match="data/term_enrichment/label_match.csv"
    output:
        umap="plots/term_enrichment_umap.png",
        accuracy_score="plots/term_enrichment_metrics.pdf",
        wholebody="data/term_enrichment/fgsea_wholebody.csv",
        tissue="data/term_enrichment/fgsea_tissue.csv"
    script:
        "code/term_enrichment/term_enrichment_metrics.R"
