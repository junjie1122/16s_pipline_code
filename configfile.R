#============================================
all_out_dir = "/mnt/e/微生物/16s_my/output/analysis_result_0927_v3/"
scripts_dir = "/mnt/e/微生物/16s_my/scripts"
sample_group_file = "/mnt/e/微生物/16s_my/need_table/sample_group.txt"
meta_data_file = sample_group_file                                               # 为NULL时，就是使用sample_group_file
fq_dir = "/mnt/e/微生物/16s_my/fq/"
fq_suffix = "_1.fastq.gz,_2.fastq.gz"


# cutadapt的参数
primer_f = "ACTCCTACGGRAGGCAGCAG"
primer_r = "GGACTACHVGGGTWTCTAAT"


#dada2步骤的参数
trim_f = 0              # 从5'端去除的碱基数 知道引物序列 就设置0
trim_r = 0              # 从3'端去除的碱基数 知道引物序列 就设置0
trunc_f = 220        
trunc_r = 220


# 是否过滤rare ASV   计算特征表的reads数之和乘以0.00005，得到的值用于过滤rare ASV  https://github.com/olabiyi/snakemake-workflow-qiime2
# min_sample_percent = 0.2 表示过滤少于在20%的样本中出现的ASV
filter_rare = "T"
min_sample_percent = 0.2

# 已有的分类器路径
classifier_file = "/mnt/e/微生物/16s_my/database/silva-138-99-515-806-nb-classifier.qza"
# classifier_file = "/mnt/e/微生物/16s_my/database/silva-138-99-nb-classifier.qza"
# classifier_file ="/mnt/e/微生物/16s_my/database/silva138_AB_V3-V4_classifier.qza"


sepp_ref_database = "/mnt/e/微生物/16s_my/database/sepp-refs-silva-128.qza"


# lefse
lda = 3


# picrust2的参数
KEGG2KO_file <- "/mnt/e/微生物/16s_my/database/KEGG_pathways_to_KO.tsv"


# kegg通路差异分析
logfc =0.001
use_adjp = "F"



# alpha rarefaction
# max_depth = 9000         # alpha rarefaction的最大深度(选择大于大部分样本redas数的深度)从denoise文件夹内的state.qzv中查看

#============================================

#============================================
# qc_raw_dirname = "quality_check_raw"
# trim_galore_dirname = "trim_galore_results"
# clean_fq_dirname = "clean_fq"
# qc_clean_dirname = "quality_check_clean"


qiime_import_dirname = "qiime_import"
qiime_cutadapt_dirname = "qiime_cutadapt"
denoise_dirname = "denoise"
taxonomy_dirname = "taxonomy"
tree_dirname = "tree"
feature_table_seq_dirname = "feature_table_seq"
beta_diversity_dirname = "beta_diversity"
alpha_diversity_dirname = "alpha_diversity"
lefse_dirname = "lefse"
picrust2_dirname = "picrust2"


plot_dirname = "plot"

#============================================