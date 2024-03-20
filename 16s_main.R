
#============================================
# all_out_dir = "/mnt/e/微生物/16s_my/output/analysis_result1/"
# scripts_dir = "/mnt/e/微生物/16s_my/scripts"
# sample_group_file = "/mnt/e/微生物/16s_my/need_table/sample_group.txt"
# meta_data_file = sample_group_file                                               # 为NULL时，就是使用sample_group_file
# fq_dir = "/mnt/e/_new/project/痛风16s数据分析/fq/"
# fq_suffix = "_1.fastq.gz,_2.fastq.gz"


# # cutadapt的参数
# primer_f = "ACTCCTACGGRAGGCAGCAG"
# primer_r = "GGACTACHVGGGTWTCTAAT"


# #dada2步骤的参数
# trim_f = 0
# trim_r = 0
# trunc_f = 220        
# trunc_r = 220


# # 是否过滤rare ASV   计算特征表的reads数之和乘以0.00005，用于过滤 rare ASV  https://github.com/olabiyi/snakemake-workflow-qiime2
# filter_rare = T


# # 已有的分类器路径
# classifier_file = "/mnt/e/微生物/16s_my/database/silva-138-99-515-806-nb-classifier.qza"
# # classifier_file = "/mnt/e/微生物/16s_my/database/silva-138-99-nb-classifier.qza"
# # classifier_file ="/mnt/e/微生物/16s_my/database/silva138_AB_V3-V4_classifier.qza"

# sepp_ref_database = "/mnt/e/微生物/16s_my/database/sepp-refs-silva-128.qza"

# # picrust2的参数
# KEGG2KO_file <- "/mnt/e/微生物/16s_my/database/KEGG_pathways_to_KO.tsv"


# # alpha rarefaction
# # max_depth = 9000         # alpha rarefaction的最大深度(选择大于大部分样本redas数的深度)从denoise文件夹内的state.qzv中查看

# #============================================

# #============================================
# # qc_raw_dirname = "quality_check_raw"
# # trim_galore_dirname = "trim_galore_results"
# # clean_fq_dirname = "clean_fq"
# # qc_clean_dirname = "quality_check_clean"


# qiime_import_dirname = "qiime_import"
# qiime_cutadapt_dirname = "qiime_cutadapt"
# denoise_dirname = "denoise"
# taxonomy_dirname = "taxonomy"
# tree_dirname = "tree"
# feature_table_seq_dirname = "feature_table_seq"
# beta_diversity_dirname = "beta_diversity"
# alpha_diversity_dirname = "alpha_diversity"
# lefse_dirname = "lefse"
# picrust2_dirname = "picrust2"


# plot_dirname = "plot"

#============================================


# 接受configfile参数
suppressPackageStartupMessages(library(argparse))
docstring = "16s_main.R"
myparser <- ArgumentParser(description=docstring, formatter_class="argparse.RawTextHelpFormatter")
myparser$add_argument("-c","--configfile", default="configfile.R", help="配置文件")
args <- myparser$parse_args()
source(args$configfile,local = T)




print("start analysis.....................")
# 加载函数
fun_fs <- list.files(paste0(scripts_dir,"/function/"),full.names = T)
for (f in fun_fs) {
  source(fun_fs)
}

library(dplyr)



if (!dir.exists(all_out_dir)) {
  dir.create(all_out_dir,recursive = T)
}


# 在all_out_dir下生成manifest.txt
generate_manifestfile_str <- paste0("Rscript ",scripts_dir,"/generate_manifestfile.R "," ",fq_dir," ",fq_suffix," -f ",sample_group_file, " -o ",all_out_dir)
record_command(command = generate_manifestfile_str,outdir = all_out_dir)
system(generate_manifestfile_str)



# qiime 导入数据
# 如果不存在qiime_import_complete.log，则运行qiime_import.sh
qiime_import_outdir <- paste0(all_out_dir,"/",qiime_import_dirname,"/")
if (!file.exists(paste0(qiime_import_outdir,"/qiime_import_complete.log"))) {

  print("import fq to qiime2.....................")
  qiime_import_str <- paste0("bash ",scripts_dir,"/qiime_import.sh ",
                             all_out_dir,"/manifest.txt "," ",qiime_import_outdir)

  record_command(command = qiime_import_str,outdir = all_out_dir)
  system(qiime_import_str)
}else {
  print("skip qiime_import.sh")
}


# 使用cutadapt去除引物 
# 再qiime_cutadapt_dirname 生成trim_demuxed.qzv
qiime_cutadapt_outdir <- paste0(all_out_dir,"/",qiime_cutadapt_dirname,"/")
if (!file.exists(paste0(qiime_cutadapt_outdir,"/qiime_cutadapt_complete.log"))) {
  print("run cutadapt.....................")
  qiime_cutadapt_str <- paste0("bash ",scripts_dir,"/qiime_cutadapt.sh ",
                             qiime_import_outdir,"/demuxed.qza ",
                             primer_f," ",primer_r," ",
                             qiime_cutadapt_outdir)

  record_command(command = qiime_cutadapt_str,outdir = all_out_dir)
  system(qiime_cutadapt_str)
}else {
  print("skip qiime_cutadapt.sh")
}

# 没有引物时的处理
# if(primer_f == "" & primer_r == ""){


# }




# 上一步生成/dada2_trunc_param.txt文件，包含自动计算得到的trunc_f,trunc_r
# source(paste0(qiime_cutadapt_outdir,"/dada2_trunc_param.txt"),local = T)
# print(trunc_f)
# print(trunc_r)


# dada2去噪 生成特征表和特征序列  使用 cuttadapt后的数据trim_demuxed.qza
denoise_input_file <- paste0(qiime_cutadapt_outdir,"/trim_demuxed.qza")
denoise_outdir <- paste0(all_out_dir,"/",denoise_dirname,"/")

if (!file.exists(paste0(denoise_outdir,"/complete.log"))) {
  print("run dada2.....................")
  denoise_str <- paste0("bash ",scripts_dir,"/qiime_denoise.sh --input ",
                       denoise_input_file," --meta ",meta_data_file," --outdir ",denoise_outdir,
                       " --trim_f ",trim_f," --trim_r ",trim_r," --trunc_f ",trunc_f," --trunc_r ",trunc_r)
  record_command(command = denoise_str,outdir = all_out_dir)
  system(denoise_str)
}else {
  print("skip denoise.sh")
}


# 计算特征表的reads数之和乘以0.00005，用于过滤 rare ASV  https://github.com/olabiyi/snakemake-workflow-qiime2

if(filter_rare =="T"){

  feature_table_tsvfile <- paste0(denoise_outdir,"/export/feature-table.tsv")
  tmpdf <- read.table(feature_table_tsvfile,comment.char = "",skip = 1,sep = "\t",header = T,row.names = 1)

  min_sample <- ceiling(ncol(tmpdf)*min_sample_percent)

  filter_freq <- round((sum(rowSums(tmpdf))) * 0.00005)
  print(paste0("Excludes rare ASV, ASVs with sequences less than 0.005% of the total number of sequences:"))
  print(paste0("0.005% of the total number of sequences 的值为: ",filter_freq))
}else{
  filter_freq <- 0
  min_sample <- 0
}


if (!file.exists(paste0(denoise_outdir,"/filter_complete.log"))) {
  print("run filetr rare ASV.....................")
  filter_rare_str <- paste0("bash ",scripts_dir,"/filter_rare_feature.sh ",
                            denoise_outdir," ",
                            filter_freq)
  record_command(command = filter_rare_str,outdir = all_out_dir)
  system(filter_rare_str)
}else {
  print("skip filter_rare_feature.sh")
}


# 后续使用filter过后的数据 如果filter_rare == F 这些filter文件还是原始的数据



# 使用feature-classifier进行物种注释
feature_seq_file <- paste0(denoise_outdir,"/feature-seq_filter.qza")
taxonomy_outdir <- paste0(all_out_dir,"/",taxonomy_dirname,"/")

if (!file.exists(paste0(all_out_dir,"/taxonomy/complete.log"))) {

  print("run taxonomy.....................")
  taxonomy_str <- paste0("bash ",scripts_dir,"/qiime_feature_classifier.sh ",
                         feature_seq_file," ",classifier_file," ",
                         taxonomy_outdir)
  record_command(command = taxonomy_str,outdir = all_out_dir)
  system(taxonomy_str)
}else {
  print("qiime_feature_classifier.sh")
}



# 修饰物种注释表，没有注释的水平用上级有注释的添加在在尾部
taxonomy_file <- paste0(taxonomy_outdir,"/export/taxonomy.tsv")
taxonomy_modified_file <- paste0(taxonomy_outdir,"/export/taxonomy_modified.tsv")

if(!file.exists(taxonomy_modified_file)){
  modified_taxonomy_str <- paste0("Rscript ",scripts_dir,"/modified_taxonomy_table.R "," --input ",taxonomy_file," --output ",taxonomy_modified_file)
  record_command(command = modified_taxonomy_str,outdir = all_out_dir)
  system(modified_taxonomy_str)
}else {
  print("skip modified_taxonomy_table.R")
}




# 生成系统发生树
feature_seq_file <- paste0(denoise_outdir,"/feature-seq_filter.qza")
sepp_ref_database <- sepp_ref_database
tree_outdir <- paste0(all_out_dir,"/",tree_dirname,"/")



if (!file.exists(paste0(all_out_dir,"/tree/complete.log"))) {
  print("run generate phylogenetic tree .....................")
  tree_str <- paste0("bash ",scripts_dir,"/generate_phylogenetic_tree.sh ",
                     feature_seq_file," ",
                     tree_outdir)
  record_command(command = tree_str,outdir = all_out_dir)             
  system(tree_str)
}else {
  print("skip generate_phylogenetic_tree.sh")
} 



# 获取样本的最低深度
feature_table_tsvfile <- paste0(denoise_outdir,"/export/feature-table_filter.tsv")
tmpdf <- read.table(feature_table_tsvfile,comment.char = "",skip = 1,sep = "\t",header = T,row.names = 1)
min_depth_for_rare <- min(colSums(tmpdf))

# 获取样本的比较高的合适深度
# max_depth_for_rare <- round(quantile(colSums(tmpdf),probs = 0.8) %/% 100) *100  %>% as.numeric()




## 使用样本的最低深度进行抽平
feature_table_tsvfile <- paste0(denoise_outdir,"/export/feature-table_filter.tsv")
feature_table_rarefile <- paste0(denoise_outdir,"/export/feature-table_rare.txt")


if(!file.exists(feature_table_rarefile)){
  rare_feature_table_str <- paste0("Rscript ",scripts_dir,"/rare_feature_table.R ","--input ",
                                  feature_table_tsvfile," --output ",feature_table_rarefile)

  record_command(command = rare_feature_table_str,outdir = all_out_dir)
  system(rare_feature_table_str)
}else {
  print("skip rare_feature_table.R")
}



# seq文件也跟过滤0的OTU表的OTU匹配  （目前不需要,因为都大于0）


# 使用vagan包的decostand函数选择hell方法 对抽屏后的特征表进行标准化
feature_table_normfile <- paste0(denoise_outdir,"/export/feature-table_norm.txt")

if(!file.exists(feature_table_normfile)){
  norm_feature_table_str <- paste0("Rscript ",scripts_dir,"/norm_feature_table.R ","--input ",
                                feature_table_rarefile," --output ",feature_table_normfile)                                
  record_command(command = norm_feature_table_str,outdir = all_out_dir)
  system(norm_feature_table_str)
}else {
  print("skip norm_feature_table.R")
}




# 抽平的的OTU表，在各个物种水平汇总每个样本中OTU的数量
# feature_table_rarefile
# taxonomy_modified_file
# feature_table_seq_outdir



summary_feature_table_dir <- paste0(denoise_outdir,"/export/summary_feature_table/")

if(!dir.exists(summary_feature_table_dir)){
  dir.create(summary_feature_table_dir,recursive = T)

  summary_feature_table_2_high_level_str <- paste0("Rscript ",scripts_dir,"/summary_feature_table_2_high_level.R ","--feature_table ",
                                feature_table_rarefile," --taxon_table ",taxonomy_modified_file," --outdir ",summary_feature_table_dir)

  record_command(command = summary_feature_table_2_high_level_str,outdir = all_out_dir)
  system(summary_feature_table_2_high_level_str)

}else {
  print("skip summary_feature_table_2_high_level.R")
}



#  β多样性分析---Beta Diversity  Anosim PCA PCOA NMDS 使用 norm的OTU表进行分析
beta_diversity_outdir <- paste0(all_out_dir,"/",beta_diversity_dirname,"/")
if (!dir.exists(beta_diversity_outdir)){
  dir.create(beta_diversity_outdir,recursive = T)
}
beta_diversity_outdir_log <- paste0(beta_diversity_outdir,"/log/")
if(!dir.exists(beta_diversity_outdir_log)){
  
  dir.create(beta_diversity_outdir_log,recursive = T)
}


## 绘制PCA
feature_table_normfile
pca_box_plot_out_file <- paste0(beta_diversity_outdir,"/","beta_diversity_pca_plot.pdf")


if (!file.exists(paste0(beta_diversity_outdir_log,"/beta_diversity_pca.complete.log"))) {

  print("run bate diversity analysis.....................")
  pca_box_plot_str <- paste0("Rscript ",scripts_dir,"/beta_diversity_pca.R ",
                            feature_table_normfile," ",
                            sample_group_file," --out_file ",pca_box_plot_out_file)
  record_command(command = pca_box_plot_str,outdir = all_out_dir)
  system(pca_box_plot_str)
}else{
  print("skip pca_box_plot.R")
}


## 绘制PCoA和NMDS
pcoa_plot_out_file <- paste0(beta_diversity_outdir,"/","beta_diversity_pcoa_plot.pdf")
nmds_plot_out_file <- paste0(beta_diversity_outdir,"/","beta_diversity_nmds_plot.pdf")

if (!file.exists(paste0(beta_diversity_outdir_log,"/beta_diversity_pcoa_nmds.complete.log"))) {
  pcoa_nmds_plot_str <- paste0("Rscript ",scripts_dir,"/beta_diversity_pcoa_nmds.R ",
                            feature_table_normfile," ",
                            sample_group_file," --out_file_pcoa ",pcoa_plot_out_file,
                            " --out_file_nmds ",nmds_plot_out_file)
  record_command(command = pcoa_nmds_plot_str,outdir = all_out_dir)
  system(pcoa_nmds_plot_str)
}else{
  print("skip beta_diversity_pcoa_nmds.R")
}


## 绘制ANOSIM箱型图
anosim_box_plot_out_file <- paste0(beta_diversity_outdir,"/","beta_diversity_anosim_plot.pdf")

if (!file.exists(paste0(beta_diversity_outdir_log,"/beta_diversity_anosim_boxplot.complete.log"))) {
  anosim_box_plot_str <- paste0("Rscript ",scripts_dir,"/beta_diversity_anosim_boxplot.R ",
                            feature_table_normfile," ",
                            sample_group_file," --out_file ",anosim_box_plot_out_file)
  record_command(command = anosim_box_plot_str,outdir = all_out_dir)
  system(anosim_box_plot_str)
}else{
  print("skip beta_diversity_anosim_boxplot.R")
}





# alpha rarefaction ！！！
# 这个命令的目的之一就是为你进行多次抽平并计算多样性，从而允许你查看不同抽平深度下多样性的变化情况。--p-min-depth和--p-max-depth参数来定义这个范围
feature_table_file <- paste0(denoise_outdir,"/feature-table_filter.qza")
meta_data_file <- meta_data_file
tree_file <- paste0(tree_outdir,"/rooted-tree.qza")
alpha_diversity_outdir <- paste0(all_out_dir,"/",alpha_diversity_dirname,"/")
# max_depth_for_rare

min_depth_for_rare

if (!file.exists(paste0(alpha_diversity_outdir,"/complete.log"))) {

  print("run alpha diversity analysis.....................")
  alpha_rarefaction_str <- paste0("bash ",scripts_dir,"/alpha_rarefaction.sh ",
                                  feature_table_file," ",meta_data_file," ",
                                  tree_file," ",min_depth_for_rare," ",alpha_diversity_outdir)
  record_command(command = alpha_rarefaction_str,outdir = all_out_dir)
  system(alpha_rarefaction_str)
}else {
  print("skip alpha_rarefaction.sh")
}



# 对alpha rarefaction结果绘制曲线图 展示抽取的reads和物种多样性的关系
sample_group_file
if (!file.exists(paste0(alpha_diversity_outdir,"/curves_plot_complete.log"))) {
  alpha_rarefaction_cruves_str <- paste0("Rscript ",scripts_dir,"/alpha_rarefaction_curves.R --alpha_rarefaction_outdir ",
                                      alpha_diversity_outdir," --groupfile ",sample_group_file," --color_by group")
  record_command(command = alpha_rarefaction_cruves_str,outdir = all_out_dir)
  system(alpha_rarefaction_cruves_str)
}else {
  print("skip alpha_rarefaction_curves.R")
}




# 对alpha rarefaction结果绘制箱线图，展示分组之间物种多样性的差异
sample_group_file
if (!file.exists(paste0(alpha_diversity_outdir,"/box_plot_complete.log"))) {
  alpha_rarefaction_boxplot_str <- paste0("Rscript ",scripts_dir,"/alpha_rarefaction_boxplot.R --alpha_rarefaction_outdir ",
                                      alpha_diversity_outdir," --groupfile ",sample_group_file)
  record_command(command = alpha_rarefaction_boxplot_str,outdir = all_out_dir)
  system(alpha_rarefaction_boxplot_str)
}else {
  print("skip alpha_rarefaction_boxplot.R")
}


# # 移动文件到plotdata文件夹
# plotdata_dir <- paste0(all_out_dir,"/plotdata/")
# if (!file.exists(paste0(all_out_dir,"/plotdata/complete.log"))) {
#   plotdata_str <- paste0("bash ",scripts_dir,"/cp_files2plotdata_dir.sh ",
#                          denoise_outdir," ",taxonomy_outdir," ",plotdata_dir)
#   print(plotdata_str)
#   system(plotdata_str)
# }else {
#   print("skip cp_files2plotdata_dir.sh")
# }



# 绘制分组的特征序列venn图  过滤0的OTU   不同水平
plot_dir <- paste0(all_out_dir,"/",plot_dirname,"/")
venn_plot_out_file <- paste0(plot_dir,"/feature_seq_venn.pdf")
plot_log_dir <- paste0(plot_dir,"/log/")
if (!dir.exists(plot_log_dir)){
dir.create(plot_log_dir,recursive = T)
}


if (!file.exists(paste0(plot_log_dir,"/venn.complete.log"))) {
  feature_seq_venn_plot_str <- paste0("Rscript ",scripts_dir,"/venn_plot.R ",
                                      feature_table_tsvfile," ",
                                      sample_group_file," -o ",venn_plot_out_file)
  record_command(command = feature_seq_venn_plot_str,outdir = all_out_dir)
  system(feature_seq_venn_plot_str)
}else {
  print("skip venn_plot.R")
}



# # 物种分类堆叠条形图
stack_bar_plot_outdir <- paste0(plot_dir,"/taxon_stack_barplot/")



if (!file.exists(paste0(stack_bar_plot_outdir,"/taxon_stack_barplot.complete.log"))) {
  stack_bar_plot_str <- paste0("Rscript ",scripts_dir,"/taxon_stack_barplot.R ",
                               feature_table_rarefile," ",
                               sample_group_file," ",
                               taxonomy_modified_file," ",
                               "--outdir ",stack_bar_plot_outdir)
  record_command(command = stack_bar_plot_str,outdir = all_out_dir)
  system(stack_bar_plot_str)
}else{
  print("skip taxon_stack_barplot.R")
}




# lefse分析
lefse_outdir <- paste0(all_out_dir,"/",lefse_dirname,"/")
lefse_inputfile <- paste0(lefse_outdir,"/lefse_inputfile.txt")

if (!file.exists(paste0(lefse_outdir,"/run_lefse_complete.log"))) {

  print("run lefse analysis.....................")
  generate_lefse_input_str <- paste0("Rscript ",scripts_dir,"/generate_lefse_input.R ",
                      feature_table_rarefile," ",
                      sample_group_file," ",
                      taxonomy_modified_file," ",
                      "--out_file ",lefse_inputfile)


  run_lefse_str <- paste0("bash ",scripts_dir,"/run_lefse.sh ",
                      lefse_inputfile," ",
                      lefse_outdir," ",
                      lda," ")


  record_command(command = generate_lefse_input_str,outdir = all_out_dir)
  record_command(command = run_lefse_str,outdir = all_out_dir)


  system(generate_lefse_input_str)
  system(run_lefse_str)
}else {
  print("skip run_lefse.sh")
}





# feature-table进行过过滤，用于PICRUST2分析
# feature_table_rarefile = "./output/analysis_result1/denoise/export/feature-table_norm.txt"
# feature_table_df <- read.table(feature_table_rarefile,header = T,row.names = 1)
# # min_sample <- ceiling(ncol(feature_table_df)*0.1)
# min_sample <- 1
# min_freq <- 10                     # 总数小于10

# if (!file.exists(paste0(denoise_outdir,"/filter_complete.log"))) {

#   feature_table_filter_str <- paste0("Rscript ",scripts_dir,"/qiime_filetr_for_picrust2.sh ",
#                                     denoise_outdir," ",
#                                     min_sample," ",
#                                     min_freq," ")
#   record_command(command = feature_table_filter_str,outdir = all_out_dir)
#   system(feature_table_filter_str)
# }else {
#    print("skip qiime_filetr_for_picrust2.sh") 

# }




# 使用picrust2进行功能预测

# 运行picrust2
feature_seq_fastafile <- paste0(denoise_outdir,"/export/feature-seq_filter.fasta")
feature_table_biomfile <- paste0(denoise_outdir,"/export/feature-table_filter.biom")
KEGG2KO_file
picrust2_out_dir <- paste0(all_out_dir,"/",picrust2_dirname,"/")

if (!file.exists(paste0(picrust2_out_dir,"/run_picrust2_complete.log"))) {

  print("run picrust2.....................")
  picrust2_str <- paste0("bash ",scripts_dir,"/run_picrust2.sh ",
                         feature_seq_fastafile," ",
                         feature_table_biomfile," ",
                         KEGG2KO_file," ",
                         picrust2_out_dir)
  record_command(command = picrust2_str,outdir = all_out_dir)
  system(picrust2_str)
}else {
  print("skip run_picrust2.sh")
}



# KEGG通路差异分析

kegg_abu_file <- paste0(picrust2_out_dir,"/KEGG_pathways_out/path_abun_unstrat_descrip.tsv")
sample_group_file
logfc
use_adjp
kegg_diff_analysis_outdir <- paste0(picrust2_out_dir,"/kegg_diff_analysis/")

if (!file.exists(paste0(kegg_diff_analysis_outdir,"/log/pathway_diff_analysis.log"))) {

  print("run kegg_diff_analysis .....................")
 pathway_diff_analysis_str <- paste0("Rscript ",scripts_dir,"/pathway_diff_analysis.R ",
                         kegg_abu_file," ",
                         sample_group_file," ",
                         "--logfc ",logfc," ",
                         "--use_adjp ",use_adjp," ",
                         "--out_path ",kegg_diff_analysis_outdir)


  record_command(command = pathway_diff_analysis_str,outdir = all_out_dir)
  system(pathway_diff_analysis_str)
}else {
  print("skip pathway_diff_analysis.R")
}




#  差异KEGG通路热图
library(forcats)
group_df <- read.table(sample_group_file,header = T,row.names = 1)
group_df$group <- fct_inorder(group_df$group)
group_level <- group_df$group %>% levels()
ref_group <- group_level[1]
case_group <- group_level[2]


# 参数
kegg_abu_file
sample_group_file
diff_analysis_res_file <- paste0(kegg_diff_analysis_outdir,"/DESeq2_deg_group_",case_group,"_vs_",ref_group,".csv")


if (!file.exists(paste0(kegg_diff_analysis_outdir,"/log/diff_pathway_heatmap.log"))) {

 diff_pathway_heatmap_str <- paste0("Rscript ",scripts_dir,"/diff_pathway_heatmap.R ",
                         kegg_abu_file," ",
                         sample_group_file," ",
                         diff_analysis_res_file," ",
                         "--out_path ",kegg_diff_analysis_outdir)


  record_command(command = diff_pathway_heatmap_str,outdir = all_out_dir)
  system(diff_pathway_heatmap_str)
}else {
  print("skip  diff_pathway_heatmap.R")
}





#  差异KEGG误差棒图


# 参数
kegg_abu_file
sample_group_file
diff_analysis_res_file <- paste0(kegg_diff_analysis_outdir,"/DESeq2_deg_group_",case_group,"_vs_",ref_group,".csv")
use_adjp


if (!file.exists(paste0(kegg_diff_analysis_outdir,"/log/diff_pathway_errorbar.log"))) {

 diff_pathway_errorbar_str <- paste0("Rscript ",scripts_dir,"/diff_pathway_errorbar.R ",
                         kegg_abu_file," ",
                         sample_group_file," ",
                         diff_analysis_res_file," ",
                         "--use_adjp ",use_adjp," ",
                         "--out_path ",kegg_diff_analysis_outdir
                         )


  record_command(command = diff_pathway_errorbar_str,outdir = all_out_dir)
  system(diff_pathway_errorbar_str)
}else {
  print("skip diff_pathway_errorbar.R")
}












                              




# 移动文件
#=====================================================================================
# all_out_dir = "./analysis_result_v3/"

result_clean <- paste0(all_out_dir,"/result_clean/")
dir.create(path = result_clean,recursive = T)






# dada2结果
# denoise_outdir = "./analysis_result_v3/denoise/"
dada2_result <- paste0(result_clean,"/01_dada2_result/")

if(!dir.exists(dada2_result)){
  dir.create(dada2_result,recursive = T)
}
if(!dir.exists(paste0(dada2_result,"/export/"))){
  dir.create(paste0(dada2_result,"/export/summary_feature_table/"),recursive = T)
}

fs <- list.files(paste0(denoise_outdir),recursive = T)
for (f in fs) {
  file.copy(from = paste0(denoise_outdir,"/",f),to = paste0(dada2_result,f),overwrite = T)
}



# 物种注释结果
# taxonomy_outdir = "./analysis_result_v3/taxonomy/"
taxonomy_result <- paste0(result_clean,"/02_taxonomy_result/")

if(!dir.exists(paste0(taxonomy_result,"/export/"))){
    # dir.create(taxonomy_result,recursive = T)
    dir.create(paste0(taxonomy_result,"/export/"),recursive = T)
}

fs <- list.files(paste0(taxonomy_outdir),recursive = T)
for (f in fs) {
  file.copy(from = paste0(taxonomy_outdir,"/",f),to = paste0(taxonomy_result,f),overwrite = T)
}



# 系统发育树
# tree_outdir = "./analysis_result_v3/tree/"
phylogenetic_tree_result <- paste0(result_clean,"/03_phylogenetic_tree_result/")
if(!dir.exists(phylogenetic_tree_result)){
  dir.create(phylogenetic_tree_result,recursive = T)
  dir.create(paste0(phylogenetic_tree_result,"/export/"),recursive = T)
}

fs <- list.files(paste0(tree_outdir),recursive = T)
for (f in fs) {
  file.copy(from = paste0(tree_outdir,"/",f),to = paste0(taxonomy_result,f),overwrite = T)
}



# alpha多样性分析结果
# alpha_diversity_outdir = "./analysis_result_v3/alpha_diversity/"


alpha_diversity_result <- paste0(result_clean,"/04_alpha_diversity_result/")

if(!dir.exists(paste0(alpha_diversity_result,"/export/"))){
    # dir.create(alpha_diversity_result,recursive = T)
    dir.create(paste0(alpha_diversity_result,"/export/"),recursive = T)
}

fs <- list.files(paste0(alpha_diversity_outdir),recursive = T,pattern = ".csv|.pdf|.log|.qzv")
for (f in fs) {
  file.copy(from = paste0(alpha_diversity_outdir,"/",f),to = paste0(alpha_diversity_result,f),overwrite = T)
}




# beta多样性分析结果
# beta_diversity_outdir = "./analysis_result_v3/beta_diversity/"

beta_diversity_result <- paste0(result_clean,"/05_beta_diversity_result/")
if(!dir.exists(beta_diversity_result)){
    dir.create(paste0(beta_diversity_result,"/log/"),recursive = T)
}


fs <- list.files(paste0(beta_diversity_outdir),recursive = T)
for (f in fs) {
  file.copy(from = paste0(beta_diversity_outdir,"/",f),to = paste0(beta_diversity_result,f),overwrite = T)
}



# lefse分析结果
# lefse_outdir = "./analysis_result_v3/lefse/"

lefse_result <- paste0(result_clean,"/06_lefse_result/")

if(!dir.exists(paste0(lefse_result))){
    dir.create(paste0(lefse_result),recursive = T)
}

fs <- list.files(paste0(lefse_outdir),recursive = T)
for (f in fs) {
  file.copy(from = paste0(lefse_outdir,"/",f),to = paste0(lefse_result,f),overwrite = T)
}






# picrust2_out_dir = "./output/analysis_result_0927_v2/picrust2/"
library(dplyr)
library(stringr)
picrust2_result <- paste0(result_clean,"/07_picrust2_result/")


if(!dir.exists(paste0(picrust2_result,"/kegg_diff_analysis/"))){
  dir.create(paste0(picrust2_result,"/kegg_diff_analysis/"),recursive = T)
}

fs <- list.files(paste0(picrust2_out_dir),recursive = T,pattern = "unstrat_descrip.tsv|.pdf|.csv")
fs_new <- fs %>% str_replace(pattern = "_metagenome_out/pred_metagenome_",replacement = "_abun_")
fs_new <- fs_new %>% str_remove(pattern = "_out/path")

for (i in 1:length(fs)) {
  file.copy(from = paste0(picrust2_out_dir,"/",fs[i]),to = paste0(picrust2_result,fs_new[i]),overwrite = T)
}












# 其他分析结果
# plot_dir = "./analysis_result_v3/plot/"
other_result <- paste0(result_clean,"/other_result/")

if(dir.exists(paste0(other_result))){
  # 存在就删除
  fs::dir_delete(paste0(other_result))
}

if(!dir.exists(paste0(other_result))){
    dir.create(paste0(other_result),recursive = T)
    fs::dir_copy(path = plot_dir,new_path = other_result,overwrite = T)
}





#=====================================================================================






rmarkdown::render(paste0(scripts_dir,'/16s_analysis_report.qmd'), output_file = paste0(result_clean,"/16s_analysis_report"))




