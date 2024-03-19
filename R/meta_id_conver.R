
#-------------------------------------------------------------

LXmetaid <- function(gene_data,meta_data,species){

#---安装R包----------------
inst_packages <- function(){
  
installed_packs <- installed.packages()[,1] 
com_packs <- c("openxlsx","dbplyr","devtools")
bio_packs <- c("clusterProfiler")
git_packs <- c("CTSgetR")

not_com <- com_packs[!com_packs %in% installed_packs]
if(length(not_com)>0){
  com_fun <- function(i){install.packages(i)}
  sapply(com_packs,com_fun)
  }

not_bio <- bio_packs[!bio_packs %in% installed_packs]
if(length(not_bio)>0){
  bio_fun <- function(i){BiocManager::install(i,update = F,ask = F)}
  sapply(bio_packs,bio_fun)
  }

not_git <- git_packs[!git_packs %in% installed_packs]
if(length(not_git)>0){
  library(devtools)
  install_github("dgrapov/CTSgetR")
  }

lib_fun <- function(i){library(i,character.only = T)}
sapply(c(com_packs,bio_packs,git_packs),lib_fun)
}

inst_packages()

# https://github.com/dgrapov/CTSgetR
# library(devtools)
# install_github("dgrapov/CTSgetR")

# -----------CTSgetR基本用法-------------------------
# CTSgetR(id, from, to, db_name = NULL, ...)
# id：需要转换的代谢物ID列表
# from：输入的代谢物ID的标识符类型
# to：需要转换成的代谢物标识符类型
# db_name：sqlite数据库存储缓存结果的字符串路径（可加快转换速度）

# --------------确认CTS API可用----------------------
# GET('https://cts.fiehnlab.ucdavis.edu/services') %>%
#  http_status(.) %>%
#  {if( .$category != 'Success'){stop('Oops looks like https://cts.fiehnlab.ucdavis.edu/services is down!') }} 

# -----------查看可用数据库（展示部分）--------------
# head(unlist(valid_from()))
# head(unlist(valid_to()))
# from_df=unlist(valid_from())
# to_df= unlist(valid_to())


#------------ 查找名称-----------------------------
# want<-'KEGG'
# from_df[grepl(want,from_df,ignore.case=TRUE)]
# to_df[grepl(want,to_df,ignore.case=TRUE)]

# want<-"Human Metabolome Database"
# from_df[grepl(want,from_df,ignore.case=TRUE)]
# to_df[grepl(want,to_df,ignore.case=TRUE)]

# want<-"compound"
# from_df[grepl(want,from_df,ignore.case=TRUE)]
# to_df[grepl(want,to_df,ignore.case=TRUE)]

# ---------示例：----------------------------------
# CTSgetR(id = c("C00026","C05381"),from = "KEGG",
#        to = c("Human Metabolome Database","PubChem CID"))

# CTSgetR(id = c("alanine",'lactic acid'),from = "Chemical Name",
#        to = c("KEGG","Human Metabolome Database","PubChem CID"))


# hmdb_list <- read.xlsx("meta_list.xlsx")
# conv_id <- CTSgetR(id = hmdb_list[,1],from = "Human Metabolome Database",to = c("KEGG","PubChem CID"))

#-----建立文件夹-------------------------------------------------------------

dir.name <- paste0("analysi_resut_",trimws(species))
if(!dir.exists(dir.name))
 dir.create(dir.name)

#------ gene 数据整理---------------------------------------------------------
gen_set0 <- read.xlsx(gene_data)
colnames(gen_set0)[c(1,3)] <- c("gene_id","pvalue") # 第一、三列重命名
gen_set0[,1] <- trimws(gen_set0[,1]) # 去掉前后空格
gen_set0 <- distinct(.data =gen_set0,gene_id,.keep_all = T ) # 对gene_id去重

#---找出组名---
group_names <- colnames(gen_set0) %>% 
  gsub("\\d+$","",.) %>%   # 去掉字符末尾的数字，等同 gsub("\\d+$","", group_gene)
  .[duplicated(.)] %>% # 选出重复项（只有组名才可能重复，因此重复即为组名）
  .[!duplicated(.)] # 去掉重复的字符串，等同 group_gene[!duplicated(group_gene)] 

group_1 <- group_names[1]  # 导出组别名称
group_2 <- group_names[2]

#---默认前面是对照组(control)，后面是干预组(treat)---
g1_col <- grep(group_1,colnames(gen_set0)) # group_1所在的列
g2_col <- grep(group_2,colnames(gen_set0)) # group_2所在的列

for(i in 1:nrow(gen_set0)){
  if(rowMeans(gen_set0[i,g2_col])/rowMeans(gen_set0[i,g1_col])<1) 
    gen_set0[i,2] <-0- abs(gen_set0[i,2]) else
      gen_set0[i,2] <-abs(gen_set0[i,2])
}

#---lgo FC 转换为 FC ---
if( any(grepl("log2", colnames(gen_set0),ignore.case = T)) )
{log2_col <- grep("log2", colnames(gen_set0)); gen_set0$FC <- 2^gen_set0[,log2_col]} else
{ if( any(grepl("log10", colnames(gen_set0),ignore.case = T)) )
{log10_col <- grep("log10", colnames(gen_set0)); gen_set$FC <- 10^gen_set0[,log10_col]} else
  gen_set$FC <- gen_set[,2]
}

#-------背景genes数据-------
#----kegg 数据库的选择---
if(species=="human")
{library(org.Hs.eg.db); Org.Db= org.Hs.eg.db} else
{ if(species=="rat")
{library(org.Rn.eg.db); Org.Db= org.Rn.eg.db} else
{library(org.Mm.eg.db); Org.Db= org.Mm.eg.db}
}

entr_id_back <- bitr(geneID = gen_set0$gene_id,fromType = "SYMBOL",
                toType = c("SYMBOL","ENTREZID"),OrgDb =Org.Db) %>% data.frame()

colnames(entr_id_back)[1] <- c("gene_id")

gen_back_lst_all <- dplyr::left_join(gen_set0,entr_id_back,by = "gene_id")
gen_back_lst <- dplyr::inner_join(gen_set0,entr_id_back,by = "gene_id") %>% na.omit()
gen_back_lst <- gen_back_lst[,c("ENTREZID","pvalue","FC")]

write.csv(gen_back_lst_all,paste0(dir.name,"/gen_back_lst_all.csv"),row.names = F)
write.table(gen_back_lst,paste0(dir.name,"/gen_back_lst (for LRpath).txt"),row.names = F,quote = F,sep ="\t")

#------差异geens数年---------------------------------------------------------
gen_set <- filter(gen_set0,pvalue<0.05) # 筛选出p<0.05的差异表达genes

#---- gene symbol 转换成 ENTREZ
entr <- bitr(geneID = gen_set$gene_id,fromType = "SYMBOL",toType = c("SYMBOL","ENTREZID"),OrgDb =Org.Db) %>% data.frame()
colnames(entr)[1] <- "gene_id"

gen_df_all <- dplyr::left_join(gen_set,entr,by = "gene_id")
gen_df <- dplyr::inner_join(gen_set,entr,by = "gene_id")

entr_col <- grep("ENTREZ",colnames(gen_df),ignore.case = T)
p_col <- grep("p",colnames(gen_df),ignore.case = T)
FC_col <- grep("FC",colnames(gen_df),ignore.case = T)

gene_lst <- gen_df[,c(entr_col,p_col,FC_col)]
colnames(gene_lst) <- c("ENTREZID","p-value","Fold Change")

write.csv(gen_df_all,paste0(dir.name,"/gene_lst_all.csv"),row.names = F)
write.csv(gene_lst,paste0(dir.name,"/gene_lst (for MetScape analysis).csv"),row.names = F)

#---------gene数据处理完成----------------------------------------------------

#-------metabolite数据处理-----------------------
meta_set <- read.xlsx(meta_data)
colnames(meta_set) <- c("HMDB_ID","pvalue","FC","highmean") # 重命名
meta_set <- distinct(.data = meta_set,HMDB_ID,.keep_all = T) # 对meta_id去重
meta_set$FC <- as.numeric(meta_set$FC)
meta_set <- filter(meta_set,pvalue<0.05) # 筛选出p<0.05

group_meta <- meta_set$highmean %>% .[!duplicated(.)] # 找出组名，并去重
g_m_1 <- group_meta[1] # 第一组名字
g_m_2 <- group_meta[2]

is.exist.group <- c(any(grepl(g_m_1,group_names,ignore.case = T)), # 判断meta和gene的组名是否一样
                   any(grepl(g_m_2,group_names,ignore.case = T)) 
                   )
                   
if(!all(is.exist.group)) # 如果meta和gene的组名不一样，则停止运行
  stop("The groups in the gene_data and meta_data are different. Please check them.")

#-----处理FC列--------
high_col <- grep("high",colnames(meta_set))
FC_col <- grep("FC",colnames(meta_set))

for(i in 1:nrow(meta_set))
   {if(grepl(meta_set[i,high_col],group_2,ignore.case = T))
      meta_set[i,FC_col] <- meta_set[i,FC_col] else
        meta_set[i,FC_col] <- 1/meta_set[i,FC_col]
      }

#----HMDB ID转换成KEGG ID-------
met_id <- CTSgetR(id = meta_set[,1],from = "Human Metabolome Database",to = c("KEGG","PubChem CID")) %>% data.frame()
colnames(met_id) <-c("HMDB_ID","KEGG_ID","PubChem_CID") 

met_df_all <- dplyr::left_join(meta_set,met_id,by="HMDB_ID")
met_kegg <- met_df_all[,c("KEGG_ID","pvalue","FC")] %>% na.omit()
met_PubChem <- met_df_all[,c("PubChem_CID","pvalue","FC")] %>% na.omit()

write.csv(met_df_all,paste0(dir.name,"/metabolite_lst_all.csv"),row.names = F)
write.csv(met_kegg,paste0(dir.name,"/metabolite_lst (for MetScape analysis).csv"),row.names = F)
#write.xlsx(met_PubChem,paste0(dir.name,"/meta_PubChem_lst (for MetScape analysis).xlsx"),colNames=T)

#-------Metabolites 处理结束---------------------------------------------------

print(paste0("The data for MetScape analsyis can seen in the folder of <",dir.name,">"))

}

