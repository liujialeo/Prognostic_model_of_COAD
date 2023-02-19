################################################
################################################
### 作者：果子
### 更新时间：2021-01-18
### 微信公众号:果子学生信
### 私人微信：guotosky
### 个人博客: https://codingsoeasy.com/

rm(list = ls())
# 02_读入数据
### 原始文件在rawdata
### 首先把所有数据读入在一个文件夹中
### 创建目的文件夹，可以右击创建，也可以dir.create函数创建
dir.create("data_in_one")
dir("rawdata/")
### 使用for循环来批量做，回顾for循环的要点
for (dirname in dir("rawdata/")){  
  ## 要查看的单个文件夹的绝对路径
  mydir <- paste0(getwd(),"/rawdata/",dirname)
  ##找到对应文件夹中的文件并提取名称，pattern表示模式，可以是正则表达式
  file <- list.files(mydir,pattern = "*.counts")
  ## 当前文件的绝对路径是
  myfile <- paste0(mydir,"/",file)
  #复制这个文件到目的文件夹
  file.copy(myfile,"data_in_one")  
}

### 在终端也可以简单实现这个效果：cp rawdata/*/*.gz data_in_one

## 文件和TCGA ID的对应关系
metadata <- jsonlite::fromJSON("metadata.cart.2022-03-03.json")
library(dplyr)
metadata_id <- metadata %>% 
  dplyr::select(c(file_name,associated_entities)) 

## 提取对应的文件
naid_df <- data.frame()
for (i in 1:nrow(metadata)){
  naid_df[i,1] <- metadata_id$file_name[i]
  naid_df[i,2] <- metadata_id$associated_entities[i][[1]]$entity_submitter_id
}

colnames(naid_df) <- c("filename","TCGA_id")

################################################################################
## 读入数据
### 1.存储文件读入的顺序
files <- dir("data_in_one")

## 构建函数
myfread <- function(files){
  data.table::fread(paste0("data_in_one/",files))$V2
}

## 测试文件以及时间，大概2秒一个
system.time(test <- myfread(files[1]))

###############################################
## lapply 的使用方法
## 1.for 循环正常速度是1222*2 大概40分钟，所以不适合
## 2.此处使用的是lapply+ function的形式，目的是为了提速
## 3.提速的实现方式是，lapply的并行化处理，使用的是future.apply这个包

library(future.apply)
plan(multisession)
### 用system.time返回时间
### 我的台式机12个线程大概需要440秒
system.time(f <- future_lapply(files,myfread))

## 列表变成数据框，回想晾衣杆的结构
expr_df <- as.data.frame(do.call(cbind,f))

## 为了把文件名称和TCGAid对应起来
rownames(naid_df) <- naid_df[,1]
naid_df <- naid_df[files,]

## 命名
colnames(expr_df) <- naid_df$TCGA_id

test <- expr_df[1:100,1:4]
## 加上一列
gene_id <- data.table::fread(paste0("data_in_one/",files[1]))$V1

expr_df <- cbind(gene_id=gene_id,expr_df)
test <- expr_df[1:100,1:4]

### 意外发现
tail(expr_df$gene_id,10)
expr_df <- expr_df[1:(nrow(expr_df)-5),]

dir.create("output")
save(expr_df,file = "./output/COAD_RNASEQ_exprdf.Rdata")

###############################################################################
##读入临床信息
rm(list = ls())
clinical = data.table::fread(file = "./TCGA-COAD.GDC_phenotype.tsv")
clinical <- as.data.frame(clinical)
tmp <- as.data.frame(colnames(clinical))
### 把要提取的条目变成向量，使用的是c()函数
index <- c( "submitter_id.samples",
            "vital_status.demographic",
            "days_to_death.demographic",
            "days_to_last_follow_up.diagnoses",
            "gender.demographic",
            "age_at_initial_pathologic_diagnosis",
            "tumor_stage.diagnoses",
            "pathologic_T",
            "pathologic_M",
            "pathologic_N",
            "history_of_neoadjuvant_treatment"
            )
clinical <- clinical[,index]
## 修改列名
colnames(clinical) <- c("TCGA_id",
                                  "vital_status",
                                  "days_to_death",
                                  "days_to_last_followup",
                                  "gender",
                                  "age",
                                  "stage",
                                  "T",
                                  "M",
                                  "N",
                                  "history_of_neoadjuvant_treatment")
## 增加两列,都是NA
clinical$fustat <- NA
clinical$futime <- NA
## NA 变成0，使用is.na 判断得到逻辑，是T的变成0
clinical[is.na(clinical)] =0

## 数据整理
## for 循环解决问题
for (i in 1:nrow(clinical)) {
  if(clinical$vital_status[i]=="Alive"){
    clinical$fustat[i]= 0
    clinical$futime[i]= clinical$days_to_last_followup[i]
  }else{
    clinical$fustat[i]= 1
    clinical$futime[i]=clinical$days_to_death[i]
  }
}

## 提取数据并且调整顺序
clinical <- clinical[,c("TCGA_id","fustat", "futime","gender", "age", "stage","T","M","N")]
## TCGA barcode 小写变成大写
clinical$TCGA_id  <- toupper(clinical$TCGA_id)

## 保存数据
save(clinical,file = "output/COAD_clinical.Rdata")
