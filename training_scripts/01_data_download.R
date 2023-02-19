################################################
################################################
### 作者：果子
### 更新时间：2020-01-04
### 微信公众号:果子学生信
### 私人微信：guotosky
### 个人博客: https://codingsoeasy.com/

### 1.从GDC下载counts数据
## https://portal.gdc.cancer.gov/repository
dir.create("rawdata")
manifest <- "gdc_manifest_20220303.txt"
rawdata <- "rawdata"
command <- sprintf("./gdc-client download -m %s -d %s",
                   manifest,
                   rawdata)
system(command = command)

### 备选方案
### 2.下载数据，使用TCGA GDC以及https://xena.ucsc.edu/
### 注意事项，下载的数据是log之后的数据，所以不能直接使用Deseq2来做差异分析。


###从GDC下载clinical数据
# dir.create("clinical")
# manifest <- "gdc_manifest.clinical.txt"
# clinical <- "clinical"
# command <- sprintf("./gdc-client download -m %s -d %s",
#                    manifest,
#                    clinical)
# system(command = command)
