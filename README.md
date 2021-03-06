# DADA2-and-phyloseq
#Hướng dẫn sử dụng DADA2
Link của nhà phát triển package dada2 https://benjjneb.github.io/dada2/tutorial.html


Để xử lý với những dữ liệu sinh học, R-Project là một sự lựa chọn khá tất yếu khi R mang lại một môi trường tính toán thống kê chính xác và linh hoạt. Trong đó Bioconductor là một kho các package của R chuyên để sử lí các số liệu sinh học đặc biệt là dữ liệu trình tự nucleotide. 
Trong nghiên cứu này package Dada2 và phyloseq được áp dụng:

###Mục đích của pipeline
Trong bài hướng dẫn này, OSD-2014 dataset được sử lý qua pipeline DADA2. Bắt đầu từ các trình tự đã được đọc bằng công nghệ đọc trình tự của Ilumina lưu dưới dạng file fastq. Sau khi sử lý qua bằng DADA2, các trình tự của từng mẫu trên được tổng hợp vào bảng ASV - amplicon sequence variant, bảng ASV chứa thông tin về số lượng OTU, và số lần OTU đó có mặt trong mẫu phân tích.

Đồng thời trong bài hướng dẫn này cũng trình bày cách đưa các trình tự này vào pipeline phyloseq và phân tícn microbiome data.

Sau khi cài các package bằng function bioClite(), áp dụng chúng vào R session bằng function library()

###Cài đặt các packages
```{r}
source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite("dada2")
biocLite('phyloseq')
biocLite("DECIPHER")
biocLite("phangorn")
biocLite("BiocStyle")
biocLite("sequencing")
biocLite("BiocGenerics")
biocLite("msa")
```
###Load các package cấn thiết khi khởi động một section
```{r}
rm(list=ls())
library("dada2"); packageVersion("dada2")
library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library("sequencing")
library("BiocGenerics")
library("msa")
```


Dữ liệu được sử dụng trong ví dụ này được cung cấp bởi Jean-Christophe, Pháp. Dữ liệu trên được thu thập từ 2 lần thu mẫu tại các biển (open ocean) và các vùng bờ biển trên thế giới. (http://ocean-microbiome.embl.de/companion.html)
Trong thí nghiệm, sự đa dạng sinh học của các vi khuẩn trong mẫu được đánh giá dựa trên trình tự nucleotide 16S rRNA. Các đoạn rRNA này được đọc trình tự dựa trên công nghệ 2x250 của ilunima. Do vậy các đoạn trình tự thu được sẽ có độ dài khoảng 250 nucleotide.


###Tải dữ liệu và dẫn đường dẫn tới tệp chứa các file dữ liệu. Và kiểm tra thành phần các files:
```{r}
path <- "/Users/nqhuynh/Documents/rdna"
list.files(path)
```


##Lọc và cắt bỏ các đoạn trình tự có chất lượng thấp
####Chuẩn bị
   Đầu tiên, chúng ta cần đọc tên các file fastq và dùng một số hàm thay đổi chuỗi kí tự để tạo 2 danh sách file trình tự F và R.
   Các file chứa trình tự sẽ có tên dạng Tên_R1.fastq.gz hoặc Tên_R2.fastq.gz trong đó R1 là các file trình tự Forward, R2 là các file trình tự Reverse.2 loại trình tự trên được nhóm vào 2 nhóm khác nhau fnFs và fnRs.
```{r}
fnFs <- sort(list.files(path, pattern="_R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq.gz", full.names = TRUE))
```

Đưa tên mẫu về dạng: OSDXXX_.fastq
```{r}
sample.names <- make.unique(sapply(strsplit(basename(fnFs), "_"), `[`, 1))
head(sample.names)
sample.namesfull <- sapply(strsplit(basename(fnFs), "_16S"), `[`, 1)
```

###Kiểm tra chất lượng trình tự
```{r fig.height=22, fig.width=20}
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])
```
Thang màu xám là bản đồ nhiệt (heat-map) thể hiện chất lượng của trình tự của từng nucleotide. Đường màu xanh biểu diễn điểm chất lượng trung bình của mỗi vị trí nucleotide, và đường màu cam biểu diễn tứ phân vị sự phân bố điểm chất lượng các nucleotide. Đường màu đỏ biểu diễn tỷ lệ số trình tự có độ dài ở vị trí nucleotide đó. Đối với ví dụ này trình tự được đọc trên máy Ilumina các đoạn đọc có độ dài như nhau ~250Nu, vì vậy mà đường màu đỏ ở ví dụ này nằm song song với trục x.

###Cắt và lọc trình tự
Gán tên file cho các trình tự sau khi cắt và lọc các file fastq.gz.
```{r}
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
```

Cắt và lọc các trình tự, tuỳ vào chất lượng của trình tự mà các thông số dưới đây có thể được thay đổi, thông thường điểm chất lượng trung bình khoảng >20 là có thể chấp nhận.
Thông số truncLen = c(F, R) dùng để cắt bỏ đoạn đuôi trình tự có chất lượng thấp, như thể hiện trên hình chất lượng trình tự, đoạn cuối hay đuôi của trình tự thường có chất lượng kém, nên tuy đoạn trình tự có độ dài là 250Nu, nhưng chỉ 220Nu được dùng cho các bước tiếp theo.
Thông số maxEE = c(F, R) dùng để khai báo sai số cho phép đi qua bước lọc. Thông số này càng nhỏ (2) thì sai số cho phép đi qua bước lọc càng nhỏ và ngược lại, thông số maxEE càng lớn (5) sai đó cho phép đi qua bước lọc càng cao. Điều này có nghĩa là khi maxEE = c(2,2) sai số của trình tự được phép đi qua bước lọc nhỏ và những trình tự có sai số lớn sẽ không thể đi qua, ngược lại khi maxEE được nới lỏng c(5,5) hầu hết các trình tự đều có thể đi qua kể cả các trình tự có sai số lớn.

```{r echo=TRUE}
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(220, 220),
                     maxN=0, truncQ=2, maxEE = c(2,2), rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) 
```
Kiển tra out, reads.in là số trình tự trước khi lọc, và reads.out là số trình tự sau khi lọc
```{r}
head(out)
```


###Tỷ lệ sai của các trình tự


```{r echo=TRUE}
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
```

```{r}
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
```

###Dereplication
Dereplicate the filtered fastq files
```{r include=FALSE}
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
```


Gán tên mẫu cho derep Forward và reverse 
```{r}
names(derepFs) <- sample.names
names(derepRs) <- sample.names
```

Sample Inference
Suy ra các biến thể trình tự trong mỗi mẫu
```{r}
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
```

Kiểm tra class của dada object
```{r}
dadaFs[[1]]
dadaRs[[1]]
```


###Ghép các trình tự đã được loại bỏ nhiễu
```{r}
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
```

Kiểm tra giá trị của trình tự đầu tiên trong mergers
```{r}
head(mergers[[1]])
```

###Xây dựng bảng trình tự
```{r}
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
```

Kiểm tra sự khác nhau của độ dài các trình tự
```{r}
table(nchar(getSequences(seqtab)))
```

###Loại bỏ Chimera
```{r}
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab) 
```
Nếu giá trị  = 1, 0% trình tự chimera bị loại bỏ, 0,9=10% trình tự chimera bị loại bỏ...

###Kiểm tra các trình tự đã được sử lí qua pipline dada2
```{r}
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)
```



###Gán taxonomy
Tải "silva_nr_v132_train_set.fa" và " silva_species_assignment_v132.fa" từ "https://benjjneb.github.io/dada2/training.html"

```{r}
taxaRC <- assignTaxonomy(seqtab.nochim, "/Users/nqhuynh/Documents/Trainset/silva_nr_v132_train_set.fa", tryRC=TRUE)
taxaSp <- addSpecies(taxaRC, "/Users/nqhuynh/Documents/Trainset/silva_species_assignment_v132.fa")
```

```{r}
taxa.print <- taxaSp
rownames(taxa.print) <- NULL
head(taxa.print)
```


###Thông tin về môi trường địa điểm thời gian thu mẫu năm 2014
xem hướng dẫn tại https://github.com/MicroB3-IS/osd-analysis/wiki/OSD-2014-environmental-data-csv-documentation# Import OSD metdata

```
# The data set itself is documented here:
# https://github.com/MicroB3-IS/osd-analysis/wiki/OSD-2014-environmental-data-csv-documentation
# If you spot any errors in this data, plese post an issue to: 
# https://github.com/MicroB3-IS/osd-analysis/issues
# First step, pull the data from the database export online.
# Note: if not using Windows, you may need to use RCurl.
setInternet2(TRUE) # Allows https access
# pull CSV file, declare that its first row is a row of headers (col names)
# that it is pipe, "|", separated, and that 
# data are quoted with double quotes
# This could be (should be?) switched to PANGAEA's tsv link...
rawContextualData <- read.csv(
	"https://owncloud.mpi-bremen.de/index.php/s/eaB3ChiDG6i9C6M/download?path=%2F2014%2Fsamples&files=OSD2014_environmental_metadata_2015-11-17.csv",
	header = T,
	sep = "|",
	quote = '"',
	colClasses = "character"
	)
# Set row names to a (sensible) unique identifier, here I chose the ENA 
# accession code, already in the data.
row.names(rawContextualData) <- rawContextualData[, "label"]
# check the dimensions of the table...
#  dim(rawContextualData)
#    [1] 203  54
# Now, partition the table into objects that contain logically aggregated
# Create a table mapping various IDs to one another...
sampleMap <- rawContextualData[
	,
	c(
		"osd_id",
		"label",
		"bioarchive_code",
		"ena_acc",
		"biosample_acc",
		"site_name"
		)
	]
	
# create a spatial data object
spatialData <- rawContextualData[, 
	c(
		"start_lat",
		"start_lon",
		"stop_lat",
		"stop_lon",
		"water_depth"
	)
]
# cast as numeric
spatialData <- sapply(spatialData, as.numeric)
# regenerate row names
row.names(spatialData) <- row.names(rawContextualData)
# create a temporal data object
# Standardise UTC times to POSIX
# R can do calculations with these, unlike simple strings
#
# for example:
# temporalData.computable$startTimes[1:4] - temporalData.computable$startTimes[5:8]
# Time differences in days
# [1] -1.239583 -3.294410 -2.386806 -2.303472
#
# row.names are lost, however, but the order is preserved
temporalData.computable <- list()
temporalData.computable$startTimes <- strptime(rawContextualData$start_date_time_utc, format="%Y-%m-%dT%H:%M:%S+00", tz = "UTC")
temporalData.computable$endTimes <- strptime(rawContextualData$end_date_time_utc, format="%Y-%m-%dT%H:%M:%S+00", tz = "UTC")
# store the textual representations of time and date.
# not useful for computation (use the list object above for that),
# but good for reference.
temporalData <- rawContextualData[, 
	c(
		"local_date",  # local at site!! Not Zulu.
		"local_start", # local at site!! Not Zulu.
		"local_end",   # local at site!! Not Zulu.
		"start_date_time_utc",
		"end_date_time_utc"
	)
]
# Create a table containing all categorical variables...
factorData <- rawContextualData[, 
	c(
		"iho_label",
		"mrgid", # Marine regions general identifer (corresponds to IHO label
		"protocol",   # 
		"biome",
		"biome_id",
		"feature",
		"feature_id",
		"material",
		"material_id"
	)
]
# Transform [biome,feature,material]_id values to valid IDs
idFields <- c("biome_id", "feature_id", "material_id")
for (i in idFields) {
  factorData[,paste(i)] <- gsub("^","http://purl.obolibrary.org/obo/ENVO_", factorData[, paste(i)])
}
rm(idFields)
# Tell R that these are "factor" data, so that they are treated like
# categorical variables...
tempNames <- names(factorData)
factorData[,tempNames] <- lapply(factorData[,tempNames], factor)
# Create a table with unstructured data (like free text descriptions)...
unstructuredFields <- rawContextualData[, 
	c(
		"objective",
		"platform", 
		"device",    
		"description"
	)
]
# create a numeric data object - some variables have only a few
# (or sometimes only one) non-empty entry!
numericData <- rawContextualData[, 
	c(
		"water_temperature",
		"salinity", 
		"ph",    
		"phosphate",
		"carbon_organic_particulate",
		"carbon_organic_dissolved_doc",
		"nitrate",
		"nitrite",
		"nano_microplankton",
		"downward_par",
		"conductivity",
		"primary_production_isotope_uptake",
		"primary_production_oxygen",
		"dissolved_oxygen_concentration",
		"nitrogen_organic_particulate_pon",
		"meso_macroplankton",
		"bacterial_production_isotope_uptake",
		"nitrogen_organic_dissolved_don",
		"ammonium",
		"silicate",
		"bacterial_production_respiration",
		"turbidity",
		"fluorescence",
		"pigment_concentration",
		"picoplankton_flow_cytometry"
	)
]
# Convert into a matrix (this is appropriate [faster and less likely to
# lead to errors] than a data frame object.
numericData <- data.matrix(numericData)
# For semantic consistency, replace NaNs with NAs
numericData[is.nan(numericData)] <- NA
combinedEnvSpatTempData <- data.frame(
  spatialData,
  temporalData,
  factorData, 
  numericData
  )
#Clean up ...
rm(list = c("tempNames", "i"))
gc()
####Kết thúc OSD 2014 environmental table###########
```
###Sau khi chạy script trong link trên lưu bảng data dưới dạng .txt 

```
write.table(combinedEnvSpatTempData, file = "/file directory/data1.txt",header=TRUE, sep="\t", dec=",")
```

Đưa bảng data1 trở lại R session
```{r}
envdata <- read.table("/Users/nqhuynh/Documents/rdna/3/data1.txt", header=TRUE, fill=TRUE, row.names=1, sep="\t", dec=",")
colnames(envdata)
```

###Tạo bảng tên đối chứng Name_list
Do tên ban đầu của dữ liệu đã seq (OSDxxx_năm-tháng-số thứ tự_16S_amala_R1.fastq) và bảng thông tin môi trường (OSDxxx_năm_tháng_ngày_độ sâu_phương pháp thu mẫu) có đôi chút khác biệt. Để đồng bộ các tên này chúng ta cần đưa chúng về một định dạng giống nhau (OSDxx.số thứ tự). Tuy vậy, trong quá trình định dạng lại tên mẫu có thể xảy ra việc gán số thứ tự khác nhau giữa tên mẫu  trong bảng đọc trình tự và bảng môi trường. Bảng Name_list được tạo ra để đối chứng các tên mẫu trong 2 bảng và tên định dạng mới của chúng.

```{r}
env_Full_names <- rownames(envdata)
env_Short_names <- make.unique(sapply(strsplit(rownames(envdata), "_"), `[`,1))
env_Names <- cbind(Short_Names = env_Short_names, env_Full_names)
Sequencing_names <- cbind(Short_Names = sample.names, S_Full_Names = sample.namesfull)
Name_list <- merge(Sequencing_names, env_Names, all = T, all.x = TRUE)
head(Name_list)
write.table(Name_list, "/file directory/Name_list.txt",sep="\t")
```


###Tạo hồ sơ di truyền "phyloseq" cho các trình tự đầu vào
Format lại tên hàng của bảng môi trường
```{r}
rownames(envdata) <- make.unique(sapply(strsplit(rownames(envdata), "_"), `[`, 1))
rownames(envdata)
DAT <- sample_data(envdata)
head(DAT)
```

Tạo bảng Physeq
```{r}
ps1 <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), tax_table(taxaSp))
ps1 <- merge_phyloseq(ps1, DAT)
ps1
```

Đổi tên trình tự thành OTU, và kiểm tra thành phẩn của ps1
```{r}
#change nom seq/OTU
summary(taxa_names(ps1))
head(taxa_names(ps1))
tail(taxa_names(ps1))
taxa_names(ps1) <- paste0("OTU", seq(ntaxa(ps1)))
```

###OUTPUT
Trích bảng OTU từ ps1
```{r}
# Extract abundance matrix from the phyloseq object
OTU1 = as(otu_table(ps1), "matrix")
# transpose if necessary
if(taxa_are_rows(ps1)){OTU1 <- t(OTU1)}
# Coerce to data.frame
OTUdf = as.data.frame(OTU1)
```

Trích bảng TAX từ ps1
```{r}
# Extract abundance matrix from the phyloseq object
TAX = as(tax_table(ps1), "matrix")
# transpose if necessary
if(taxa_are_rows(ps1)){TAX <- t(TAX)}
# Coerce to data.frame
TAXdf = as.data.frame(TAX)
```

Trích bảng DAT từ ps1
```{r}
# Extract abundance matrix from the phyloseq object
DAT1 = as(sample_data(ps1), "matrix")
# transpose if necessary
if(taxa_are_rows(ps1)){DAT1 <- t(DAT1)}
```

Ghi các bảng OTU TAX DAT1 ra file .txt
```{r}
DAT1df = as.data.frame(DAT1)
write.table(OTU1, "/Users/nqhuynh/Documents/1_OTU1.txt", sep="\t")
write.table(TAX, "/Users/nqhuynh/Documents/1_TAX.txt", sep="\t")
write.table(DAT1, "/Users/nqhuynh/Documents/1_DAT.txt", sep="\t")
```
