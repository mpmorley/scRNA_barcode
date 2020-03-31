library(ggplot2)
library(dplyr)
library(gridExtra)
library(stringdist)
library(readr)
library(optparse)

#############################################################################
# This code was modify from a series scripts written by Yogesh Goyal with Shweta Ramdas. Michale Morley moerged teh code 
# into a single and added the ability to run via cmd line parameters. 
#
#################################################################################



###########################################################
# Edit for Starcode params
###############################################################
samplerRun<- list()
samplerRun[['sub50_d8']] <- list(subSampleSize=50, d=8)
samplerRun[['sub40_d8']] <- list(subSampleSize=40, d=8)
samplerRun[['sub30_d8']] <- list(subSampleSize=30, d=8)
samplerRun[['sub30_d6']] <- list(subSampleSize=30, d=6)


subSampleSize <- 4000
  
########################################################################################
option_list = list(
  make_option(c("-s", "--shavedReads"), default=NULL, 
              help="ShavedReadsFile file name"),
  make_option(c("-w", "--whitelist"), default=NULL, 
              help="10x whitelist file [default= %default]"),
  make_option(c("-o", "--outdir"), type="character", default="./", 
              help="output dir [default= %default]"),
  make_option(c("-t", "--threads"), type="integer", default=4, 
              help="# of threads for starcode [default= %default]")
  
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#testing 
opt$whitelist<-'CellIDsWhitelist.tsv' 
opt$shavedReads <- 'EEM018/uniqueShavedReads.txt'
opt$outdir <- 'EEM018'
opt$threads <- 4
##################################################################################
#    DEfine some functions 
#
#
################################################################################
MakeLVPlot <- function(m){
  lvdistSamp.plot <- bind_rows(LVdist) %>% 
    group_by(subsamNum, lvdist) %>% 
    summarise(length(lvdist)) %>%
    group_by(subsamNum) %>% 
    mutate(totalNum = sum(`length(lvdist)`), fracLVDist = `length(lvdist)`/totalNum) %>% 
    ggplot(., aes(hamdist, fracLvDist)) +
    geom_bar(width = 0.5, stat = 'identity') +
    facet_grid(rows=vars(subsamNum)) +
    theme_classic()
  
}


###########################################################################
# Step one merge files with 10x Barcode 
#
#
###########################################################################

#whitelist <- read_table(file = opt$whitelist,col_names = 'cellID') 
#shavedReads <- read_table(opt$shavedReads, col_names = c('cellID','UMI','V3'))
shavedReads.whitelist = inner_join(read_table(opt$shavedReads, col_names = c('cellID','UMI','V3')), 
                                   read_table(file = opt$whitelist,col_names = 'cellID') , 
                                   by = "cellID")



############################################################################
# Step 3
#
############################################################################
temp3 <- shavedReads.whitelist %>% pull(V3) %>% unique() %>% substring(.,1,70)

n=10000
LVdist <- list()
for(r in 1:3){
  LVdist[[r]] = data.frame(lvdist=as.integer(stringdistmatrix(sample(temp3,n), method = "lv")))  
  LVdist[[r]]$subsamNum<-r
}
OrigBarcode.plot <- MakeLVPlot(LVdist)
ggsave(OrigBarcode.plot,file=paste0(opt$outdir,'OringalBarcodesSubsampleHammingPlots.pdf'), width = 12, height = 9)

###########################################################################################################################################
#                                       Step Four         
# This will runb starcode using the run param define at the top of the script. 
#
#
###########################################################################################################################################


lapply(samplerRun, function(s){
  shavedReads.whitelist %>% pull(V3) %>% substring(.,1,s[['subSampleSize']]) %>% as.data.frame(.) %>%
    write_tsv(., path= paste0(opt$outdir,'Barcodes_subSample_',s[['subSampleSize']],'.txt'),col_names=F)
  starfile <-paste0(opt$outdir,"Barcodes_subSample_",s[['subSampleSize']],".txt")
  if(file.exists(starfile)){
    print(paste0(starfile, ' exists, skipping starcode run. Please remove if you wish to re-run'))
  }else{
    system(paste0("starcode -i ",starfile," -d", s[['d']]," -o ",opt$outdir,"StarcodeBarcodes_",s[['subSampleSize']],"_d",s[['d']] ," --seq-id -s -t ", opt$threads))
    }
}
  )

###########################################################################################################################################
#         This will parse the starcode output and create matrix of barcodes for each starcode run
#
###########################################################################################################################################


res.list <- lapply(samplerRun, function(s){
  res <- read.table(paste0(opt$outdir,"StarcodeBarcodes_",s[['subSampleSize']],"_d",s[['d']] ),col.names =c('barcode','n','index'), stringsAsFactors = F)
  name=paste0('sub',s[['subSampleSize']],"_d",s[['d']])
  tmp.res <- data.frame(cellID=shavedReads.whitelist$cellID) 
  for(i in 1:nrow(res)){
    inds <- as.integer(unlist(strsplit(res[i,'index'], split=",")))
    tmp.res[inds,name] <- res[i,'barcode'] 
  }
  tmp.res %>% select(-cellID)
}
)

FinalBarcodes <- bind_cols(shavedReads.whitelist,res.list)



#################################################
# Step five 
# Sample barcodes from Starcode output and plot the Levenshtein 
####################################################


set.seed(2059)

#loop over the names from the samplerRun and write each starcode result to a file 
Starcode.LVdist <- list()
for(name in names(samplerRun)){
  
  col <- quo(name)
  samp <- FinalBarcodes %>% pull( !!col) %>% unique() %>% sample(.,subSampleSize)
  Starcode.LVdist[[name]] = data.frame(hamdist=as.integer(stringdistmatrix(samp), method = "lv"))  
  Starcode.LVdist[[name]]$subsamNum<-name
  

  finalResults <- 
    FinalBarcodes %>% select(cellID,!!col) %>% 
    distinct() %>%
    group_by(cellID) %>% 
    summarise(barcode = paste0(sub50_d8, collapse = ""), nbarcodes=n()) %>%
    select(cellID,barcode,nbarcodes) %>%
    ungroup()
  cloneid <- finalResults %>% ungroup() %>% select(barcode) %>% distinct() %>% mutate(cloneid = paste0('Clone_',row_number()))
  finalResults <- inner_join(finalResults, cloneid)
  write_tsv(finalResults, paste0(opt$outdir,'FinalBarcode_',name,'.tsv'))
  

}


Starcode.plot <- MakeLVPlot(Starcode.LVdist)
ggsave(file=paste0(opt$outdir,'StarcodeBarcodesHammingPlots.pdf'), width = 12, height = 9)



