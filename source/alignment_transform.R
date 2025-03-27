args <- commandArgs(trailingOnly = TRUE)
fileName <- args[1]
fileOut <- args[2]


con=file(fileName,open="r")
line=readLines(con)
close(con)
stk <- data.frame()
stk[1,1:2] <- 0
colnames(stk) <- c("ID","seq")
count <- 1
for(i in 1:length(line)){
  if(line[i] == ""){next}
  if(substr(line[i],1,3) == "ASV"){
    word <- strsplit(line[i],split=" ",fixed=T)[[1]]
    for(j in 1:length(word)){
      if(substr(word[j],1,3) == "ASV"){
        stk[count,1] <-  word[j]
      }
      if(substr(word[j],1,1) == "-"){
        stk[count,2] <-  word[j]
        count <- count + 1
      }
    }
  }
}

fileConn<-  file(fileOut)
for(i in c(1:dim(stk)[1])){
  word <- gsub("\\.", replacement="-", stk[i,2])
  word <- gsub("U", replacement="T", word)         
  write(paste0(">",stk[i,1]),file="ssu_align.fasta",append=TRUE)
  write(word, file="ssu_align.fasta", append=TRUE)
}
close(fileConn)