## This script is modified from https://gitmemory.com/issue/vegandevs/vegan/378/668767001
ggrare <- function(physeq_object, step = 10, label = NULL, colour = NULL, plot = TRUE, parallel = FALSE, se = TRUE) {
  x <- methods::as(phyloseq::otu_table(physeq_object), "matrix")
  if (phyloseq::taxa_are_rows(physeq_object)) { x <- t(x) }
  ## This script is adapted from vegan `rarecurve` function
  tot <- rowSums(x)
  S <- rowSums(x > 0)
  nr <- nrow(x)
  rarefun <- function(i) {
    #cat(paste("rarefying sample", rownames(x)[i]), sep = "\n")
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) {
      n <- c(n, tot[i])
    }
    y <- vegan::rarefy(x[i, ,drop = FALSE], n, se = se)
    if (nrow(y) != 1) {
      rownames(y) <- c(".S", ".se")
      return(data.frame(t(y), Size = n, Sample = rownames(x)[i]))
    } else {
      return(data.frame(.S = y[1, ], Size = n, Sample = rownames(x)[i]))
    }
  }
  if (parallel) {
    out <- parallel::mclapply(seq_len(nr), rarefun, mc.preschedule = FALSE)
  } else {
    out <- lapply(seq_len(nr), rarefun)
  }
  df <- do.call(rbind, out)
  # Get sample data
  if (!is.null(phyloseq::sample_data(physeq_object, FALSE))) {
    sdf <- methods::as(phyloseq::sample_data(physeq_object), "data.frame")
    sdf$Sample <- rownames(sdf)
    data <- merge(df, sdf, by = "Sample")
    labels <- data.frame(x = tot, y = S, Sample = rownames(x))
    labels <- merge(labels, sdf, by = "Sample")
  }
  # Add, any custom-supplied plot-mapped variables
  if ( length(colour) > 1 ) {
    data$colour <- colour
    names(data)[names(data) == "colour"] <- deparse(substitute(colour))
    colour <- deparse(substitute(colour))
  }
  if ( length(label) > 1 ) {
    labels$label <- label
    names(labels)[names(labels) == "label"] <- deparse(substitute(label))
    label <- deparse(substitute(label))
  }
  p <- ggplot2::ggplot(data = data,
                       ggplot2::aes_string(x = "Size",
                                           y = ".S",
                                           group = "Sample",
                                           colour = colour))
  p <- p + ggplot2::labs(x = "Sequence Sample Size", y = "Species Richness")
  if (!is.null(label)) {
    p <- p + ggplot2::geom_text(data = labels,
                                ggplot2::aes_string(x = "x",
                                                    y = "y",
                                                    label = label,
                                                    linetype = linetype),
                                size = 4, hjust = 0)
  }
  p <- p + ggplot2::geom_line()
  if (se) { ## add standard error if available
    p <- p +
      ggplot2::geom_ribbon(ggplot2::aes_string(ymin = ".S - .se",
                                               ymax = ".S + .se",
                                               colour = NULL,
                                               fill = colour),
                           alpha = 0.2)
  }
  invisible(p)
}

DADA2Adapter <- function(seqtab, taxa, out_path){
  n_ASVs <- ncol(seqtab)
  magnitude <- ceiling(log10(n_ASVs))
  ASVID <- sprintf(paste0("ASV%0", magnitude, "d"), 1:n_ASVs)
  
  ASVseq <- colnames(seqtab)
  colnames(seqtab) <- ASVID
  seqtab <- t(seqtab)
  
  rownames(taxa) <- ASVID[match(rownames(taxa), ASVseq)]
  taxa[is.na(taxa)] <- "unidentified"
  
  fileConn <- file(paste0(out_path,"/ASV.fasta"))
  writeLines(c(rbind(paste0(">", ASVID), ASVseq)), fileConn)
  close(fileConn)
  
  write.table(seqtab,file=paste0(out_path,"/ASV_table.txt"), quote=F, sep="\t", row.names=T)
  write.table(taxa,file=paste0(out_path,"/ASV_taxa_table.txt"), quote=F, sep="\t", row.names=T)
}

get_composition <- function(rawdata, taxa, level, metadata, grouping){
  for(i in 1:ncol(rawdata)){
    rawdata[,i] <- rawdata[,i]/colSums(rawdata)[i]*100
  }
  
  uniqtaxa <- unique(taxa[,level])
  taxaRA <- c()
  for(i in 1:length(uniqtaxa)){
    taxaRA <- c(taxaRA, sum(rawdata[rownames(taxa)[taxa[,level] == uniqtaxa[i]],]))
  }
  names(taxaRA) <- uniqtaxa
  sort_list <- sort(taxaRA, decreasing = T)
  sort_list <- sort_list[!names(sort_list) %in% c("unidentified","uncultured")]
  sort_list <- names(sort_list[1:10])
  compair <- data.frame()
  compair[1,1:4] <- 0
  colnames(compair) <- c("sample","taxa","percentage","group")
  count <- 1
  for(i in 1:ncol(rawdata)){
    for(x in sort_list){
      compair[count, 1] <- colnames(rawdata)[i]
      compair[count, 2] <- x
      compair[count, 3] <- sum(rawdata[rownames(taxa)[taxa[,level] == x],i])
      compair[count, 4] <- metadata[i, grouping]
      count <- count + 1
    }
    compair[count, 1] <- colnames(rawdata)[i]
    compair[count, 2] <- "Others"
    compair[count, 3] <- sum(rawdata[rownames(taxa)[!taxa[,level] %in% sort_list],i])
    compair[count, 4] <- metadata[i, grouping]
    count <- count + 1
  }
  compair$taxa <- factor(compair$taxa, levels = c(sort_list,"Others"))
  sorder <- compair[compair$taxa == sort_list[1],]
  sorder <- sorder$sample[order(sorder$percentage, decreasing = T)]
  compair$sample <- factor(compair$sample, levels = sorder)
  return(compair)
}

phyloseq_to_edgeR = function(physeq, group, method="RLE", ...){
  if( !taxa_are_rows(physeq) ){ physeq <- t(physeq) }
  x = as(otu_table(physeq), "matrix")
  # Add one to protect against overflow, log(0) issues.
  x = x + 1
  # Check `group` argument
  if( identical(all.equal(length(group), 1), TRUE) & nsamples(physeq) > 1 ){
    # Assume that group was a sample variable name (must be categorical)
    group = get_variable(physeq, group)
  }
  # Define gene annotations (`genes`) as tax_table
  taxonomy = tax_table(physeq, errorIfNULL=FALSE)
  if( !is.null(taxonomy) ){
    taxonomy = data.frame(as(taxonomy, "matrix"))
  } 
  # Now turn into a DGEList
  y = edgeR::DGEList(counts=x, group=group, genes=taxonomy, remove.zeros = TRUE, ...)
  # Calculate the normalization factors
  z = edgeR::calcNormFactors(y, method=method)
  # Check for division by zero inside `calcNormFactors`
  if( !all(is.finite(z$samples$norm.factors)) ){
    stop("Something wrong with edgeR::calcNormFactors on this data,
         non-finite $norm.factors, consider changing `method` argument")
  }
  # Estimate dispersions
  return(edgeR::estimateTagwiseDisp(edgeR::estimateCommonDisp(z)))
}

track_ab <- function(curr_taxa, rank ,taxa){
  ASVgroup <- rownames(taxa)[taxa[,rank] == curr_taxa]
  return(ASVgroup)
}

track_ta <- function(curr_taxa, taxa){
  new_id <- c()
  for(i in 1:length(curr_taxa)){
    header <- substr(tolower(colnames(taxa)[i]),1,1)
    if(i == 7){
      tail <- paste0(c(curr_taxa[6], curr_taxa[7]), collapse = '_')
      new_id <- c(new_id, paste0(c(header,'_',tail), collapse = ''))
    } else {
      new_id <- c(new_id, paste0(c(header,'_',curr_taxa[i]), collapse = ''))
    }
  }
  new_id <- paste0(new_id, collapse = '|')
  return(new_id)
}

LEfSe_preparation <- function(data, taxa, meta, factor_, pvl_filter = 0.05){
  
  relative_data <- data
  for(i in 1:ncol(relative_data)){
    relative_data[,i] <- relative_data[,i]/colSums(relative_data)[i]*100
  }
  
  exdata <- data
  exdata[exdata > 0] <- 1
  keep <- rownames(exdata)[rowSums(exdata) > ncol(exdata)*pvl_filter]
  relative_data <- relative_data[keep, ]
  taxa <- taxa[keep, ]
  
  transfromed <- data.frame()
  factorID <- unname(meta[,factor_])
  samID <- rownames(meta)
  transfromed <- rbind(factorID,samID)
  newID <- c()
  for(i in 1:7){
    ready_taxa <- unique(taxa[,i])
    ready_taxa <- ready_taxa[! ready_taxa %in%  c("unidentified","uncultured")]
    for(j in 1:length(ready_taxa)){
      temp <- taxa[taxa[,i] == ready_taxa[j],]
      if(is.null(nrow(temp))){word = temp[1:i]} else {word = temp[1,1:i]}
      word <- track_ta(word, taxa)
      newID <- c(newID, word)
      downstream_ASVid <- track_ab(ready_taxa[j], i, taxa)
      if(length(downstream_ASVid) == 1 ){
        append <-  unname(relative_data[downstream_ASVid,])
        append <- as.numeric(append)
        transfromed <- rbind(transfromed, append)
        #write(paste0(c(paste0(word,collapse = "|"), append )  , collapse = "    "),file="for_lefse.tsv",append=TRUE)
      } else {
        append <-  unname(colSums(relative_data[downstream_ASVid,]))
        transfromed <- rbind(transfromed, append)
        #write(paste0(c(paste0(word,collapse = "|"), append )  , collapse = "    "),file="for_lefse.tsv",append=TRUE)
      }
    }
  }
  colnames(transfromed) <- NULL
  rownames(transfromed) <- c(factor_, "sampleID", newID)
  return(transfromed)
}

GMMviz <- function(module, meta, control_){
  case_ <- unique(meta$group)[unique(meta$group) != control_]
  sigdf <- data.frame()
  boxdf <- data.frame()
  covdf <- data.frame()
  for(i in 1:nrow(module)){
    T3 <- as.numeric(module[i,rownames(meta)[meta$group == case_]])
    T1 <- as.numeric(module[i,rownames(meta)[meta$group == control_]])
    tmp <- data.frame('group'=c(rep(case_,length(T3)),rep(control_,length(T1))), 'value'=c(T3,T1))
    tmp$module <- rownames(module)[i]
    boxdf <- rbind(boxdf, tmp)
    res <- wilcox.test(T1, T3)
    sigdf <- rbind(sigdf, data.frame('module'=rownames(module)[i],
                                     'log2FC'=log(mean(T3)+1,2) - log(mean(T1)+1,2),
                                     'p'=res$p.value))
    T3 <- as.numeric(coverage[i,rownames(meta)[meta$group == case_]])*100
    T1 <- as.numeric(coverage[i,rownames(meta)[meta$group == control_]])*100
    tmp <- data.frame('group'=c(control_, case_), 'mean'=c(mean(T1), mean(T3)), 'std'=c(sd(T1),sd(T3)))
    tmp$module <- rownames(module)[i]
    covdf <- rbind(covdf, tmp)
  }
  sigdf$adj_p <- p.adjust(sigdf$p, 'fdr')
  sigdf$group <- ifelse(sigdf$log2FC < 0, control_, case_)
  return(list('significance_df'=sigdf,'relative_abundance_df'=boxdf,'coverage_df'=covdf))
}

make_edgetable <- function(cor_matrix, data){
  edge_table <- edge(cor_matrix, rownames(data))
  if(nrow(edge_table) == 0){
    print("Do not had edge in this condition")
    return(edge_table)
  }
  edge_table$ID <- 0:(dim(edge_table)[1]-1)
  edge_table$Type <- "undirected"
  edge_table$pn <- ifelse(edge_table$weight > 0, 1, 0)
  return(edge_table)
}

make_nodetable <- function(edge_table, data, taxa){
  tag <- unique(c(edge_table$Source,edge_table$Target))
  node_table <- data.frame()
  node_table[1,1:8] <- 1
  colnames(node_table) <- c("ID","label","Phylum","Class","Order","Family","Genus","size")
  
  if(length(tag) == 0){
    print("Do not had edge in this condition")
    return(node_table)
  }
  count <- 1
  for(x in tag){
    node_table[count,1] <- x
    node_table[count,2] <- taxa$Order[rownames(taxa) == x]
    node_table[count,3] <- taxa$Phylum[rownames(taxa) == x]
    node_table[count,4] <- taxa$Class[rownames(taxa) == x]
    node_table[count,5] <- taxa$Order[rownames(taxa) == x]
    node_table[count,6] <- taxa$Family[rownames(taxa) == x]
    node_table[count,7] <- taxa$Genus[rownames(taxa) == x]
    node_table[count,8] <- sum(data[rownames(data) == x,]) / sum(data) * 100
    count <- count + 1
  }
  return(node_table)
}

make_net <- function(df){
  name = unique(c(as.character(df$Source),as.character(df$Target)))
  mat <- matrix(0, nrow = length(name), ncol = length(name))
  rownames(mat) <- name
  colnames(mat) <- name
  for(i in c(1:dim(df)[1])){
    anchor1 <- as.character(df$Source[i])
    anchor2 <- as.character(df$Target[i])
    weigt <- df$weight[i]
    mat[anchor1, anchor2] <- weigt
    mat[anchor2, anchor1] <- weigt
  }
  return(as.matrix(mat))
}
assignModul <- function(nodeT,modul,Method,cutoff){
  nodeT$modul <- "uncluster"
  for(x in nodeT$ID){
    if(x %in% modul$names){
      nodeT[nodeT$ID == x,"modul"] <- paste0(c("modul",membership(modul)[x]),collapse = "_")
    } else {
      next
    }
  }
  sortModul <- sort(table(nodeT$modul),decreasing = T)
  count <- 1
  for(x in names(sortModul)){
    if(x == "uncluster"){next}
    nodeT$modul[nodeT$modul == x] <- paste0(c(Method,count),collapse = "_")
    count <- count + 1
  }
  sortModul <- sort(table(nodeT$modul),decreasing = T)
  nodeT$modul[nodeT$modul %in% names(sortModul[sortModul < cutoff])] <- "uncluster"
  colnames(nodeT)[ncol(nodeT)] <- Method
  return(nodeT)
}