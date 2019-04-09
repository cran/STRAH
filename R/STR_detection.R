#' Detection of short tandem repeats (STRs) in a given region of any reference genome
#' @description  This function searches for short tandem repeats (STRs) in a specified region of any reference genome. The parameters of the regions under study can be directly given in the function arguments or read in via either a BED-file or a position matrix. We recommend to search for STRs of minimum length 6. Options to save the output or usage of any reference genome are provided.
#' @param seqName A character string which is the name of the given sequence file under study. Can also be set to "" in order to analyze a defined sequence from any reference genome such as the package BSgenome.Hsapiens.UCSC.hg19 for humans.
#' @param chrs A string reflecting the chromosome under study (starting with "chr" and adding either the integers from 1-22 or "X" respectively "Y" for the human genome). This argument can also be a vector of strings to study several chromosomes.
#' @param start.position An integer value reflecting the start position of the region to be analyzed. If set to \code{NA} the analysis starts from the beginning of the chromosome.
#' @param end.position An integer value reflecting the end position of the region to be analyzed. If set to \code{NA} the analysis is performed until the end of the chromosome.
#' @param nr.STRs An integer value as the minimum length of STRs to be detected.
#' @param nr.mismatch An integer value reflecting the allowed number of mismatches of the short tandem repeats. By defaults set to 0.
#' @param reverse.comp A logical value by default \code{FALSE}. If set to \code{TRUE} then the reverse complement of the sequence is analyzed.
#' @param STR A character string for the nucleotide to be searched for. By default one searches for poly-As, hence set to "A".
#' @param bed_file A bed file containing the chromosomes, start, and end positions of the region(s) that should be analyzed.
#' @param pos_matrix A matrix or dataframe containing the chromosomes, start, and end positions of the region(s) that should be analyzed.
#' @param species The human genome (version 19) is default but an alternative genome can be provided. For chimpanzees the parameter has to be BSgenome.Ptroglodytes.UCSC.panTro5 (given that the data is installed).
#' @param translated_regions A logical value by default \code{FALSE}. If set to \code{TRUE} then the function assumes that the parameters start.position and end.position were translated by some tool (e.g. liftOver) from one species to another. The untranslated and translated positions are included in the output.
#' @param output_file The default is an empty string and does not save an output-file. The output will be saved if the parameter is changed to a user defined string excluding the extension (by default .bed).
#'
#' @return The output of the function is a list with the following content:
#' \item{Sequence name}{The chromosome with the start and end position of the region under study is provided. If \code{translated_region} is set to \code{TRUE}, the interval will be the translated region.}
#' \item{Sequence name (untranslated)}{Only if \code{translated_region} is set to \code{TRUE}, then the untranslated region (chromosome with the starting and end position) is provided.}
#' \item{Reverse complement}{An indicator whether the reverse complement was considered}
#' \item{Number of allowed mismatches}{The number of allowed mismatches is provided.}
#' \item{Minimum length}{The minimum length of the STR to be extracted is provided.}
#' \item{Number of matches}{The total number of STR matches of the region is provided.}
#' \item{Length of STR stretch in bp}{A vector containing the length of STRs per match is provided.}
#' \item{Start positions}{The starting positions of the STRs are provided.}
#' \item{Matched segments}{The matched segments of the STRs are provided.}
#' @return A BED file with chromosomes, start, and end position of the STRs, length of the STR stretch, the matched segments, and the specified region (untranslated and translated) that was analyzed are given as columns.
#' @author Philipp Hermann, \email{philipp.hermann@@jku.at}, Monika Heinzl, \email{monika.heinzl@@edumail.at}
#' Angelika Heissl, Irene Tiemann-Boege, Andreas Futschik
#' @seealso \code{\link{getflank2}}, \code{\link{STR_analysis}}
#' @examples
#' data(chr6_1580213_1582559)
#' STR_detection(seqName = chr6_1580213_1582559, chrs = "chr6", start.position = 1580213,
#' end.position = 1582559, nr.STRs = 10, nr.mismatch = 0, reverse.comp = FALSE, STR = "A",
#' species = BSgenome.Hsapiens.UCSC.hg19::Hsapiens, translated_regions=FALSE)
#' \donttest{
#' STR_detection(chrs = "chr22", start.position = 30000000, end.position = 31000000,
#' nr.STRs = 10, nr.mismatch = 0, reverse.comp = FALSE, STR = "A",
#' species=BSgenome.Hsapiens.UCSC.hg19::Hsapiens, translated_regions=FALSE)
#' # If you want to use the function with a different reference genome
#' # make your choice and install it before:
#' if(requireNamespace("BSgenome.Ptroglodytes.UCSC.panTro5")) {
#' STR_detection(chrs = "chr1", start.position =222339618, end.position = 222339660,
#' nr.STRs = 10, nr.mismatch = 0, reverse.comp = FALSE, STR = "A",
#' species = BSgenome.Ptroglodytes.UCSC.panTro5::BSgenome.Ptroglodytes.UCSC.panTro5)
#' }
#' }
#' @references Heissl, A., et al. (2018) Length asymmetry and heterozygosity strongly influences the evolution of poly-A microsatellites at meiotic recombination hotspots. doi: https://doi.org/10.1101/431841
#'
#' Pratto, F., et al. (2014). Recombination initiation maps of individual human genomes. Science, 346(6211).
#'
#' Kuhn RM, et al. (2013) The UCSC genome browser and associated tools, Brief. Bioinform., 14, 144-161.
#'
#' @keywords datasets array list methods univar
#' @export

STR_detection = function(seqName, chrs, start.position = NA, end.position = NA, bed_file,
                         pos_matrix, nr.STRs,  nr.mismatch = 0, reverse.comp = F, STR = "A",
                         species=BSgenome.Hsapiens.UCSC.hg19::Hsapiens, translated_regions=F, output_file) {
  min.nr = nr.STRs
  start.position = max(start.position, 1, na.rm = T)
  df <- ""
  if(!missing(seqName)) {
    if(is(seqName, "DNAStringSet")) {
      sequence = seqName
    } else {
      sequence = Biostrings::readDNAStringSet(seqName, format = "fasta")
    }
  }
  else if(!missing(bed_file)){
    bed <- read.table(bed_file,header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
    chrs <- bed[,1]
    start.position <- bed[,2] + 1
    end.position <- bed[,3]
    if (translated_regions == TRUE){
      original_region <- bed[,4]
    }
    sequence <- Biostrings::DNAStringSet(getflank2(species = species, chrs = chrs, start.position = start.position, end.position = end.position))
  }
  else if(!missing(pos_matrix)){
    chrs <- pos_matrix[,1]
    start.position <- pos_matrix[,2]
    end.position <- pos_matrix[,3]
    sequence <- Biostrings::DNAStringSet(getflank2(species = species, chrs = chrs, start.position = start.position, end.position = end.position))
  }
  else {
    sequence <- Biostrings::DNAStringSet(getflank2(species = species, chrs = chrs, start.position = start.position, end.position = end.position))
    }
  if(reverse.comp) {
    sequence <- Biostrings::reverseComplement(s)
    STR <- toString(Biostrings::reverseComplement(Biostrings::DNAString(STR)))
  }


  seq_name_list <- list()
  nr.matches <- list()
  start.position_list <- list()
  nr.STRs.c_list <- list()
  matches_list <- list()
  original_region_list <- list()

  for (s in 1:(length(sequence))){
    Sys.sleep(0.2)
    if(missing(seqName)) {
      if(missing(bed_file) && missing(pos_matrix)){
        name <- paste0(chrs[s], ":", formatC(start.position[s], format = "fg"), "-", formatC(end.position[s], format = "fg"))
      }
      else{
        name <- paste0(chrs[s], ":", formatC(start.position[s], format = "fg"), "-", formatC(end.position[s], format = "fg"))
      }
    } else {
      name = names(seqName)
    }
    if(missing(seqName)) {
      message(paste0(name, " of ", toString(species), " is under study!"),"\r",appendLF=TRUE)
    } else {
        message(paste0(name, " is under study!"),"\r",appendLF=TRUE)
    }

    tp = unlist(strsplit(Biostrings::toString(sequence[s,]), split  ="")) == STR
    #motif <- paste(rep(STR, nr.STRs), collapse = "")
    #motif <- paste(STR, "{", nr.STRs, ",}")
    # tp_2 <- gregexpr(paste(motif, "+", sep=""), Biostrings::toString(s))

    #tp_2 <- Biostrings::matchPattern(motif, Biostrings::toString(s), max.mismatch = nr.mismatch)
    #print(tp_2)

    # tp_2 <- agrep( paste0(motif, '+'), Biostrings::toString(s),
    #                max = list(sub = nr.mismatch),  ignore.case = T, fixed=FALSE, useBytes = TRUE)
    # print(tp_2)

    if (length(which(tp)) == 0 | (length(which(tp))-(nr.STRs-1)) <= 0){
      message(paste0("No STRs are contained in ", name),"\r",appendLF=TRUE)
      message("","\r",appendLF=TRUE)
      next
    }
    tp2 = which(
      sapply(1:(length(tp)-(nr.STRs-1)), function(i) {
        sum(tp[i:(i+(nr.STRs-1))])
      }
      )>= (nr.STRs-nr.mismatch))
    start.pos = tp2[tp[tp2]]
    if(length(start.pos)==0) {
      message(paste0("No STRs are contained in ", name),"\r",appendLF=TRUE)
      next
    }

    add = c(ave(start.pos, cumsum(c(F, diff(start.pos) > 1)), FUN=seq_along) - 1,0)
    nr.STRs.c = add[which(add==0)-1]+nr.STRs


    #### remove duplicates based on the difference in the starting position ####
    ind.double = which(diff(start.pos) == 1)+1
    if (length(ind.double) != 0){ # only if there are more than one start positions of STRs
      start.pos <- start.pos[-ind.double]
    }
    matches = sapply(1:length(start.pos), function(k) {
    # if STR ends at last position of interval
    print(sequence)
    if ((start.pos[k]+nr.STRs.c[k]-1) >= Biostrings::width(sequence[s,])){
      end_pos_STR <- width(sequence[s,]) # end of interval is end of STR
    }
    else{
      end_pos_STR <- start.pos[k]+nr.STRs.c[k]-1
    }

      temp = unlist(Biostrings::subseq(sequence[s,], start.pos[k], end_pos_STR))
      ### we do not want to end a stretch with anything unequal to STR so as long as the last character is unequal to it, this character will be removed ####
      while(Biostrings::toString(temp[nchar(Biostrings::toString(temp))]) != STR) {
        nr.STRs.c[k] = nr.STRs.c[k]-1
        temp = unlist(Biostrings::subseq(temp, 1, nchar(Biostrings::toString(temp))-1))
      }
      # if STR starts at first position of interval
      if (start.pos[k] <= 5){
        start <- start.pos[k]
      }
      else{
        start <- start.pos[k]-5
      }
      # if STR ends with last or less than 5nt before end position of interval
      if ((end_pos_STR+5)>= width(sequence[s,])){
        end <- width(sequence[s,])
      }
      else{
        end <- end_pos_STR+5
      }
      temp.2 = unlist(Biostrings::subseq(sequence[s,], start, end))

      #print(temp.2)
      return(Biostrings::toString(temp.2))
    })
    print(matches)
    message("","\r",appendLF=TRUE)

    seq_name_list <- append(seq_name_list, name)
    nr.matches<- append(nr.matches, length(start.pos))
    nr.STRs.c_list <- append(nr.STRs.c_list,list(nr.STRs.c))
    start.position_list <- append(start.position_list, list(start.pos))
    matches_list <- append(matches_list, list(matches))

    if (translated_regions){
      original_region_list <- append(original_region_list, original_region[s])
      header <- c("chr", "start_STR", "end_STR", "len_STR", "sequence_STR", "translated chr_start_stop", "untranslated chr_start_stop")
    }
    else{
      header <- c("chr", "start_STR", "end_STR", "len_STR", "sequence_STR", "chr_start_stop")
    }

    flush.console()
    if(!missing(output_file)){
      if (is.data.frame(df) == FALSE){
        write.table(rbind(header),paste0(output_file, ".bed"), sep = "\t", col.names = FALSE, row.names = FALSE, quote=FALSE)
      }
      if (translated_regions){
        df <- data.frame(chr_str = rep(chrs[s], length(start.pos)),start_str = start.pos-1, end_str = start.pos + nr.STRs.c-1, len_str = nr.STRs.c, sequence_str = matches, chr_name = rep(name, length(start.pos)), chr_name_original = rep(original_region[s], length(start.pos)))
      }
      else{
        df <- data.frame(chr_str = rep(chrs[s], length(start.pos)),start_str = start.pos-1, end_str = start.pos + nr.STRs.c-1, len_str = nr.STRs.c, sequence_str = matches, chr_name = rep(name, length(start.pos)))
      }
      write.table(df, paste0(output_file, ".bed"), col.names=FALSE, sep="\t", row.names = FALSE, append=TRUE, quote=FALSE)

    }
  }
  if (translated_regions){
    output <- list("Sequence name (translated)" = seq_name_list, "Sequence name (untranslated)" = original_region_list, "Reverse complement" = reverse.comp, "Number of allowed mismatches" = nr.mismatch, "Minimum length" = min.nr, "Number of matches" = nr.matches,
                   "Length of STR stretch in bp" = nr.STRs.c_list, "Start positions" = start.position_list, "Matched segments" = matches_list)
  }
  else{
    output <- list("Sequence Name" = seq_name_list, "Reverse complement" = reverse.comp, "Number of allowed mismatches" = nr.mismatch, "Minimum length" = min.nr, "Number of matches" = nr.matches,
                   "Length of STR stretch in bp" = nr.STRs.c_list, "Start positions" = start.position_list, "Matched segments" = matches_list)
  }

    output <- lapply(output,FUN=function(x) {
    if(length(x) == length(unlist(seq_name_list)) & length(x) != 1){
      names(x) <- unlist(seq_name_list)
    }
    return(x)
  })
  closeAllConnections()
  return(output)
  }

