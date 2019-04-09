#' Analysis of short tandem repeats (STRs) in a given region of any reference genome
#' @description  This function separates detected short tandem repeats (STRs) into different zones. These zones are either the hotspot zone defined by the double strand break maps of Pratto et al. (2014) or adjacent flanking zones (greyzones) left and right of the hotspots of user specified lengths. The parameters of the regions under study can be directly given in the function arguments or read in via either a BED-file or a position matrix.
#' @param seqName A character string which is the name of the given sequence file under study. Can also be set to "" in order to analyze a defined sequence from any reference genome such as the package BSgenome.Hsapiens.UCSC.hg19 for humans.
#' @param nr.STRs An integer value reflecting the minimum length of STRs to be searched for.
#' @param nr.mismatch An integer value reflecting the allowed number of mismatches of the short tandem repeats. By default it set to 0.
#' @param chrs A string reflecting the chromosome under study (starting with "chr" and adding either the integers from 1-22 or "X" respectively "Y" for the human genome). This argument can also be a vector of strings to study several chromosomes.
#' @param start.position An integer value reflecting the start position of the region to be analyzed. If set to \code{NA} the analysis starts from the beginning of the chromosome.
#' @param end.position An integer value reflecting the end position of the region to be analyzed. If set to \code{NA} the analysis is performed until the end of the chromosome.
#' @param reverse.comp A logical value by default \code{FALSE}. If set to \code{TRUE} then the reverse complement of the sequence is analyzed.
#' @param STR A character string for the nucleotide to be searched for. By default one searches for poly-As, hence set to "A".
#' @param lens.grey An integer value which is by default a vector of 6 integer values. These values represent the greyzones to be studied left and right from the hotspot regions.
#' @param bed_file A bed file containing the chromosomes, start, and end positions of the region(s) that should be analyzed.
#' @param pos_matrix A matrix or dataframe containing the chromosomes, start, and end positions of the region(s) that should be analyzed.
#' @param output_file The default is an empty string and does not save an output-file. The output will be saved if the parameter is changed to a user defined string excluding the extension (by default .bed).
#' @param species The human genome (version 19) is default but an alternative genome can be provided. For chimpanzees the parameter has to be BSgenome.Ptroglodytes.UCSC.panTro5 (given that the data is installed).
#' @param dsb_map The DSB map of the human genome (version 19) is default but an alternative DSB map from a different genome can be provided. This parameter needs to be a data frame with at least 3 columns that contains the chromosome, start, and end position of the DSB. The DSB map for chimpanzees is included in the package.
#'
#' @return The output of the function is a list with the following content:
#' \item{Sequence Name}{The chromosome with the starting and end position of the region under study is provided.}
#' \item{Reverse Complement}{An indicator whether the reverse complement was considered}
#' \item{Number of allowed Mismatches}{The number of allowed mismatches is provided.}
#' \item{Minimum Length}{The minimum length of the STR to be extracted is provided.}
#' \item{Number of Matches}{The total number of STR matches of the region is provided.}
#' \item{Length of STR stretch in bp}{A vector containing the length of STRs per match is provided.}
#' \item{Start positions}{The starting positions of the STRs are provided.}
#' \item{Zone}{The zones where the STR is found are provided. 1 reflects within a hotspot, the last integer reflects that it is outside, and the integers between these two reflect the given flanking regions starting with 2 as the next closest region to the hotspot.}
#' @return A BED file with the chromosomes, start, and end position of the STRs, length of the STR stretch, the zone where the STR was found, and the specified region that was analyzed are given as columns.

#' @author Philipp Hermann, \email{philipp.hermann@@jku.at}, Monika Heinzl, \email{monika.heinzl@@edumail.at}
#' Angelika Heissl, Irene Tiemann-Boege, Andreas Futschik
#' @seealso \code{\link{getflank2}}, \code{\link{STR_detection}}
#' @references Heissl, A., et al. (2018) Length asymmetry and heterozygosity strongly influences the evolution of poly-A microsatellites at meiotic recombination hotspots. doi: https://doi.org/10.1101/431841
#'
#' Pratto, F., et al. (2014). Recombination initiation maps of individual human genomes. Science, 346(6211).
#'
#' Kuhn RM, et al. (2013) The UCSC genome browser and associated tools, Brief. Bioinform., 14, 144-161.
#' @examples
#' data(chr6_1580213_1582559)
#' STR_analysis(seqName = chr6_1580213_1582559, nr.STRs = 10, nr.mismatch = 0, chrs = "chr6",
#' STR = "A", lens.grey = 0:1*100, start.position = 1580213, end.position = 1582559,
#' reverse.comp = FALSE,
#' species = BSgenome.Hsapiens.UCSC.hg19::Hsapiens, dsb_map = STRAH::dsb_map)
#' \donttest{
#' STR_analysis(nr.STRs = 10, nr.mismatch = 0, chrs = "chr22", STR = "A", lens.grey = 0:1*100,
#' start.position = 30000000, end.position = 31000000, reverse.comp = FALSE,
#' species = BSgenome.Hsapiens.UCSC.hg19::Hsapiens, dsb_map = STRAH::dsb_map)
#' # If you want to use the function with a different reference genome
#' # make your choice and install it before:
#' if(requireNamespace("BSgenome.Ptroglodytes.UCSC.panTro5")) {
#' STR_analysis(nr.STRs = 10, nr.mismatch = 0, chrs = "chr22", STR = "A", lens.grey = 0:5*1000,
#' start.position = 30000000, end.position = 31000000, reverse.comp = FALSE,
#' species = BSgenome.Ptroglodytes.UCSC.panTro5::BSgenome.Ptroglodytes.UCSC.panTro5,
#' dsb_map = STRAH::dsb_map_chimp_full)
#' }
#' }
#' @keywords datasets array list methods univar
#' @export
#'
STR_analysis = function(seqName, nr.STRs = 10, nr.mismatch = 0, chrs, STR = "A", lens.grey = 0:5*1000, start.position = NA,
                        end.position = NA, reverse.comp = FALSE, bed_file, pos_matrix, output_file,
                        species=BSgenome.Hsapiens.UCSC.hg19::Hsapiens, dsb_map = STRAH::dsb_map) {
 if(all(missing(seqName), missing(chrs), missing(bed_file), missing(pos_matrix))) {stop("Please only provide one of the parameters seqName, chrs, bed_file, or pos_matrix!")}
 if(!missing(seqName) & missing(chrs)) {stop("Please also provide chrs when you use a DNAStringSet-object!")}
 ind = length.As = pos.As = seq_name_list = nr_matches = original_region_list = list()
 df <- ""
 if(!missing(bed_file)){
   bed <- read.table(bed_file,header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
   chrs <- bed[,1]
   start.position <- bed[,2] +1
   end.position <- bed[,3]
 }
 else if(!missing(pos_matrix)) {
   # if((is.data.frame(pos_matrix) | is.matrix(pos_matrix))){
   chrs <- pos_matrix[,1]
   start.position <- pos_matrix[,2]
   end.position <- pos_matrix[,3]
   # }
 }

 index_chr_no_str <- vector(mode="logical", length=0)
 index_str = 1

 for(index in 1:length(chrs)) {
   if(any(!missing(bed_file), !missing(pos_matrix))){
     assign(paste0("results.", chrs[index]),
            STRAH::STR_detection(chrs = chrs[index], nr.mismatch = nr.mismatch, nr.STRs = nr.STRs, STR = STR,
                          start.position = start.position[index], end.position = end.position[index],
                          translated_regions=F, output_file = output_file, species = species)) # start.position = 1, end.position = chr.lengths[index],
   } else if(missing(seqName)) {
    assign(paste0("results.", chrs[index]),
           STRAH::STR_detection(chrs = chrs[index], nr.mismatch = nr.mismatch, nr.STRs = nr.STRs, STR = STR,
                                 start.position = start.position, end.position = end.position,
                                 translated_regions=F, output_file = output_file, species = species)) # start.position = 1, end.position = chr.lengths[index],
  } else {
    assign(paste0("results.", chrs[index]),
           STRAH::STR_detection(seqName = seqName, chrs = chrs[index], nr.mismatch = nr.mismatch, nr.STRs = nr.STRs, STR = STR,
                                 start.position = NA, end.position = NA, translated_regions=F, output_file = output_file,
                                 species = species)) # start.position = 1, end.position = chr.lengths[index],
   }
  temp2 <- which(unlist(get(paste("results.", chrs[index], sep = ""))[[6]])>= nr.STRs)
  if(length(temp2)==0) {
    index_chr_no_str <- append(index_chr_no_str,FALSE) # add index of chr where no STR is present
    next
  }
  else{
    index_chr_no_str <- append(index_chr_no_str,TRUE)
  }

  ind[[index_str]] = temp2
  length.As[[index_str]] = unname(unlist(get(paste("results.", chrs[index], sep = ""))[[6]][temp2]))
  pos.As[[index_str]] = unname(unlist(get(paste("results.", chrs[index], sep = ""))[[7]][temp2]))
  nr_matches[[index_str]] = length(length.As[[index_str]])
  index_str = index_str + 1
 }
if(length(which(index_chr_no_str == FALSE)) == length(chrs)){
  stop("Analysis aborted since no STRs were found in the regions!")
}
 else{ # only some chrs contain contain STRs or all
   length_chr <- which(index_chr_no_str == T)
 }
### We manually filtered the DSB Map to contain only those where AA-type is present ####
 # dsb_map = read.csv("DSB-map_Pratto_man.csv", header = T, sep = ";")
 pos.chr = sapply(length_chr, function(i) { # choose only those chromosomes where STRs are present
 #pos.chr = sapply(1:length(chrs), function(i) {
   which(dsb_map$chrom == as.character(chrs[i]))
   })

 if(length(length_chr) == 1) {
   pos.chr <- list(c(pos.chr))
 }

 if (is.data.frame(pos.chr) == TRUE | is.matrix(pos.chr) == TRUE){
   pos.chr <- lapply(apply(pos.chr, 2, list), unlist)
 }
 start.dsb.int    = lapply(1:length(pos.chr), function(j) {dsb_map[pos.chr[[j]],2]})
 end.dsb.int      = lapply(1:length(pos.chr), function(j) {dsb_map[pos.chr[[j]],3]})

#### For being within it has to be larger than the start minus 500bp and smaller than the end position plus 500bp. ####
#### The greyzone is restructured in 5 segments of length 2 kb increasing in the distance to the hotspot ####
#### Moreover it is not allowed to contain a hotspot ####
 within = lapply(1:length(pos.chr), function(j) {
   sapply(1:length(ind[[j]]), function(i) {
     return(any(pos.As[[j]][i] >= (start.dsb.int[[j]]-500) & pos.As[[j]][i] <= (end.dsb.int[[j]]+500)))
     })
   })
 Sys.sleep(0.2)

  for(len.grey in 2:length(lens.grey)) {
    assign(paste("greyzone",lens.grey[len.grey-1]/1000,lens.grey[len.grey]/1000,sep="_"), lapply(1:length(pos.chr), function(j) {
      !within[[j]] & sapply(1:length(ind[[j]]), function(i) {
        return(any((pos.As[[j]][i] < (start.dsb.int[[j]]-500-lens.grey[len.grey-1]) & (pos.As[[j]][i] >= (start.dsb.int[[j]]-500-lens.grey[len.grey]))) | (pos.As[[j]][i] <= (end.dsb.int[[j]]+500+lens.grey[len.grey-1]) & pos.As[[j]][i] > (end.dsb.int[[j]]+500+lens.grey[len.grey]))))
        })
      }))
    message(paste("Number of greyzones finished: ", len.grey-1, " (of ", length(lens.grey)-1, ")", sep = ""))


 }## end lens.grey

 message("")
 flush.console()

 name.greyzones = unlist(strsplit(paste("greyzone", lens.grey[-length(lens.grey)]/1000, lens.grey[-1]/1000, collapse = " ", sep = "_"), split = " "))
 code.zone = lapply(1:length(get(name.greyzones[1])), function(x){
   return(sapply(1:length(get(name.greyzones[1])[[x]]), function(y) {
    if(within[[x]][y]) {
      return(1)
      } # "within"
    else for(k in 1:length(name.greyzones)) {
      if(get(name.greyzones[k])[[x]][[y]]) {return(k+1)}
    }
    return(length(name.greyzones)+2)}))
   })

 counter = 1
 for (s in length_chr) {
   if(missing(bed_file) && missing(pos_matrix)){
     name_file <- paste0(chrs[s], ":", start.position, "-", end.position)
  }
  else{
   name_file <- paste0(chrs[s], ":", start.position[s], "-", end.position[s])
  }
  seq_name_list <- append(seq_name_list, name_file)
  header <- c("chr", "start_STR", "end_STR", "len_STR", "zones","chr_start_stop")


  if(!missing(output_file)){
    if (is.data.frame(df) == FALSE){
      write.table(rbind(header), paste0(output_file, ".bed"), sep = "\t", col.names = FALSE, row.names = FALSE, quote=FALSE)
    }
    df <- data.frame(chr_str = rep(chrs[s], length(pos.As[[counter]])),start_str = pos.As[[counter]]-1, end_str = pos.As[[counter]] + length.As[[counter]]-1, len_str = length.As[[counter]], zones = code.zone[[counter]], chr_name = rep(name_file, length(pos.As[[counter]])))

    write.table(df, paste0(output_file, ".bed"), col.names=FALSE, sep="\t", row.names = FALSE, append=TRUE, quote=FALSE)
   }
  counter = counter + 1
 }

  output <- list("Sequence Name" = seq_name_list,
                  "Reverse Complement" = reverse.comp, "Number of allowed Mismatches" = nr.mismatch, "Minimum Length" = nr.STRs,
                  "Number of Matches" = nr_matches, "Length of STR stretch in bp" = length.As, "Start positions" = pos.As, "Zone" = code.zone)
  output <- lapply(output,FUN=function(x) {
   if(length(x) == length(unlist(seq_name_list)) & length(x) != 1){
     names(x) <- unlist(seq_name_list)
   }
 return(x)
 })
 closeAllConnections()
 return(output)
  }

