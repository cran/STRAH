#' DNA-motif detection in a given DNAStringSet, a given DNA-sequence in fasta-format, or a specified region of any genome
#' @description This function searches for a given "DNA-motif" in a DNA-sequence. The argument \code{seqName} can be either a DNAStringSet object or it refers to a fasta-file. Additionally, we provide the option to specify a species, a chromosome, a start, and a stop position for a region of any reference genome to be analyzed. By default a region of the human genome is analyzed. Optionally, one can also specify the number of mismatches of the DNA-motif and whether the reverse complement has to be searched.
#' @param seqName A character string which can either be the name of a DNAStringSet object or a sequence name referring to a fasta-file to be analyzed. This argument can only be ignored if \code{chr} and \code{start.position} and \code{end.position} are specified.
#' @param chrs A character string reflecting the chromosome under study (starting with "chr" and adding either the integers from 1-22 or "X" respectively "Y" for the human chromosome). This argument can also be a vector of strings to study several chromosomes.
#' @param start.position An integer value reflecting the start position of the region to be analyzed. If set to \code{NA} the analysis starts from the beginning of the chromosome.
#' @param end.position An integer value reflecting the end position of the region to be analyzed. If set to \code{NA} the analysis is performed until the end of the chromosome.
#' @param motif A character string reflecting the specified DNA-motif to be searched for in the DNA-sequence.
#' @param nr.mismatch This integer specifies the number of allowed mismatches when searching for the specified DNA-motif.
#' @param reverse.comp A logical value, by default \code{FALSE}, which enables to search the reverse complement of the sequence if set to \code{TRUE}.
#' @param species The human genome (version 19) is default but an alternative genome can be provided. For chimpanzees the parameter has to be BSgenome.Ptroglodytes.UCSC.panTro5 (given that the data is installed).
#' @param print.status A logical value reflecting whether the current status of the worked sequence (relative to the sequence length) is printed (\code{TRUE}) or not (\code{FALSE}).
#' @author Philipp Hermann, \email{philipp.hermann@@jku.at}, Monika Heinzl, \email{monika.heinzl@@edumail.at}
#' Angelika Heissl, Irene Tiemann-Boege, Andreas Futschik
#' @return The output of the function is a list with the following content:
#' \item{Species}{The name of the species under study}
#' \item{Sequence Name}{The name of the region under study}
#' \item{Reverse Complement}{Indicator whether the reverse complement was searched}
#' \item{Number of Matches}{The frequency of found DNA-motifs in the region under study}
#' \item{Start Positions of Matches}{The start positions of the found DNA-motifs}
#' \item{Number of allowed Mismatches}{The number of allowed mismatches when searching for the DNA-motif}
#' \item{Matched Segments}{The list of the segments containing the DNA-motif}
#' @keywords datasets array list methods univar
#' @examples
#' data(chr6_1580213_1582559)
#' motif_detection(seqName = chr6_1580213_1582559, start.position = NA, end.position = NA,
#' motif = "CCNCCNTNNCCNC", nr.mismatch = 1, reverse.comp = FALSE, print.status = FALSE)
#'
#' \donttest{
#' motif_detection(chrs = "chr6", start.position = 1580213, end.position = 1582559,
#' motif = "CCNCCNTNNCCNC", nr.mismatch = 1, reverse.comp = FALSE, print.status = FALSE)
#' # If you want to use the function with a different reference genome
#' # make your choice and install it before:
#' if(requireNamespace("BSgenome.Ptroglodytes.UCSC.panTro5")) {
#' motif_detection(chrs = "chr1", start.position =222339618, end.position = 222339660,
#' motif = "A", nr.mismatch = 0, reverse.comp = FALSE, print.status = FALSE,
#' species = BSgenome.Ptroglodytes.UCSC.panTro5::BSgenome.Ptroglodytes.UCSC.panTro5)
#' }
#' }
#' @seealso \code{\link{getflank2}}
#' @references Heissl, A., et al. (2018) Length asymmetry and heterozygosity strongly influences the evolution of poly-A microsatellites at meiotic recombination hotspots. doi: https://doi.org/10.1101/431841
#' @export
motif_detection = function(seqName, chrs, start.position, end.position, motif, nr.mismatch = 0, reverse.comp = F,
                           print.status = T, species=BSgenome.Hsapiens.UCSC.hg19::Hsapiens) {
  if(missing(seqName) & missing(chrs)) {stop("Please provide either seqName or chrs!")}
  if(!missing(seqName)) {
    if(!is(seqName, "DNAStringSet")) {
      temp = Biostrings::readDNAStringSet(seqName, format = "fasta")
    } else {temp = seqName}
  } else {
    temp = Biostrings::DNAStringSet(getflank2(species = species, chrs = chrs, start.position = start.position, end.position = end.position))
    seqName = paste(chrs, start.position, end.position, sep = "_")

  }
  if(reverse.comp) {temp = Biostrings::reverseComplement(temp); motif = Biostrings::toString(Biostrings::reverseComplement(Biostrings::DNAString(motif)))}
  pos.n = which(unlist(strsplit(motif, "")) != "N")
  checks = start.pos = motifs = c()
  for(i in 1:(length(temp[[1]])-(nchar(motif)-1))) {
    sub = Biostrings::toString(Biostrings::subseq(temp, start = i, end=i+nchar(motif)-1))
    temps = c()
    for(ind in pos.n) {
      temps = c(temps, Biostrings::substr(sub,ind,ind) == Biostrings::substr(motif,ind,ind))
    }
    #### combn creates all possible combination of position comparisons between motif and sequence ####
    #### these combinations take into account the number of mismachtes such that only vectors of size e.g. -1 for 1 mismatch are generated ####
    #### for these vectors all possible subsamples are drawn with combn ####
    #### the colSums function counts the correct matches between sample(sequence) and motif #####
    #### this sum is compared to the length of the motif - number of mismatches #####
    #### we aim to have any of these to be exactly of desired length wherefore the resulting value will search for any T in the set of subsamples #####
    temp2 = any(colSums(combn(x = temps, m = length(temps) - nr.mismatch)) == (length(temps)-nr.mismatch))
    if(print.status) {print(paste(round(i/(length(temp[[1]])-nchar(motif)+1)*100,1), "%", sep = ""))}
    checks = c(checks, temp2)
    if(temp2) {start.pos = c(start.pos, i); motifs = c(motifs, Biostrings::toString(Biostrings::subseq(temp, i, i+nchar(motif)-1)))}

  }
  return(list("Species" = toString(species), "Sequence Name" = seqName, "Reverse Complement" = reverse.comp, "Number of Matches" = sum(checks), "Start Positions of Matches" = start.pos, "Number of allowed Mismatches" = nr.mismatch, "Matched Segments" = motifs))
}


