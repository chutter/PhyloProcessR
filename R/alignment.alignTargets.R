#' @title alignTargets
#'
#' @description Function for batch trimming a folder of alignments, with the various trimming functions available to select from
#'
#' @param targets.to.align path to a folder of sequence alignments in phylip format.
#'
#' @param output.directory available input alignment formats: fasta or phylip
#'
#' @param min.taxa contigs are added into existing alignment if algorithm is "add"
#'
#' @param subset.start available output formats: phylip
#'
#' @param subset.end algorithm to use: "add" add sequences with "add.contigs"; "localpair" for local pair align. All others available
#'
#' @param threads path to a folder of sequence alignments in phylip format.
#'
#' @param memory give a save name if you wnat to save the summary to file.
#'
#' @param overwrite path to a folder of sequence alignments in phylip format.
#'
#' @param resume contigs are added into existing alignment if algorithm is "add"
#'
#' @param quiet algorithm to use: "add" add sequences with "add.contigs"; "localpair" for local pair align. All others available
#'
#' @param mafft.path algorithm to use: "add" add sequences with "add.contigs"; "localpair" for local pair align. All others available
#'
#' @return an alignment of provided sequences in DNAStringSet format. Also can save alignment as a file with save.name
#'
#' @examples
#'
#' your.tree = ape::read.tree(file = "file-path-to-tree.tre")
#' astral.data = astralPlane(astral.tree = your.tree,
#'                           outgroups = c("species_one", "species_two"),
#'                           tip.length = 1)
#'
#'
#' @export

alignTargets = function(targets.to.align = NULL,
                        target.file = NULL,
                        output.directory = "alignments",
                        algorithm = c("localpair", "globalpair"),
                        min.taxa = 4,
                        removal.threshold = 0.35,
                        subset.start = 0,
                        subset.end = 1,
                        adjust.direction = TRUE,
                        remove.reverse.tag = TRUE,
                        threads = 1,
                        memory = 1,
                        overwrite = FALSE,
                        quiet = TRUE,
                        mafft.path = NULL) {

  # #Debug setup
  # setwd("/Users/chutter/Dropbox/SharewithCarl")
  # targets.to.align = "Venom-Markers-Nov23_to-align.fa"
  # target.file = "venom_loci_updated_Mar12_cdhit95_duplicate_exons_renamed_Feb2023_FINAL.fa"
  # output.directory = "alignments"
  #
  # #Main settings
  # subset.start = 0
  # subset.end = 1
  # min.taxa = 4
  # threads = 4
  # memory = 8
  # overwrite = TRUE
  # resume = FALSE
  # quiet = TRUE
  # adjust.direction = TRUE
  # algorithm = "localpair"
  # remove.reverse.tag = TRUE
  # removal.threshold = 0.35
  #
  # #program paths
  # mafft.path = "/Users/chutter/Bioinformatics/miniconda3/envs/PhyloProcessR/bin"

  #Same adds to bbmap path
  if (is.null(mafft.path) == FALSE){
    b.string = unlist(strsplit(mafft.path, ""))
    if (b.string[length(b.string)] != "/") {
      mafft.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { mafft.path = NULL }

  #Initial checks
  if (is.null(targets.to.align) == T){ stop("A fasta file of targets is needed for alignment.") }

  if (dir.exists(output.directory) == TRUE) {
    if (overwrite == TRUE){
      system(paste0("rm -r ", output.directory))
      dir.create(output.directory)
    }
  } else { dir.create(output.directory) }

  #Load sample file from probe matching step
  all.data = Biostrings::readDNAStringSet(targets.to.align)   # loads up fasta file
  locus.names = unique(gsub("_\\|_.*", "", names(all.data)))

  target.seqs = Biostrings::readDNAStringSet(target.file)   # loads up fasta file

  #Checks for alignments already done and removes from the to-do list
  if (overwrite == FALSE){
    done = list.files(output.directory)
    locus.names = locus.names[!locus.names %in% gsub(".phy$", "", done)]
  }
  if (length(locus.names) == 0){ quit() }

  #Figures out start and end of subset
  sub.start = floor(subset.start * length(locus.names))
  if (sub.start == 0){ sub.start = 1}
  sub.end = floor(subset.end * length(locus.names))

  #Loops through each locus and writes each species to end of file
  for (i in sub.start:sub.end) {

    #Match probe names to contig names to acquire data
    match.data = all.data[grep(pattern = paste0(locus.names[i], "_"), x = names(all.data))]

    #STEP 1: Throw out loci if there are too few taxa
    if (length(names(match.data)) <= min.taxa){
      print(paste0(locus.names[i], " had too few taxa"))
      next
    }

    #STEP 2: Sets up fasta for aligning
    names(match.data) = gsub(pattern = ".*_\\|_", replacement = "", x = names(match.data))

    #Gets reference locus
    ref.locus = target.seqs[grep(pattern = paste(locus.names[i], "$", sep = ""), x = gsub("_\\|_.*", "", names(target.seqs) ) )]
    names(ref.locus) = paste("Reference_Locus")
    final.loci = append(match.data, ref.locus)

    # Checks for duplicates
    dup.names = match.data[duplicated(names(match.data)), ]
    if (length(dup.names) != 0) {
      print(paste0(locus.names[i], " did not successfully align. Duplicate samples, likely paralog."))
      next
    } # end if

    #STEP 3: Runs MAFFT to align
    #Aligns and then reverses back to correction orientation
    alignment = runMafft(sequence.data = final.loci,
                         save.name = paste0(output.directory, "/", locus.names[i]),
                         algorithm = algorithm,
                         adjust.direction = adjust.direction,
                         threads = threads,
                         cleanup.files = T,
                         quiet = quiet,
                         mafft.path = mafft.path)

    #Checks for failed mafft run
    if (length(alignment) == 0){
      print(paste0(locus.names[i], " did not successfully align."))
      next }

    #Removes the reverse name
    if (remove.reverse.tag == TRUE){
      names(alignment) = gsub(pattern = "^_R_", replacement = "", x = names(alignment))
    }#end if


    #Gets the divergence to make sure not crazy
    diff = pairwiseDistanceTarget(alignment, "Reference_Locus")
    bad.seqs = names(diff)[which(diff >= removal.threshold)]
    rem.align = alignment[!names(alignment) %in% bad.seqs]

    # Moves onto next loop in there are no good sequences
    if (length(rem.align) <= as.numeric(min.taxa)){
      #Deletes old files
      print(paste(locus.names[i], " had too few taxa", sep = ""))
      next }

    ### realign if bad seqs removed
    if (length(bad.seqs) != 0){

      alignment = runMafft(sequence.data = rem.align,
                           save.name = paste0(output.directory, "/", locus.names[i]),
                           algorithm = algorithm,
                           adjust.direction = adjust.direction,
                           threads = threads,
                           cleanup.files = T,
                           quiet = quiet,
                           mafft.path = mafft.path)

      #Checks for failed mafft run
      if (length(alignment) == 0){
        print(paste0(locus.names[i], " did not successfully align."))
        next }

      #Removes the reverse name
      if (remove.reverse.tag == TRUE){
        names(alignment) = gsub(pattern = "^_R_", replacement = "", x = names(alignment))
      }#end if

    } # end bad.seqs if

    #Saves prelim exon file
    red.align = alignment[!names(alignment) %in% "Reference_Locus"]
    new.align = strsplit(as.character(red.align), "")
    aligned.set = as.matrix(ape::as.DNAbin(new.align))

    #readies for saving
    writePhylip(alignment = aligned.set,
                file=paste0(output.directory, "/", locus.names[i], ".phy"),
                interleave = F,
                strict = F)

    #Deletes old files
    print(paste0(locus.names[i], " alignment saved."))

  }# end big i loop


}# end function


#END SCRIPT

