#' @title dropboxDownload
#'
#' @description Function for downloading data off your dropbox account
#'
#' @param sample.spreadsheet spreadsheet path or table of file names to download
#'
#' @param dropbox.directory the directory that contains these target files
#'
#' @param dropbox.token your dropbox token path, created the first time you run the function but needs interactivity for log-in. Can be done once on your computer and the token can be moved around.
#'
#' @param output.directory the directory to save the files you are downloading
#'
#' @param skip.not.found the directory to save the files you are downloading
#'
#' @param overwrite whether to overwrite or not
#'
#' @return saves dropbox files to output directory
#'
#' @examples
#'
#' your.tree = ape::read.tree(file = "file-path-to-tree.tre")
#' astral.data = astralPlane(astral.tree = your.tree,
#'                           outgroups = c("species_one", "species_two"),
#'                           tip.length = 1)
#'
#' @export

dropboxDownload = function(sample.spreadsheet = NULL,
                           dropbox.directory = NULL,
                           dropbox.token = NULL,
                           output.directory = NULL,
                           skip.not.found = FALSE,
                           overwrite = FALSE){

  # # ### Example usage
  # sample.spreadsheet = "/Users/chutter/Dropbox/Research/2_Master-Databases/seqcap_datasets/1_FrogCap_Read-Database_June12-2022.csv"
  # output.directory = "/Users/chutter/Dropbox/Research/3_Sequence-Database/Raw-Reads-Published/Hutter_Esselstyn_2021/raw-reads-frogs"
  # dropbox.directory = "/Research/3_Sequence-Database/Raw-Reads"
  # dropbox.token = "/home/c111h652/dropbox-token.RDS"
  # dropbox.token = "/Users/chutter/Dropbox/dropbox-token.RDS"
  # overwrite = FALSE

  #Sets directory and reads in  if (is.null(output.dir) == TRUE){ stop("Please provide an output directory.") }
  if (dir.exists(output.directory) == F){ dir.create(output.directory) } else {
    if (overwrite == TRUE){
      system(paste0("rm -r ", output.directory))
      dir.create(output.directory)
    }
  }#end else

  all.reads = rdrop2::drop_dir(path = dropbox.directory,
                               recursive = TRUE)$path_display

  all.reads = all.reads[grep("fastq.gz$|fq.gz$", all.reads)]

  sample.data = read.csv(sample.spreadsheet)

  #Resumes file download
  if (overwrite == FALSE){
    done.files = list.files(output.directory)
    done.files = gsub("_L0.*", "", done.files)
    done.files = unique(done.files)
    sample.data = sample.data[!sample.data$Sample %in% done.files,]
  }

  if (nrow(sample.data) == 0){ return("no samples remain to analyze.") }

  save.miss = c()
  for (i in 1:nrow(sample.data)){

    sample.reads = all.reads[grep(as.character(sample.data$File[i]), all.reads)]

    #For checking if reads are present
    if (length(sample.reads) >= 3){ stop("Problem with reads, row ", i, " File column matches to more than 1 sample. Check to ensure sample spreadsheet has multiple entries for samples with more than 1 lane of data.") }

    #Skip not found or crash
    if (length(sample.reads) == 0) {
      if (skip.not.found == FALSE){
        stop(paste0("Error: sample reads for ", sample.data$Sample[i], " not found!"))
      } else { next }
    }#end if

    if (length(sample.reads) != 2) {
      if (skip.not.found == FALSE){
        stop(paste0("Error: too many sample reads for ", sample.data$Sample[i], " found!"))
      } else { next }
    }#end if

    #Save the read files with the new names in the new directory
    rdrop2::drop_download(path = sample.reads[1],
                          local_path = paste0(output.directory, "/", sample.data$Sample[i], "_L001_READ1.fastq.gz"),
                          overwrite = TRUE)

    rdrop2::drop_download(path = sample.reads[2],
                          local_path = paste0(output.directory, "/", sample.data$Sample[i], "_L001_READ2.fastq.gz"),
                          overwrite = TRUE)
  }#end i loop

}#end function
