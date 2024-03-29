#' @title createContaminantDB
#'
#' @description Function for removing contamination from other organisms from adaptor trimmed Illumina sequence data using BWA
#'
#' @param decontamination.list path to a folder of adaptor trimmed reads in fastq format.
#'
#' @param output.directory the new directory to save the adaptor trimmed sequences
#'
#' @param include.human the new directory to save the adaptor trimmed sequences
#'
#' @param include.mouse directory of genomes contaminants to scan samples
#'
#' @param overwrite system path to samtools in case it can't be found
#'
#' @return a new directory of adaptor trimmed reads and a summary of the trimming in logs/
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

createContaminantDB = function(decontamination.list = NULL,
                               output.directory = "contaminant-references",
                               include.human = TRUE,
                               include.univec = TRUE,
                               include.genbank = NULL,
                               include.fasta = NULL,
                               overwrite = FALSE) {

  # #Debug
  # setwd("/Users/chutter/Dropbox/Research/0_Github/Test-dataset")
  # output.directory = "contaminant-references"
  # decontamination.list = "decontamination_database.csv"
  # overwrite = FALSE
  # quiet = TRUE
  # include.human = TRUE
  # include.univec = TRUE
  # include.genbank = NULL
  # include.fasta = NULL
  # # map.match = 0.99

  #Quick checks
  if (is.null(decontamination.list) == TRUE){ stop("Please provide list of sequences and genbank numbers.") }
  if (file.exists(decontamination.list) == F){ stop("Input list not found.") }

  #Sets directory and reads in  if (is.null(output.dir) == TRUE){ stop("Please provide an output directory.") }
  if (dir.exists(output.directory) == F){ dir.create(output.directory) } else {
    if (overwrite == TRUE){
      system(paste0("rm -r ", output.directory))
      dir.create(output.directory)
    } else { return("Contaminant database already exists, skipping.") }
  }#end else

  sample.data = read.csv(file = decontamination.list)

  if (include.human == FALSE){ sample.data = sample.data[!sample.data$Genome %in% c("Human", "Homo", "Sapien", "Homo_Sapien", "Homo Sapien"),] }

  if (nrow(sample.data) == 0){ stop("no samples remain to analyze.") }

  if (is.null(include.fasta) != TRUE){
    if (file.exists(include.fasta) == F){ stop("include.fasta file not found.") }
    system(paste0("cp ", include.fasta, " ", output.directory, "/manually-included-data.fa"))
  }#end if

  if (is.null(include.genbank) != T){

    for (i in 1:length(include.genbank)){

      biomartr::getGenome(db = "genbank",
                          organism = include.genbank[i],
                          path = output.directory,
                          reference = FALSE)

      new.name = paste0("include-genbank_", include.genbank[i])

      file.list = list.files(output.directory)
      sample.files = file.list[grep(include.genbank[i], file.list)]
      old.name = sample.files[grep(".fna.gz|.fa$", sample.files)]

      system(paste0("mv ", output.directory, "/", old.name, " ", output.directory, "/", new.name))
    }#end i

  }#end is.null

  if (include.univec == TRUE){
   system(paste0("wget -c https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec -O ", output.directory, "/UniVec.fa"))
  }
  sample.data = sample.data[!sample.data$Genome %in% c("univec", "UniVec", "UNIVEC", "Univec"),]


  for (i in 1:nrow(sample.data)){
    biomartr::getGenome(db = "genbank",
                        organism = sample.data$GenBank_Accession[i],
                        path = output.directory,
                        reference = FALSE)

    new.name = paste0(sample.data$Genome[i], "-", sample.data$GenBank_Accession[i], ".fna.gz")

    file.list = list.files(output.directory)
    sample.files = file.list[grep(sample.data$GenBank_Accession[i], file.list)]
    old.name = sample.files[grep(".fna.gz|.fa$", sample.files)]

    system(paste0("mv ", output.directory, "/", old.name, " ", output.directory, "/", new.name))

  }#end i loop

  file.list = list.files(output.directory)
  del.files = file.list[grep(".fna.gz|.fa$", file.list, invert = T)]
  system(paste0("rm ", paste0(output.directory, "/", del.files, collapse = " ")))

  final.files = list.files(output.directory)

  if (length(final.files) != 0){ print("Creation of contamination database was successful.") }

} #end function

