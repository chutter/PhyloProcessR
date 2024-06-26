#' @title summary.missingData
#'
#' @description Function for removing contamination from other organisms from adaptor trimmed Illumina sequence data using BWA
#'
#' @param input.reads path to a folder of adaptor trimmed reads in fastq format.
#'
#' @param output.directory the new directory to save the adaptor trimmed sequences
#'
#' @param decontamination.path directory of genomes contaminants to scan samples
#'
#' @param samtools.path system path to samtools in case it can't be found
#'
#' @param bwa.path system path to bwa in case it can't be found
#'
#' @param threads number of computation processing threads
#'
#' @param mem amount of system memory to use
#'
#' @param resume TRUE to skip samples already completed
#'
#' @param overwrite TRUE to overwrite a folder of samples with output.dir
#'
#' @param quiet TRUE to supress screen output
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

missingData = function(alignment.directory = NULL,
                                      target.file = NULL,
                                      output.directory = "sample-sensitivity",
                                      threads = 1,
                                      memory = 1,
                                      mafft.path = NULL,
                                      overwrite = FALSE,
                                      quiet = FALSE) {


    # alignment.directory = "/Users/chutter/Dropbox/VenomCap_test_data/untrimmed_all-markers"
    # target.file = "/Users/chutter/Dropbox/VenomCap_test_data/venom_loci_updated_Mar12_cdhit95_duplicate_exons_renamed_Feb2023_FINAL.fa"
    # output.directory = "/Volumes/LaCie/VenomCap/data-analysis/genetic-distance"
    # mafft.path = "/Users/chutter/Bioinformatics/miniconda3/envs/PhyloProcessR/bin"
    # quiet = TRUE
    # overwrite = FALSE
    # threads = 6
    # memory = 6

  ####################################################################
  ##### Required program path check
  ####################################################################

  require(foreach)

  #Same adds to bbmap path
  if (is.null(mafft.path) == FALSE){
    b.string = unlist(strsplit(mafft.path, ""))
    if (b.string[length(b.string)] != "/") {
      mafft.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { mafft.path = "" }

  ####################################################################
  ##### Input path check
  ####################################################################

  if (is.null(alignment.directory) == TRUE){ stop("Please provide input directory.") }

  #Sets directory and reads in  if (is.null(output.dir) == TRUE){ stop("Please provide an output directory.") }
  if (file.exists(output.directory) == F){ dir.create(output.directory) } else {
    if (overwrite == TRUE){
      system(paste0("rm -r ", output.directory))
      dir.create(output.directory)
    }
  }#end else

  ####################################################################
  ##### Raw data collection for each alignment and sample
  ####################################################################

  #Gather alignments
  alignment.files = list.files(alignment.directory, full.names = T)
  alignment.files = alignment.files[grep(".phy$", alignment.files)]
  target.markers = Biostrings::readDNAStringSet(target.file)
  target.markers = target.markers[duplicated(names(target.markers)) != T]

  marker.names = gsub(".*\\/", "", alignment.files)
  marker.names = gsub(".phy$", "", marker.names)

  target.markers = target.markers[names(target.markers) %in% marker.names]
  marker.names = marker.names[marker.names %in% names(target.markers)]

  #Finds taxa
  taxa.temp = c()
  for (i in 1:length(alignment.files)){
    align = Biostrings::DNAStringSet(Biostrings::readAAMultipleAlignment(file = alignment.files[i], format = "phylip"))
    taxa.temp = unique(append(taxa.temp, names(align)))
  }

  sample.names = unique(taxa.temp)

  #Sets up multiprocessing
  cl = parallel::makeCluster(threads, outfile = "")
  doParallel::registerDoParallel(cl)
  mem.cl = floor(memory/threads)

  #Loops through each locus and does operations on them
  #all.data = c()
  all.data = foreach::foreach(i=1:length(marker.names),  .combine = rbind, .packages = c("PhyloProcessR", "foreach", "Biostrings","data.table", "ape", "stringr")) %dopar% {
  #Loops through each alignment
  #for (i in 1:length(marker.names)){

    #START HERE
    align.file = alignment.files[grep(paste0(marker.names[i], ".phy$"), alignment.files)]
    align = Biostrings::DNAStringSet(Biostrings::readAAMultipleAlignment(file = align.file, format = "phylip"))

    temp.target = target.markers[names(target.markers) == marker.names[i],]
    temp.target = temp.target[1]
    names(temp.target) = "Reference_Marker"

    #Aligns and then reverses back to correction orientation
    alignment <- PhyloProcessR::runMafft(
      sequence.data = align,
      add.contigs = temp.target,
      algorithm = "add",
      adjust.direction = TRUE,
      threads = 1,
      cleanup.files = T,
      quiet = quiet,
      mafft.path = mafft.path
    )

    #Aligns and then reverses back to correction orientation
    reversed = names(alignment)[grep(pattern = "^_R_", names(alignment))]
    if (length(reversed[grep(pattern = "Reference_Marker", reversed)]) == 1){ alignment = Biostrings::reverseComplement(alignment) }
    names(alignment) = gsub(pattern = "^_R_", replacement = "", x = names(alignment))

    # #Removes the edge gaps
    ref.aligned <- as.character(alignment["Reference_Marker"])
    not.gaps <- stringr::str_locate_all(ref.aligned, pattern = "[^-]")[[1]][, 1]
    ref.start <- min(not.gaps)
    ref.finish <- max(not.gaps)

    # Finds weird gaps to fix
    temp.gaps <- as.numeric(1)
    for (k in 1:length(not.gaps) - 1) {
        temp.gaps <- append(temp.gaps, not.gaps[k + 1] - not.gaps[k])
    }
    temp.gaps <- temp.gaps - 1
    names(temp.gaps) <- not.gaps
    bad.gaps <- which(temp.gaps >= 30)
    front.gaps <- bad.gaps[bad.gaps <= length(not.gaps) * 0.10]
    end.gaps <- bad.gaps[bad.gaps >= length(not.gaps) * 0.90]

    # Fix big gaps if there are any
    if (length(front.gaps) != 0) {
        temp.change <- (max(as.numeric(names(front.gaps)) - ref.start)) - (max(front.gaps) - 1)
        ref.start <- ref.start + temp.change
    } # end gap if

    # Fix big gaps if there are any
    if (length(end.gaps) != 0) {
        add.bp <- length(temp.gaps) - min(end.gaps)
        # add.bp<-(ref.finish-min(as.numeric(names(end.gaps))))
        min.gaps <- temp.gaps[min(end.gaps)]
        temp.change <- as.numeric(names(min.gaps)) - as.numeric(min.gaps)
        ref.finish <- temp.change + add.bp
    } # end gap if

    target.region <- Biostrings::subseq(alignment, ref.start, ref.finish)

    #Removes the edge gaps
    ref.aligned = as.character(target.region['Reference_Marker'])
    not.gaps = stringr::str_locate_all(ref.aligned, pattern = "[^-]")[[1]][,1]
    ref.start = min(not.gaps)
    ref.finish = max(not.gaps)
    trim.align = Biostrings::subseq(target.region, ref.start, ref.finish)

    trim.cols = PhyloProcessR::trimAlignmentColumns(
      alignment = trim.align,
      min.gap.percent = 50
    )

    #Gets lengths
    write.temp = strsplit(as.character(trim.cols), "")
    aligned.set = as.matrix(ape::as.DNAbin(write.temp) )
    len.temp = as.character(as.list(aligned.set))
    len.loci = lapply(len.temp, function (x) x[x != "-"])
    spp.len = unlist(lapply(len.loci, function (x) length(x)))
    spp.sens = spp.len/Biostrings::width(temp.target)
    spp.sens[spp.sens > 1] = 1
    spp.sens = spp.sens[names(spp.sens) != "Reference_Marker"]
    spp.sens = 1 - spp.sens

    #Sets up data to collect
    header.all = c("sample" , "target_marker", "target_length", "missing_basepair_data", "missing_marker_data")
    sample.data = data.table::data.table(matrix(as.numeric(-1), nrow = length(sample.names), ncol = length(header.all)))
    data.table::setnames(sample.data, header.all)
    sample.data[, sample:=as.character(sample.names)]
    sample.data[, target_marker:=as.character(target_marker)]
    sample.data[, target_marker:=as.character(marker.names[i])]
    sample.data[, target_length:=Biostrings::width(temp.target)]

    data.table::set(sample.data, i = match(names(spp.sens), sample.data$sample), j = match("missing_basepair_data", header.all), value = spp.sens)

    data.table::set(sample.data, i = match(names(spp.sens), sample.data$sample), j = match("missing_marker_data", header.all), value = spp.sens)

    sample.data[sample.data$missing_basepair_data == -1]$missing_basepair_data = NA
    sample.data[sample.data$missing_marker_data >= 0]$missing_marker_data = 1
    sample.data[sample.data$missing_marker_data == -1]$missing_marker_data = 0

    #Writes sample table
    #all.data = rbind(all.data, sample.data)
    print(data.frame(sample.data))

  }#end i loop

  parallel::stopCluster(cl)

  #Writes final table
  write.table(all.data, file = paste0(output.directory, "/missing_basepair_data_raw-data.txt"), sep = "\t", row.names = F)

  print(paste0("Saved raw data at ", output.directory, "/missing_basepair_data_raw-data.txt"))

  ####################################################################
  ##### Summarizes the raw data previously collected into a summary table
  ####################################################################

  #Get stats
  temp.min = unlist(lapply( split(all.data$missing_basepair_data, all.data$sample), min, na.rm = TRUE))
  temp.max = unlist(lapply(split(all.data$missing_basepair_data, all.data$sample), max, na.rm = TRUE))
  temp.mean = unlist(lapply(split(all.data$missing_basepair_data, all.data$sample), mean, na.rm = TRUE))
  temp.median = unlist(lapply(split(all.data$missing_basepair_data, all.data$sample), median, na.rm = TRUE))
  temp.sd = unlist(lapply(split(all.data$missing_basepair_data, all.data$sample), sd, na.rm = TRUE))

  #Sets up data to collect
  header.summ = c("sample", "mean_missing_basepair_data", "median_missing_basepair_data",
                 "min_missing_basepair_data", "max_missing_basepair_data", "sd_missing_basepair_data")
  summary.data = data.table::data.table(matrix(as.numeric(0), nrow = length(sample.names), ncol = length(header.summ)))
  data.table::setnames(summary.data, header.summ)
  summary.data[, sample:=as.character(sample)]
  summary.data[, sample:=as.character(sample.names)]

  data.table::set(summary.data, i = match(names(temp.min), summary.data$sample), j = match("min_missing_basepair_data", header.summ), value = temp.min)
  data.table::set(summary.data, i = match(names(temp.max), summary.data$sample), j = match("max_missing_basepair_data", header.summ), value = temp.max)
  data.table::set(summary.data, i = match(names(temp.mean), summary.data$sample), j = match("mean_missing_basepair_data", header.summ), value = temp.mean)
  data.table::set(summary.data, i = match(names(temp.median), summary.data$sample), j = match("median_missing_basepair_data", header.summ), value = temp.median)
  data.table::set(summary.data, i = match(names(temp.sd), summary.data$sample), j = match("sd_missing_basepair_data", header.summ), value = temp.sd)

  write.table(summary.data, file = paste0(output.directory, "/missing_basepair_data_summary.txt"), sep = "\t", row.names = F)

  print(paste0("Saved raw data at ", output.directory, "/missing_basepair_data_summary.txt"))

# Sets up data to collect for missing marker data
##################################################

temp.total <- unlist(lapply(split(all.data$missing_marker_data, all.data$sample), sum, na.rm = TRUE))

header.summ <- c(
    "sample", "total_markers", "missing_markers", "proportion_missing"
)

summary.data <- data.table::data.table(matrix(as.numeric(0), nrow = length(sample.names), ncol = length(header.summ)))
data.table::setnames(summary.data, header.summ)
summary.data[, sample := as.character(sample)]
summary.data[, sample := as.character(sample.names)]

data.table::set(summary.data, i = match(names(temp.min), summary.data$sample), j = match("total_markers", header.summ), value = temp.total)
data.table::set(summary.data, i = match(names(temp.max), summary.data$sample), j = match("missing_markers", header.summ), value = length(alignment.files) - temp.total)
data.table::set(summary.data, i = match(names(temp.mean), summary.data$sample), j = match("proportion_missing", header.summ), value = (length(alignment.files) - temp.total) / length(alignment.files))

write.table(summary.data, file = paste0(output.directory, "/missing_marker_data_summary.txt"), sep = "\t", row.names = F)

print(paste0("Saved raw data at ", output.directory, "/missing_marker_data_summary.txt"))


}#end function


#
#
#     # #Sets up the save data
#     # header.data<-c("ProbeSet", "Scale", "Sample")
#     # loci.headers<-gsub(".phy$", "", locus.names)
#     # new.header<-append(header.data, loci.headers)
#     # collect.data<-data.table(matrix(as.numeric(0), nrow = length(sample.names$Sample), ncol = length(locus.names)+length(header.data) ))
#     # setnames(collect.data, new.header)
#     # collect.data[, ProbeSet:=as.character(sample.names$Probe_Set)]
#     # collect.data[, Scale:=as.character(rep(s.name, times = length(sample.names$Sample)) )]
#     # collect.data[, Sample:=as.character(sample.names$Sample)]
#     #
#     # dist.data<-data.table(matrix(as.numeric(), nrow = length(sample.names$Sample), ncol = length(locus.names)+length(header.data) ))
#     # setnames(dist.data, new.header)
#     # dist.data[, ProbeSet:=as.character(sample.names$Probe_Set)]
#     # dist.data[, Scale:=as.character(rep(s.name, times = length(sample.names$Sample)) )]
#     # dist.data[, Sample:=as.character(sample.names$Sample)]
#
#
#     }#end j loop
#
#     ### Saves all the scale level data now
#     ##########
#     setwd("/Users/chutter/Dropbox/Research/2_WIP/Anura_Seqcap/Paper_statistics")
#     #Saves raw data
#     write.table(collect.data, file = paste0(s.name, "_species_sensitivity_raw.txt"), row.names = F)
#     write.table(dist.data, file = paste0(s.name, "_species_distances_raw.txt"), row.names = F)
#
#     #Saves summary data for each species
#     species.collect<-rowMeans(collect.data[,4:(ncol(collect.data))], na.rm = T)
#     species.dist<-rowMeans(dist.data[,4:(ncol(dist.data))], na.rm = T)
#     species.write<-data.frame(Sample = collect.data$Sample, Sensitivity = species.collect, Distances = species.dist)
#     write.table(species.write, file = paste0(s.name, "_species_sens-dist_summary.txt"), row.names = F)
#
#     #Saves summary data for each locus
#     locus.collect<-colMeans(collect.data[,4:(ncol(collect.data))], na.rm = T)
#     locus.dist<-colMeans(dist.data[,4:(ncol(dist.data))], na.rm = T)
#     locus.write<-data.frame(Locus = colnames(collect.data)[4:ncol(collect.data)],
#                             Sensitivity = locus.collect, Distances = locus.dist)
#     write.table(locus.write, file = paste0(s.name, "_locus_sens-dist_summary.txt"), row.names = F)
#
#     #Collects together summary data per dataset
#     scale.summary$Scale[i]<-s.name
#     scale.summary$S_Mean[i]<-mean(rowMeans(collect.data[,4:(ncol(collect.data))], na.rm = T))
#     scale.summary$S_Median[i]<-median(rowMeans(collect.data[,4:(ncol(collect.data))], na.rm = T))
#     scale.summary$S_Min[i]<-min(rowMeans(collect.data[,4:(ncol(collect.data))], na.rm = T))
#     scale.summary$S_Max[i]<-max(rowMeans(collect.data[,4:(ncol(collect.data))], na.rm = T))
#     scale.summary$S_SD[i]<-sd(rowMeans(collect.data[,4:(ncol(collect.data))], na.rm = T))
#     #Distances
#     scale.summary$D_Mean[i]<-mean(rowMeans(dist.data[,4:(ncol(dist.data))], na.rm = T))
#     scale.summary$D_Median[i]<-median(rowMeans(dist.data[,4:(ncol(dist.data))], na.rm = T))
#     scale.summary$D_Min[i]<-min(rowMeans(dist.data[,4:(ncol(dist.data))], na.rm = T))
#     scale.summary$D_Max[i]<-max(rowMeans(dist.data[,4:(ncol(dist.data))], na.rm = T))
#     scale.summary$D_SD[i]<-sd(rowMeans(dist.data[,4:(ncol(dist.data))], na.rm = T))
#
#   }#end i loop
#
#   setwd("/Users/chutter/Dropbox/Research/2_WIP/Anura_Seqcap/Paper_statistics")
#   write.table(scale.summary, file = "Scale_sensitivity_distance_summary.txt", row.names = F)
#
#
#
#
#
#   ###### OLDDDDDD
#
#
#   #  #Read in sample data **** sample is run twice?!
#   reads = list.files(read.directory, recursive = T, full.names = T)
#   if (is.null(sub.directory) != TRUE) {
#     reads = reads[grep(paste0(sub.directory, "/"), reads)]
#     sample.names = gsub(paste0("/", sub.directory, "/.*"), "", reads)
#     sample.names = unique(gsub(paste0(read.directory, "/"), "", sample.names))
#   } else {
#     sample.names = list.files(read.directory, recursive = F, full.names = F)
#   }
#
#   if (length(sample.names) == 0){ return("no samples remain to analyze.") }
#
#   #Read in group data and the target markers
#   group.data = read.csv(sample.groups, header = T)
#   target.markers = Biostrings::readDNAStringSet(target.file)
#
#   #Specificity refers to the percentage of cleaned reads that can be mapped back to the target markers.
#
#   #Sets up data to collect
#   header.all = c("group","sample", "mapped_reads", "unmapped_reads","unmapped_mates",
#                   "total_reads", "percent_mapped", "median_marker_rpkm", "mean_marker_rpkm")
#   all.data = data.table::data.table(matrix(as.numeric(0), nrow = length(sample.names), ncol = length(header.all)))
#   data.table::setnames(all.data, header.all)
#   all.data[, sample:=as.character(sample)]
#   all.data[, group:=as.character(group)]
#
#   #Loops through each locus and does operations on them
#   for (i in 1:length(sample.names)){
#
#     #Reads in the contig file to check
#     dir.create(paste0(output.directory, "/", sample.names[i]))
#
#     #Sets up data to collect
#     header.data = c("group","sample", "marker", "mapped_reads", "unmapped_mates", "total_read_pairs", "percent_mapped", "marker_rpkm")
#     collect.data = data.table::data.table(matrix(as.numeric(0), nrow = length(target.markers), ncol = length(header.data)))
#     data.table::setnames(collect.data, header.data)
#     collect.data[, sample:=as.character(sample)]
#     collect.data[, group:=as.character(group)]
#     collect.data[, marker:=as.character(marker)]
#
#     #################################################
#     ### Part A: prepare for loading and checks
#     #################################################
#     sample.reads = reads[grep(pattern = paste0(sample.names[i], "_"), x = reads)]
#
#     #Checks the Sample column in case already renamed
#     if (length(sample.reads) == 0){ sample.reads = reads[grep(pattern = sample.names[i], x = reads)] }
#
#     sample.reads = unique(gsub("_1.f.*|_2.f.*|_3.f.*|-1.f.*|-2.f.*|-3.f.*|_R1_.*|_R2_.*|_R3_.*|_READ1_.*|_READ2_.*|_READ3_.*|_R1.f.*|_R2.f.*|_R3.f.*|-R1.f.*|-R2.f.*|-R3.f.*|_READ1.f.*|_READ2.f.*|_READ3.f.*|-READ1.f.*|-READ2.f.*|-READ3.f.*|_singleton.*|-singleton.*|READ-singleton.*|READ_singleton.*|_READ-singleton.*|-READ_singleton.*|-READ-singleton.*|_READ_singleton.*", "", sample.reads))
#
#     #Returns an error if reads are not found
#     if (length(sample.reads) == 0 ){
#       stop(sample.names[i], " does not have any reads present for files ")
#     } #end if statement
#
#     for (j in 1:length(sample.reads)){
#
#       lane.reads = reads[grep(pattern = paste0(sample.reads[j], "_"), x = reads)]
#
#       #Checks the Sample column in case already renamed
#       if (length(lane.reads) == 0){ lane.reads = reads[grep(pattern = sample.reads[j], x = reads)] }
#       #Returns an error if reads are not found
#       if (length(lane.reads) == 0 ){
#         stop(sample.reads[j], " does not have any reads present for files ")
#       } #end if statement
#
#       #Copies ortholog file to act as a reference
#       system(paste0("cp ", target.file, " ", output.directory, "/", sample.names[i], "/sample-reference.fa"))
#       system(paste0(bwa.path, "bwa index ", output.directory, "/", sample.names[i], "/sample-reference.fa"),
#              ignore.stdout = quiet, ignore.stderr = quiet)
#
#       #Creates a bam alignment file of reads mapped to reference
#       system(paste0(bwa.path, "bwa mem -M -t ", threads, " ",
#                     output.directory, "/", sample.names[i], "/sample-reference.fa ",
#                     lane.reads[1], " ", lane.reads[2],
#                     " | ", samtools.path, "samtools sort -@", threads, " -O BAM",
#                     " -o ", output.directory, "/", sample.names[i], "/paired.bam  -"),
#              ignore.stdout = quiet, ignore.stderr = quiet)
#
#       system(paste0(picard.path, "picard -Xmx", memory, "g",
#                     " MarkDuplicates INPUT=", output.directory, "/", sample.names[i], "/paired.bam",
#                     " OUTPUT=", output.directory, "/", sample.names[i], "/dedup_paired.bam",
#                     " METRICS_FILE=", output.directory, "/", sample.names[i], "/metrics.txt",
#                     " USE_JDK_DEFLATER=true USE_JDK_INFLATER=true"),
#              ignore.stdout = quiet, ignore.stderr = quiet)
#
#       system(paste0("mv ", output.directory, "/", sample.names[i], "/metrics.txt ",
#                     output.directory, "/", sample.names[i], "/duplication_stats.txt"))
#
#       #HERe again
#       system(paste0(picard.path, "picard -Xmx", memory, "g",
#                     " BuildBamIndex INPUT=", output.directory, "/", sample.names[i], "/dedup_paired.bam",
#                     " USE_JDK_DEFLATER=true USE_JDK_INFLATER=true"),
#              ignore.stdout = quiet, ignore.stderr = quiet)
#
#       #Gets the genome coverage using bedtools
#       system(paste0(samtools.path, "samtools idxstats --threads ", threads, " ",
#                     output.directory, "/", sample.names[i], "/dedup_paired.bam > ",
#                     output.directory, "/", sample.names[i], "/samtools_idxstats.txt"))
#
#       system(paste0("rm ", output.directory, "/", sample.names[i], "/paired.bam"))
#
#       #IDX data
#       i.headers = c("marker", "length", "mapped_reads", "unmapped_mates")
#       idx.stats = data.table::fread(file = paste0(output.directory, "/", sample.names[i], "/samtools_idxstats.txt"))
#       data.table::setnames(idx.stats, i.headers)
#       idx.stats = idx.stats[idx.stats$marker != "*",]
#
#       #Gets the number of reads that mapped to the reference, also total reads
#       #Total number of un-mapped reads
#       mapped.all = sum(idx.stats$mapped_reads)
#       unmapped.all = as.numeric(system(paste0("samtools view -c -f 4 ",
#                                               output.directory, "/", sample.names[i], "/dedup_paired.bam"), intern = T))
#
#       #Goes through each locus and calculates stats
#       locus.names = unique(idx.stats$marker)
#       group.name = group.data[group.data[,2] %in% sample.names[i],][,1]
#       locus.rpkm = idx.stats$mapped_reads / (idx.stats$length/1000 * mapped.all/1000000)
#
#       #Starter stats
#       data.table::set(collect.data, j = match("group", header.data), value = group.name )
#       data.table::set(collect.data, j = match("sample", header.data), value = sample.names[i] )
#       data.table::set(collect.data, j = match("marker", header.data), value = locus.names )
#       data.table::set(collect.data, i = match(locus.names, idx.stats$marker), j = match("mapped_reads", header.data), value = idx.stats$mapped_reads )
#       data.table::set(collect.data, i = match(locus.names, idx.stats$marker), j = match("unmapped_mates", header.data), value = idx.stats$unmapped_mates )
#       data.table::set(collect.data, i = match(locus.names, idx.stats$marker), j = match("total_read_pairs", header.data), value = idx.stats$mapped_reads+idx.stats$unmapped_mates )
#       data.table::set(collect.data, i = match(locus.names, idx.stats$marker), j = match("percent_mapped", header.data), value = idx.stats$mapped_reads/(idx.stats$mapped_reads+idx.stats$unmapped_mates) )
#       data.table::set(collect.data, i = match(locus.names, idx.stats$marker), j = match("marker_rpkm", header.data), value = locus.rpkm )
#
#     } #end j loop
#
#     write.table(collect.data, file = paste0(output.directory, "/", sample.names[i], "/species_specificity_data.txt"), sep = "\t", row.names = F)
#
#     #Collects the overall summary data
#     data.table::set(all.data, i = as.integer(i), j = match("group", header.all), value = group.name )
#     data.table::set(all.data, i = as.integer(i), j = match("sample", header.all), value = sample.names[i] )
#     data.table::set(all.data, i = as.integer(i), j = match("mapped_reads", header.all), value = sum(collect.data$mapped_reads) )
#     data.table::set(all.data, i = as.integer(i), j = match("unmapped_reads", header.all), value = unmapped.all )
#     data.table::set(all.data, i = as.integer(i), j = match("unmapped_mates", header.all), value = sum(collect.data$unmapped_mates) )
#     data.table::set(all.data, i = as.integer(i), j = match("total_reads", header.all), value = sum(collect.data$mapped_reads)+unmapped.all )
#     data.table::set(all.data, i = as.integer(i), j = match("percent_mapped", header.all), value = sum(collect.data$mapped_reads)/(sum(collect.data$mapped_reads)+unmapped.all) )
#     data.table::set(all.data, i = as.integer(i), j = match("mean_marker_rpkm", header.all), value = mean(collect.data$marker_rpkm) )
#     data.table::set(all.data, i = as.integer(i), j = match("median_marker_rpkm", header.all), value = median(collect.data$marker_rpkm) )
#
#     print(paste0(sample.names[i], " Completed!"))
#
#   }# end i loop
#
#   write.table(all.data, file = paste0(output.directory, "/sample-specificity_summmary.txt"), sep = "\t", row.names = F)
#
#   #### Print some textual summary here
#
# }# end function

### END SCRIPT
