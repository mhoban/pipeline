#!/usr/bin/env littler

sampleFastq <- function(infile, outfile, size=opt$sample_size, seed=F) {
  if (seed != F) {
    set.seed(ifelse(seed==T,123,seed)) 
    #this allows us to sample the same locations from multiple files
    #useful for paired-end reads
    #do we need this actually?
  }
  
  fq <- FastqSampler(infile,size) #create the sampler object
  print(sprintf("Subsampling %s...",infile))
  pct <- proc.time() #time this, cuz it's slow
  
  fq.sub <- yield(fq)       #grab our samples
  writeFastq(fq.sub,outfile)   #write 'em to the output
  
  tt <- proc.time() - pct
  print(sprintf("Done in %fs. Saved as %s",tt['elapsed'],outfile))
  close(fq)
}

# suppressMessages(require('docopt'))
# 
# doc <- "Usage: sample_fastq.r [options]
# 
# -d DIR --input=DIR            Input directory [default: ./]
# -p PATTERN --pattern=PATTERN  Match input files to regex [default: *.fastq.gz]
# -o OUTPUT --output=OUTPUT     Directory in which to place subsampled fastq files [default: ./]
# -n SIZE --size=SIZE           Number of records to sample from source fastq files [default: 10000]"
# 
# opt <- docopt(doc)

print("Loading required packages (takes forever)...")
suppressMessages(require('ShortRead')) #ShortRead contains the fastq sampling routine

#canned options
opt <- list(dir="~/code/pipeline/data/shark/fastq/",search_pattern="*R1.*fastq.gz$",output_dir="~/code/pipeline/data/shark/fastq/subsample/",sample_size=10000)


for (f in list.files(path=opt$dir, pattern=opt$search_pattern)) {
  outf <- sub("//","/",sprintf("%s/sub%d_%s",opt$output_dir,opt$sample_size,f)) #build our output filename
  
  #do the forward read (R1)
  sampleFastq(f,outf,seed=123)
  
  #do the reverse read (R2)
  f <- sub("R1","R2",f)
  outf <- sub("//","/",sprintf("%s/sub%d_%s",opt$output_dir,opt$sample_size,f)) #rebuild our output filename
  sampleFastq(f,outf,seed=123)
}

