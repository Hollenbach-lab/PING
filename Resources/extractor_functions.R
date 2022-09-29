
extractor.bowtie2_align <- function(bowtie2_command, threads, currentSample, extractedFastqDirectory){
  currentSample$samPath <- file.path(extractedFastqDirectory,paste0(currentSample$name,'.sam'))
  optionsCommand <- c('-x Resources/extractor_resources/reference/output',
                     '-5 3','-3 7','-L 20','-i S,1,0.5','--score-min L,0,-0.187',
                     '-I 75','-X 1000', '--no-unal',
                     paste0('-p ',threads),
                     paste0('-1 ',currentSample$rawfastq1path),
                     paste0('-2 ',currentSample$rawfastq2path),
                     paste0('-S ',currentSample$samPath),
                     paste0('--al-conc-gz ',
                            file.path(extractedFastqDirectory,paste0(currentSample$name,'_KIR_%.fastq.gz'))
                            )
                     )
  cat('\n\n',bowtie2_command, optionsCommand)
  output.sampleAlign <- system2(bowtie2_command, optionsCommand, stdout=T, stderr=T)
  
  if(!is.null(attributes(output.sampleAlign))){
    cat('\nBowtie2 alignment failed, retrying...')
    output.sampleAlign <- system2(bowtie2_command, optionsCommand, stdout=T, stderr=T)
  }
  
  if(!is.null(attributes(output.sampleAlign))){
    cat('\nBowtie2 alignment failed, marking sample as failed.')
    currentSample$failed <- T
    currentSample[['failureMessage']] <- paste0(currentSample$name,': Failed KIR extraction.')
    cat(currentSample[['failureMessage']], file = failureLog, sep = "\n", append = TRUE)
  }
  
  check.system2_output(output.sampleAlign, 'Bowtie2 KIR extraction alignment failed.')
  
  cat('\n',paste0(output.sampleAlign, collapse='\n'))
  
  currentSample$kirfastq1path <- normalizePath(file.path(extractedFastqDirectory,paste0(currentSample$name,'_KIR_1.fastq.gz')))
  currentSample$kirfastq2path <- normalizePath(file.path(extractedFastqDirectory,paste0(currentSample$name,'_KIR_2.fastq.gz')))

  cat('\n\nSuccessfully extracted KIR reads for',currentSample$name)
  
  cat('\nCleaning up alignment files.')
  #file.remove('delete.me')
  file.remove(currentSample$samPath)
  
  return(currentSample)
}

extractor.run <- function(sampleList, threads, extractedFastqDirectory, forceRun=F){
  
  # Create the directory if it does not yet exist
  if(!file.exists(extractedFastqDirectory)){
    dir.create(extractedFastqDirectory,recursive=T)
    cat('\n',extractedFastqDirectory,'created.\n')
  }else{
    cat('\n',extractedFastqDirectory,'found.\n')
  }
  
  for(currentSample in sampleList){
    cat('\n\nProcessing',currentSample$name,'----------')
    
    ## forceRun=T forces the alignment to run, even if there are previous results
    if(!forceRun){
      previousResultVect <- grep(currentSample$name, list.files(extractedFastqDirectory,full.names=T),value=T)
      
      ## If there are two files found and they both have _KIR_ in them, 
      #  add the file names to the sample object and skip alignment
      if(length(previousResultVect) == 2 & length(grep('_KIR_',previousResultVect,fixed=T)) == 2){
        cat('\nExtracted files found, skipping alignment.')
        currentSample$kirfastq1path <- grep('KIR_1.fastq.gz',previousResultVect,value=T)
        currentSample$kirfastq2path <- grep('KIR_2.fastq.gz',previousResultVect,value=T)
        
        next
      }
    }
    
    currentSample <- extractor.bowtie2_align(bowtie2, threads, currentSample, extractedFastqDirectory)
  }
  
  cat("\n\n----- PING2_extractor is complete. Extracted reads are deposited in",extractedFastqDirectory,'-----\n')
  
  return(sampleList)
}
