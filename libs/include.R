## Just include all the files in this dir
## Need to have a global variable called
## code.dir before you source this. It's 
## not an ideal solution...

dir <- getwd()
setwd(paste(code.dir ,"/libs/", sep=""))
for( file in list.files() ){
  if(file != "include.R" & substring(file, nchar(file)-1)==".R"){
    source(file)
  }
}
try(dyn.load("inference.so"))            #Try to load the shared lib for density sampling
setwd(dir)

 
