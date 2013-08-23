## Just include all the files in this dir
## Need to have a global variable called
## code.dir before you source this. It's 
## not an ideal solution...

dir <- getwd()
setwd(paste(code.dir ,"/libs/", sep=""))
for( file in list.files() ){
  if(file != "include.R"){
    source(file)
  }
}
setwd(dir)

 
