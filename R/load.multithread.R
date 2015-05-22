load.multithread <-function(cores)
{
### Set number of cores default to 1
### MULTITHREADING
### If running on multiple cores, check for packages "foreach" and "doMC" (linux) or "foreach" and "doSNOW" (windows)
    OS<-tolower(.Platform$OS.type)
    if(cores>1)
    {
        # for linux based operating systems:
        if(OS=="unix")
        {
            if(any(installed.packages()[,1]=="doMC")==F) 
            {
                print ("Package 'doMC' is not installed. You can install this with the function install.packages()")
                con<-file("answer.txt")
                cat("Do you want to do it now? (y/n)")
                answer3 <- readLines(,1)
                close(con)
                if(answer3=="y")
                {
                    install.packages("doMC")
                }else
                {
                    stop("Stopped. Package 'doMC' is critical for script to run on multiple processors. Run script with cores=1.")
                }
            }
            library("doMC")                                    
            registerDoMC(cores)                                 
        }
        
        # for windows based operating systems:
        if(OS=="windows")
        {
            if(any(installed.packages()[,1]=="doSNOW")==F) 
            {
                print ("Package 'doSNOW' is not installed. You can install this with the function install.packages()")
                con<-file("answer.txt")
                cat("Do you want to do it now? (y/n)")
                answer3 <- readLines(,1)
                close(con)
                if(answer3=="y")
                {
                    install.packages("doSNOW")
                }else
                {   
                    stop("Stopped. Package 'doSNOW' is critical for script to run on multiple processors. Run script with cores=0.")
                }
            }
        library("doSNOW")                                  
        cluster<-makeCluster(cores,type="SOCK")             
        registerDoSNOW(cluster)                              
        }
    }    
}
