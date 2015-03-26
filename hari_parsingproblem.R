http://stackoverflow.com/questions/28284946/how-can-i-skip-increments-in-r-for-loop

x <- rnorm(100000000, mean=0, sd=4) #100million
x <- rnorm(100, mean=10, sd=10)
x <- c(40,41,42,43,44,45,46,47,48,49,50,50.5, 1,2,3,51,52,8,6,53,54,55,57,68,44,47,48,49,50,1,2,3,51,52,8,6,53,54,55,57,68,44,45,46,47,48,49,50,1,2,3,51,52,8,6,53,54,55,57,68,44,45,46,47,48,49,50,1,2,3,51,52,8,6,53,54,44,45,46,47,48,49,50,1,2,3,51,52,8,6,53,54,55,57,68,44,47,48,49,50,1,2,3,51,52,8,6,53,54,55,57,68,44,45,46,47,48,49,50,1,2,3,51,52,8,6,53,54,55,57,68,44,45,46,47,48,49,50,1,2,3,51,52,8,6,53,54,44,45,46,47,48,49,50,1,2,3,51,52,8,6,53,54,55,57,68,44,47,48,49,50,1,2,3,51,52,8,6,53,54,55,57,68,44,45,46,47,48,49,50,1,2,3,51,52,8,6,53,54,55,57,68,44,45,46,47,48,49,50,1,2,3,51,52,8,6,53,54)
# data.frame(x)
# length(which(x>10))
# proc.time()
# length(which(x<10))
# proc.time()


cutoff=0
window=10
position=1
i=1
start=NULL
end=NULL
counter=1
storage <- matrix(,nrow=0, ncol=4)
colnames(storage) <- c("start", "end", "length", "checkpoint")
above <- aboveCutoff(obj)
# aboveCutoff <- function(object){
#      x <- object
     while(i<length(x)){
          print(paste0("ROUNDi: ", i))
        #IF START IS POSITIVE
        if(x[i]>cutoff){ 
             start.save <- i
             end <- i+window
             # AND END POSITION IS POSITIVE
             # PATTERN: POSITIVE START, POSITIVE END
             if(x[end]>cutoff){
                  print(paste0(i, " - ", x[i], " :postive start; ", end, " - ", x[end], " :positive end"))
                  i <- i+1
                  
                  
                  
                  while(i<end){ # CHECK EACH POSITION IN BETWEEN
                       # IF POSITION IS POSITIVE, THEN CHECK NEXT POSITION
                       if(x[i]>cutoff){ 
                            print(paste0(i, ":position positive - ", x[i]))
                            i <- i+1
                            }
                       # IF POSITION IS NEGATIVE, THEN CHECK POSITION BEFORE END POSITION FOR START OF REGION
                       else{ 
                            print(paste0(i, ":position negative - ", x[i]))
                            b <- end-1
                            while(x[b]>cutoff){
                                 b <- b-1
                                 print(paste0(b, ":ADD position positive1 - ", x[b]))
                            }
                            start <- b+1
                         
                            c <- end
                            while(x[c]>cutoff){
                                 c <- c+1
                                 print(paste0(c, ":ADD position positive2 - ", x[c]))
                            }
                            end <- c-1
                            if((end-start)>window-1){ #If region is greater than window, then save start and end positions
                                 print(paste0("REGION FOUND1. START: ", start, "  END: ", end))
                                 storage <- rbind(storage, cbind(start, end, end-start, 1))
                                 i <- end+1
                                 print(paste0("NEW START POSITION: ", i))
                            }else{
                                 i <- i+1
                                 break;
                            }
                       }
                  }
                  # IF ALL POSITION IN BETWEEN START/END ARE POSITIVE
                  if(i==end){ 
                       print(paste0(i, ":while loop done - ", x[i]))
                       # check to see if points after end position are positive
                       start <- start.save
                       while(x[end]>cutoff){ # THEN CHECK TO SEE HOW MANY POSITIONS AFTER END POSITION ARE POSITIVE
                            end <- end+1
                            print(paste0(i, " - ", x[i], ":end extension - ", end, " - ", x[end]))
                       }
                       end <- end-1
                       print(paste0("REGION FOUND2. START: ", start, "  END: ", end))
                       storage <- rbind(storage, cbind(start, end, end-start, 2))
#                        storage <- rbind(storage, cbind(start, end, end-start))
                       i <- end+1
                       print(paste0("NEW START POSITION: ", i))
                       #SAVE start position and check positions after end position for positives
                  }
                  #IF ALL POSITIONS IN BETWEEN START/END ARE NOT POSITIVE
                  else{ 
                       print(paste0(i, ":while loop broken - ", x[i]))
                       i <- end+1 # SET NEW START POSITION
                  }
             }
             else{ # IF END POSITION IS NEGATIVE, THEN NEW START END+1
                  print(paste0(i, "-", x[i], ":end position negative: ", end, " - ", x[end]))
                  i <- end
             }
        }
        #IF START IS NEGATIVE THEN CHECK IF END POSITION IS POSITIVE
        # PATTERN: NEGATIVE START, ? END
        else{
             start.save <- i
             end <- i+window
             if(x[end]>cutoff){
                  print(paste0(i, " - ", x[i], " :negative start; ", end, " - ", x[end], " :positive end"))
                  j <- end
#                   j <- end-1
                  while(x[j]>cutoff){ # FIND START POSITION OF REGION
                       j<-j-1
                  }
                  start <- j+1

                  k <- end
#                   k <- end+1
                  while(x[k]>cutoff){ # FIND END POSITION OF REGION
                       k<-k+1
                  }
                  end <- k+1
                  if((end-start)>window-1){ #If region is greater than window, then save start and end positions
                       print(paste0("REGION FOUND3. START: ", start, "  END: ", end))
                       storage <- rbind(storage, cbind(start, end, end-start, 3))
                       i <- end+1
                       print(paste0("NEW START POSITION: ", i))
                  }else{
                       i <- end+1
                  }
                  #IF END IS POSITIVE, THEN POTENTIAL START POSITION; CHECK EVERY POSITION BEFORE END POSITION FOR START OF REGION
             }else{ # IF END IS NEGATIVE, THEN NEW START POSITION AFTER END
                  print(paste0(i, " :negative start; ", x[end], " :negative end"))
                  i <- end+1
             }
        }
      }
     # proc.time()
     rownames(storage) <- seq(1:nrow(storage))
     storage
#      return storage
# }



