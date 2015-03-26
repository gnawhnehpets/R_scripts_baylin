xdat <- savedx
length(xdat)

xdat <-  c(7,40,41,42,43,44,45,46,47,48,49,50,50.5, 1,2,3,51,52,8,6,53,54,55,57,68,44,47,48,49,50,1,2,3,51,52,8,6,53,54,55,57,68,44,45,46,47,48,49,50,1,2,3,51,52,8,6,53,54,55,57,68,44,45,46,47,48,49,50,1,2,3,51,52,8,6,53,54,44,45,46,47,48,49,50,1,2,3,51,52,8,6,53,54,55,57,68,44,47,48,49,50,1,2,3,51,52,8,6,53,54,55,57,68,44,45,46,47,48,49,50,1,2,3,51,52,8,6,53,54,55,57,68,44,45,46,47,48,49,50,1,2,3,51,52,8,6,53,54,44,45,46,47,48,49,50,1,2,3,51,52,8,6,53,54,55,57,68,44,47,48,49,50,1,2,3,51,52,8,6,53,54,55,57,68,44,45,46,47,48,49,50,1,2,3,51,52,8,6,53,54,55,57,68,44,45,46,47,48,49,50,1,2,3,51,52,8,6,53,54)
xdat <-  rnorm(1000000, mean=0, sd=2)
xdat <-  as.integer(x)
# save(xdat, file="./Desktop/xdata100k.rda")
write.table(xdat, file="./Desktop/xdat.txt", quote=FALSE, row.names=TRUE, col.names=FALSE, sep="\t")

# upperCutoff=10
# windowSize=9
# abc <-  aboveCutoff(xdat, upperCutoff, windowSize)

# abc <-  aboveCutoff(xdat, 1.5, 100)
# def <-  belowCutoff(xdat, -1.5, 100)

abc
def

# aboveCutoff <-  function(xd, y, z){
# #      x <-  xd
# #      cutoff <-  y
# #      window <-  z
     
     x <-  xdat
     cutoff <-  1
     window <-  5
     window <- window-1
     position=1
     i=1
     start=NULL
     end=NULL
     counter=1
     storage <-  matrix(,nrow=0, ncol=4)
     colnames(storage) <-  c("start", "end", "length", "checkpoint")
     
     while(i<length(x)){
          print(paste0("ROUNDi: ", i))
          start.save <-  i
          end <-  i+window
          if(!is.na(x[end])){
               counter=0
               flag=0
               #positive start
               if(x[i]>cutoff){ 
                    
                    # positive end
                    if(x[end]>cutoff){      # PATTERN: POSITIVE START, POSITIVE END
                         i <-  i+1 # if start and end positive, then check in between
                         
                         print(paste0(i, " - ", x[i], " :POSITIVE start; ", end, " - ", x[end], " :POSITIVE end"))
                         print(paste0("end: ", end))
                         
                         # check each position between start&end
                         while(i<end){ 
                              # if position is positive, move to next position
                              if(x[i]>cutoff){ 
                                   counter=counter+1
                                   print(paste0("i: ", i))
                                   print(paste0("counter: ", counter))
                                   # IF ALL POSITION IN BETWEEN START/END ARE POSITIVE
                                   if(counter>=window&&flag==0){
                                        #                          if(flag==0){
                                        print(paste0("HIT: ", counter))
                                        
                                        print(paste0(i, ":while loop done - ", x[i]))
                                        # check to see if points after end position are positive
                                        start <-  start.save
                                        while(x[end]>cutoff){ # THEN CHECK TO SEE HOW MANY POSITIONS AFTER END POSITION ARE POSITIVE
                                             end <-  end+1
                                             print(paste0(i, ":end extension - ", x[end]))
                                        }
                                        end <-  end-1
                                        print(paste0("REGION FOUND2. START: ", start, "  END: ", end))
                                        storage <-  rbind(storage, cbind(start, end, end-start+1, 2))
                                        print(paste0("add: ", start, "; end: ", end))
                                        #                        storage <-  rbind(storage, cbind(start, end, end-start))
                                        print(paste0("NEW START POSITION: ", i))
                                        
                                   }
                                   i <-  i+1
                              }
                              # if position is negative, find start of positive region from the end position, then extend past end position to find end of region
                              else{ 
                                   counter=0 #RESET COUNTER
                                   flag=1
                                   #go upstream of end position until position that doesn't meet cutoff is met
                                   b <-  end-1
                                   while(x[b]>cutoff){ 
                                        b <-  b-1
                                   }
                                   start <-  b+1
                                   #go downstream of end position until position that doesn't meet cutoff is met
                                   c <-  end+1
                                   while(x[c]>cutoff){ 
                                        c <-  c+1
                                   }
                                   end <-  c-1
                                   
                                   #If region is greater than window, then save start and end positions
                                   print(paste0("end: ", end))
                                   if((end-start)>=window-1){ 
                                        storage <-  rbind(storage, cbind(start, end, end-start+1, 1))
                                        print(paste0("REGION FOUND1. START: ", start, "  END: ", end))
                                        print(paste0("add: ", start, "; end: ", end))
                                        break
                                   }
                                   i <-  i+1
                              }
                         
                         }
                         # BUG ######################
                         
                         
                         
                    }
                    # negative end
                    else{ # IF END POSITION IS NEGATIVE, THEN NEW START END+1
                         print(paste0(i, " - ", x[i], " :end position negative: ", x[end]))
                         print("BREAK!")
                    }
               } 
               
               
# negative start
               else{ 
                    start.save <-  i
                    end <-  i+window    
                    # positive end
                    if(x[end]>cutoff){
                         print(paste0(i, " - ", x[i], " :NEGATIVE start; ", end, " - ", x[end], " :POSITIVE end"))
                         j <-  end-1
                         while(x[j]>cutoff){ # FIND START POSITION OF REGION
                              j<- j-1
                         }
                         start <-  j+1
                         
                         k <-  end+1
                         while(x[k]>cutoff){ # FIND END POSITION OF REGION
                              k<- k+1
                         }
                         end <-  k-1
                         if((end-start)>=window-1){ #If region is greater than window, then save start and end positions
                              print(paste0("REGION FOUND3. START: ", start, "  END: ", end))
                              storage <-  rbind(storage, cbind(start, end, end-start+1 , 3))
                              print(paste0("NEW START POSITION: ", i))
                         }
                         #IF END IS POSITIVE, THEN POTENTIAL START POSITION; CHECK EVERY POSITION BEFORE END POSITION FOR START OF REGION
                    }
                    # negative end
                    else{ # IF END IS NEGATIVE, THEN NEW START POSITION AFTER END
                         print(paste0(i, " - ", x[i], ":NEGATIVE start; ", end, " - ", x[end], " :NEGATIVE end"))
                    }
                    
               } 
               i <-  end+1 # NEW START POSITION
               print(paste0("THIS IS NEW START: ", i))
          }
          else{ #IF NO END, RUN OUT INCREMENTS
               i<- i+1
          }
     }
     # proc.time()
     rownames(storage) <-  seq(1:nrow(storage))
     storage

#      return(storage)
# }

storage1 <- storage[storage[,3]>window,]
indices <- NULL
for(i in 1:nrow(storage1)){
     indices <- c(indices, seq(storage1[i,1], storage1[i,2]))
}
which(indices==283131)
length(which(storage1[,3]<window))
length(indices)
length(unique(indices))
length(indices)==sum(storage1[,3])
length(which(x[indices]>0))
length(which(x[indices]>1))
length(which(x[indices]>1.01))
rownames(storage1) <- seq(1:nrow(storage1))
storage1
