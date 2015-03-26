x <- c(7,40,41,42,43,44,45,46,47,48,49,50,50.5, 1,2,3,51,52,8,6,53,54,55,57,68,44,47,48,49,50,1,2,3,51,52,8,6,53,54,55,57,68,44,45,46,47,48,49,50,1,2,3,51,52,8,6,53,54,55,57,68,44,45,46,47,48,49,50,1,2,3,51,52,8,6,53,54,44,45,46,47,48,49,50,1,2,3,51,52,8,6,53,54,55,57,68,44,47,48,49,50,1,2,3,51,52,8,6,53,54,55,57,68,44,45,46,47,48,49,50,1,2,3,51,52,8,6,53,54,55,57,68,44,45,46,47,48,49,50,1,2,3,51,52,8,6,53,54,44,45,46,47,48,49,50,1,2,3,51,52,8,6,53,54,55,57,68,44,47,48,49,50,1,2,3,51,52,8,6,53,54,55,57,68,44,45,46,47,48,49,50,1,2,3,51,52,8,6,53,54,55,57,68,44,45,46,47,48,49,50,1,2,3,51,52,8,6,53,54)
x <- rnorm(100, mean=10, sd=10)
x <- as.integer(x)

cutoff=10
window=9
position=1



i=1
start=NULL
end=NULL
counter=1
storage <- matrix(,nrow=0, ncol=4)
colnames(storage) <- c("start", "end", "length", "checkpoint")
above <- aboveCutoff(obj)
while(i<length(x)){
     print(paste0("ROUNDi: ", i))
     start.save <- i
     end <- i+window
     if(x[i]>cutoff){ 
          if(x[end]>cutoff){      # PATTERN: POSITIVE START, POSITIVE END
               i <- i+1
               counter=0
               print(paste0(i, " - ", x[i], " :POSITIVE start; ", end, " - ", x[end], " :POSITIVE end"))
               print(paste0("end: ", end))
               while(i<end){ # CHECK EACH POSITION IN BETWEEN
                    
                    if(x[i]>cutoff){ 
                         i <- i+1
                         counter=counter+1
                         print(paste0("i: ", i))
                         print(paste0("counter: ", counter))
                    }
                    # IF POSITION IS NEGATIVE, THEN CHECK POSITION BEFORE END POSITION FOR START OF REGION
                    else{ 
                         
                         b <- end-1
                         while(x[b]>cutoff){ #go upstream of end position until position that doesn't meet cutoff is met
                              b <- b-1
                         }
                         start <- b+1

                         c <- end+1
                         while(x[c]>cutoff){ #go downstream of end position until position that doesn't meet cutoff is met
                              c <- c+1
                         }
                         end <- c-1
#                          print(paste0("start: ", start))
                         print(paste0("end: ", end))
                         if((end-start)>=window){ #If region is greater than window, then save start and end positions
                              storage <- rbind(storage, cbind(start, end, end-start+1, 1))
                              print(paste0("REGION FOUND1. START: ", start, "  END: ", end))
                              print(paste0("add: ", start, "; end: ", end))
                              break
                         }
                         i <- i+1
                    }
               
               }




               # IF ALL POSITION IN BETWEEN START/END ARE POSITIVE
               if(counter==(window-1)){
                    print(paste0("HIT: ", counter))

                    print(paste0(i, ":while loop done - ", x[i]))
                    # check to see if points after end position are positive
                    start <- start.save
                    while(x[end]>cutoff){ # THEN CHECK TO SEE HOW MANY POSITIONS AFTER END POSITION ARE POSITIVE
                         end <- end+1
                         print(paste0(i, ":end extension - ", x[end]))
                    }
                    end <- end-1
                    print(paste0("REGION FOUND2. START: ", start, "  END: ", end))
                    storage <- rbind(storage, cbind(start, end, end-start+1, 2))
                    print(paste0("add: ", start, "; end: ", end))
                    #                        storage <- rbind(storage, cbind(start, end, end-start))
                    print(paste0("NEW START POSITION: ", i))
                    #SAVE start position and check positions after end position for positives
               }
               #IF ALL POSITIONS IN BETWEEN START/END ARE NOT POSITIVE
               else{ 
                    print(paste0(i, ":while loop broken - ", x[i]))
                    
               }
          }
          else{ # IF END POSITION IS NEGATIVE, THEN NEW START END+1
               print(paste0(i, " - ", x[i], " :end position negative: ", x[end]))
               print("BREAK!")
#                break
          }
     }
     
     # PATTERN: NEGATIVE START, ? END
     else{
          start.save <- i
          end <- i+window          
          if(x[end]>cutoff){
               print(paste0(i, " - ", x[i], " :NEGATIVE start; ", end, " - ", x[end], " :POSITIVE end"))
               j <- end-1
               while(x[j]>cutoff){ # FIND START POSITION OF REGION
                    j<-j-1
               }
               start <- j+1
               
               k <- end+1
               while(x[k]>cutoff){ # FIND END POSITION OF REGION
                    k<-k+1
               }
               end <- k-1
               if((end-start)>=window-1){ #If region is greater than window, then save start and end positions
                    print(paste0("REGION FOUND3. START: ", start, "  END: ", end))
                    storage <- rbind(storage, cbind(start, end, end-start+1 , 3))
                    print(paste0("NEW START POSITION: ", i))
               }
               #IF END IS POSITIVE, THEN POTENTIAL START POSITION; CHECK EVERY POSITION BEFORE END POSITION FOR START OF REGION
          }else{ # IF END IS NEGATIVE, THEN NEW START POSITION AFTER END
               print(paste0(i, " - ", x[i], ":NEGATIVE start; ", end, " - ", x[end], " :NEGATIVE end"))
          }
     }
     i <- end+1 # NEW START POSITION
     print(paste0("THIS IS NEW START: ", i))
}
# proc.time()
rownames(storage) <- seq(1:nrow(storage))
storage
#      return storage
# }

> data.frame(x)
x
1    7.0
2   40.0
3   41.0
4   42.0
5   43.0
6   44.0
7   45.0
8   46.0
9   47.0
10  48.0
11  49.0
12  50.0
13  50.5
14   1.0
15   2.0
16   3.0
17  51.0
18  52.0
19   8.0
20   6.0
21  53.0
22  54.0
23  55.0
24  57.0
25  68.0
26  44.0
27  47.0
28  48.0
29  49.0
30  50.0
31   1.0
32   2.0
33   3.0
34  51.0
35  52.0
36   8.0
37   6.0
38  53.0
39  54.0
40  55.0
41  57.0
42  68.0
43  44.0
44  45.0
45  46.0
46  47.0
47  48.0
48  49.0
49  50.0
50   1.0
51   2.0
52   3.0
53  51.0
54  52.0
55   8.0
56   6.0
57  53.0
58  54.0
59  55.0
60  57.0
61  68.0
62  44.0
63  45.0
64  46.0
65  47.0
66  48.0
67  49.0
68  50.0
69   1.0
70   2.0
71   3.0
72  51.0
73  52.0
74   8.0
75   6.0
76  53.0
77  54.0
78  44.0
79  45.0
80  46.0
81  47.0
82  48.0
83  49.0
84  50.0
85   1.0
86   2.0
87   3.0
88  51.0
89  52.0
90   8.0
91   6.0
92  53.0
93  54.0
94  55.0
95  57.0
96  68.0
97  44.0
98  47.0
99  48.0
100 49.0
101 50.0
102  1.0
103  2.0
104  3.0
105 51.0
106 52.0
107  8.0
108  6.0
109 53.0
110 54.0
111 55.0
112 57.0
113 68.0
114 44.0
115 45.0
116 46.0
117 47.0
118 48.0
119 49.0
120 50.0
121  1.0
122  2.0
123  3.0
124 51.0
125 52.0
126  8.0
127  6.0
128 53.0
129 54.0
130 55.0
131 57.0
132 68.0
133 44.0
134 45.0
135 46.0
136 47.0
137 48.0
138 49.0
139 50.0
140  1.0
141  2.0
142  3.0
143 51.0
144 52.0
145  8.0
146  6.0
147 53.0
148 54.0
149 44.0
150 45.0
151 46.0
152 47.0
153 48.0
154 49.0
155 50.0
156  1.0
157  2.0
158  3.0
159 51.0
160 52.0
161  8.0
162  6.0
163 53.0
164 54.0
165 55.0
166 57.0
167 68.0
168 44.0
169 47.0
170 48.0
171 49.0
172 50.0
173  1.0
174  2.0
175  3.0
176 51.0
177 52.0
178  8.0
179  6.0
180 53.0
181 54.0
182 55.0
183 57.0
184 68.0
185 44.0
186 45.0
187 46.0
188 47.0
189 48.0
190 49.0
191 50.0
192  1.0
193  2.0
194  3.0
195 51.0
196 52.0
197  8.0
198  6.0
199 53.0
200 54.0
201 55.0
202 57.0
203 68.0
204 44.0
205 45.0
206 46.0
207 47.0
208 48.0
209 49.0
210 50.0
211  1.0
212  2.0
213  3.0
214 51.0
215 52.0
216  8.0
217  6.0
218 53.0
219 54.0

