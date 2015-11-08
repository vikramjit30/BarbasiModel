dir <- getwd()
setwd(dir)
library(igraph)

same_positive <<- 0
same_negative <<- 0

diff_positive <<- 0
diff_negative <<- 0

diff_comm <<- 0
same_comm <<- 0

## Reading all graph (.edges) files from the folder
files = list.files(pattern = "txt$") 

## Define the percentage of dir
percentage_diff_value = 28
difference_value = 0

file_name = paste("final_detection_","percentage_Diff",percentage_diff_value,".txt",sep="")

## function to detect communities
final_resuls <- array()
detect_Communities <- function(data,difference_value)
{
  
  #print(length(data))
  eachArea <- array() 
  eachColor <- array()
  eachDetect <- array()
  
  
  for (i in 1:length(data))
          {
            if (i==1 ) { node = data[i]
            }else if (i %% 2 == 0)  { eachArea <- c(eachArea,data[i])
                                      
            }else { eachColor <- c(eachColor,data[i])   } 
          }
          
          ## find the mean value 
          eachArea <- as.numeric(eachArea[-1])
          eachColor <- eachColor[-1]
          mean_value = mean(eachArea)
          
          difference_value = (percentage_diff_value * mean_value) / 100
         
          
          ## Classify them as same or different
          for (j in 1 : length(eachArea))
              {
                diff = abs(eachArea[j] - mean_value)
                if ( diff > difference_value ){
                      eachDetect <- c(eachDetect,"diff_comm")
                      
                }else{
                  eachDetect <- c(eachDetect,"same_comm")
                  }
              }
          
          eachDetect <- eachDetect[-1]
          
          ## Write the data into the final resutls file
          for (k in 1 : length(eachDetect))
          {
              result <- paste(eachColor[k]," ",eachDetect[k],sep="")
              # print(result)
              write(result, file = file_name,append=TRUE,sep="")
          }
          
        }



result <- paste("Actual_Communities "," ","Detected_Communitie ",sep="")
write(result, file = file_name,append=TRUE,sep="")

for ( areaFile in files)
{
  conn=file(areaFile,open="r")
  linn=readLines(conn)
  for (i in 1:length(linn))
    {
    
      ## Split each line based on "space"
      data = as.array(strsplit(linn[i], " "))[[1]]
      detect_Communities(data,difference_value)
      
    }
      close(conn)
}


## To find correctness of the results

## Read the files

conn=file(file_name,open="r")
linn=readLines(conn)
for (i in 1:length(linn))
{
  
  if (i == 1) next
  ## Split each line based on "space"
  ## <<- is for global values
  ## Read each line of the code and run the program   
  
  data = as.array(strsplit(linn[i], " "))[[1]]
  if ((data[1] == "black" && data[2] == "same_comm"))
      {
        same_positive <<- same_positive + 1
      }
  
  if ( (data[1]=="black" && data[2]=="diff_comm"))
      {
        same_negative  <<- same_negative + 1
      }
  
  
  if ((data[1] != "black" && data[2]=="diff_comm"))
      {
        diff_positive  <<- diff_positive + 1
      }
  
       
  if ((data[1] != "black" && data[2]=="same_comm"))
      {
        diff_negative <<- diff_negative + 1
      }
 
}

## Close the file 
close(conn)

write("\n ###Analysis###", file = file_name,append=TRUE,sep="")

total_cases = same_positive+same_negative+diff_positive+diff_negative

# Calculate the percentage value
correctnnes = (same_positive+diff_positive)*100/ total_cases

result = paste("same_positive: ",same_positive," same_negative: ",same_negative," diff_positive: ",diff_positive," diff_negative: ",diff_negative," total_cases: ",total_cases,sep="")
write(result, file = file_name,append=TRUE,sep="")


# Put the results on the same file and display at the end of the file [The final conclusion

result = paste("Correctness Value is: ",correctnnes,"%",sep="")
write(result, file = file_name,append=TRUE,sep="")



write("", file = file_name,append=TRUE,sep="")

###  End of the program  

