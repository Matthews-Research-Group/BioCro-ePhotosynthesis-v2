library(stringr)
# Function to read the second-to-last line from a text file
read_last_two_lines <- function(file_path) {
  # Read all lines from the file into a character vector
  lines <- readLines(file_path)
  
  # Check if the file has at least two lines
  if (length(lines) < 2) {
    stop("The file does not contain enough lines.")
  }
  
  # Return the second-to-last line
  second_to_last_line <- lines[length(lines) - 1]
  # Return the last line
  last_line <- lines[length(lines)]
  return(c(second_to_last_line,last_line))
}
# Function to extract the number after "Penalty is" from a line
extract_penalty <- function(line) {
  # Use regular expression to find the number after "Penalty is"
  match <- regexpr("Penalty is\\s*([0-9.]+)", line)
  if (match != -1) {
    # Extract the matched substring and convert to numeric
    penalty <- regmatches(line, match)
    # Extract the number part and convert to numeric
    penalty_value <- as.numeric(sub("Penalty is\\s*", "", penalty))
    return(penalty_value)
  } else {
    #stop("No match found for 'Penalty is'.")
    return(NaN) 
  }
}

#ephoto_single.exe is linked to the Conda's libstdc++.so.6.
#somehow in R session, the OS's libstdc++.so.6 is prioritized. The two differ in version
#you can check that by running "ldd ephoto_single.exe" to see its dependencies
#here we use LD_PRELOAD specifically for R sessions to prioritize libraries from Conda
Sys.setenv(LD_PRELOAD = "/home/n-z/yufenghe/.conda/envs/ephoto/lib/libstdc++.so.6")

reps = 1:5 #name the case
rerun= FALSE #is it a rerun?
for (rep in reps){
  #get the inputs of parameters of scaling factors
  if(rerun){
    para_best_all = readRDS(paste0("para_best_reRun_Rep",rep,".rds"))
  }else{
    para_best_all = readRDS(paste0("para_best_all_Rep",rep,".rds"))
  }
  output_matrix = c()
  Actl = c()
  for (i in 1:nrow(para_best_all)){
      print(paste("processing rep",rep,"case",i))
      para_best = para_best_all[i,]
  
      r = para_best[1]
      t = para_best[2]
      c = para_best[3]

#run optimized case
      opt_scaling_factors = para_best[4:28]
      write(opt_scaling_factors, 
                  file = "optimized_enzyme_scaling_factors.txt",
            ncolumns = length(opt_scaling_factors), sep = " ")
      cmd = paste0("./ephoto_single.exe -r ",r," -t ",t," -c ",c," >& log")
      print(cmd)
      system(cmd)
      Sys.sleep(2)
      two_lines = read_last_two_lines("log")
      penalty_line = two_lines[1] 
      penalty      = extract_penalty(penalty_line)
      Avalue       = as.numeric(two_lines[2])
      last_metabolite = read.table("last_data.txt")
      if(is.nan(penalty)) last_metabolite[] = NaN
      output = c(r,t,c,penalty,last_metabolite)
      output_matrix = rbind(output_matrix,output)
#run control
#control needs to be run only once without rerun
      if(!rerun){
       opt_scaling_factors[] = 1 
       write(opt_scaling_factors, 
                   file = "optimized_enzyme_scaling_factors.txt",
             ncolumns = length(opt_scaling_factors), sep = " ")
       cmd = paste0("./ephoto_single.exe -r ",r," -t ",t," -c ",c," >& log")
       system(cmd)
       Sys.sleep(2)
       two_lines = read_last_two_lines("log")
       Actl       = c(Actl,as.numeric(two_lines[2]))
      }
  }
  colnames(output_matrix)[1] = "Q"
  colnames(output_matrix)[2] = "T"
  colnames(output_matrix)[3] = "Ci"
  colnames(output_matrix)[4] = "Penalty"
  if(rerun){
    saveRDS(output_matrix,paste0("steady_state_metabolite_rep",rep,"_rerun.rds"))
  }else{
    saveRDS(output_matrix,paste0("steady_state_metabolite_rep",rep,".rds"))
    para_best_all_with_control = cbind(para_best_all,Actl)
    saveRDS(para_best_all_with_control,paste0("para_best_all_Rep",rep,"_withCTL.rds"))
  }
}
