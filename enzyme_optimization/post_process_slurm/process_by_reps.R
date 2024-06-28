# Load the stringr package
library(stringr)
#slurm_patterns = c("slurm-871062","slurm-871253","slurm-871280","slurm-871307")
slurm_patterns = c("slurm-873279")
results_path = "../slurm_results_test"
repeats = 1 #the nth repetition
for (ptn in slurm_patterns){
  all_files = dir(results_path,pattern=ptn,full.names=TRUE)
  all_files_ordered = all_files[order(as.numeric(gsub("\\D", "", all_files)))]
  print((all_files_ordered))
  #first, we have to remove some lines with the CVODE warning messages
  #i use a vim script to do this
  if(TRUE){
    env_conditions = c() 
    for (file in all_files_ordered){
     first_line <- readLines(file, n = 1)
    # Extract numbers from the string
     numbers <- str_extract_all(first_line, "\\d+")[[1]]
     numbers <- as.numeric(numbers)
     env_conditions = rbind(env_conditions,numbers)
     cmd = paste("vi -c 'source scripts.vim'",file) 
     # Call the system command
     system(cmd)
     Sys.sleep(2)
    }
    colnames(env_conditions) = c("Q","T","Ci")
    saveRDS(env_conditions,"env_conditions.rds")
  }else{
    env_conditions = readRDS("env_conditions.rds")
  }
  print("starting second step......")
  #second, we get the last para values 
  para_best_all = c()
  for (i in 1:length(all_files_ordered)){
    file = all_files_ordered[i]
    print(file)
    x = read.table(file)
    enzyme_info = read.table('../ProteinContentCal_simple.txt')
    enzyme_names = enzyme_info$V1
    enzyme_names = as.character(enzyme_names)
  
    para_best = x[nrow(x),]
    para_best = as.numeric(para_best[5:30])
    names(para_best) = c(enzyme_names,"Aopt")
    para_best_all = rbind(para_best_all,para_best)
  
    totalE = sum(enzyme_info$V3*para_best[1:25]/enzyme_info$V5*enzyme_info$V6/1000)
    print(totalE)
  }
  para_best_all = cbind(env_conditions,para_best_all)
  saveRDS(para_best_all,paste0("para_best_all_",repeats,".rds"))
  repeats = repeats + 1
}
