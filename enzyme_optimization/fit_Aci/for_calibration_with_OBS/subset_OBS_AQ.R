energy_farm_path="/Users/yufeng/Desktop/UIUC/Research/Github/energy-farm-biocro/"

bb_data = read.csv(paste0(energy_farm_path,'ball_berry_curves_2022/bb_soybean_ld11_2022.csv'))
A_Q = bb_data[,c("A","Qin","Ci","Tleaf","curve_identifier")]
A_Q = A_Q[-c(1,2),]

require(ggplot2)
library(data.table)
A_Q = data.frame(A_Q)
A_Q$A = as.numeric(A_Q$A)
A_Q$Qin = as.numeric(A_Q$Qin)
A_Q$Ci = as.numeric(A_Q$Ci)
A_Q$Tleaf = as.numeric(A_Q$Tleaf)

A_Q_217 = A_Q#[A_Q$curve_identifier %like% "217",]
# ggplot(A_Q_217, aes(x=Qin, y=A,colour=curve_identifier)) + geom_line()

write.csv(A_Q_217,"ACi_AQ_data/AQ_DOY217_2022.csv",row.names = FALSE)

# Ci = as.numeric(A_Q_217$Ci)
# Ci = matrix(Ci, nrow=7)
# Ci_mean = rowMeans(Ci)
# 
# Q = as.numeric(A_Q_217$Qin)
# Q = matrix(Q,nrow=7)
# Q_mean = rowMeans(Q)
# 
# write.csv(cbind(Ci_mean,Q_mean),"AQ_data/AQ_DOY217_2021_mean.csv",row.names = FALSE)
