energy_farm_path="/Users/yufeng/Desktop/UIUC/Research/Github/energy-farm-biocro/"

aci_data = read.csv(paste0(energy_farm_path,'aci_curves_2022/aci_soybean_ld11_2022.csv'))
A_Ci = aci_data[,c("A","Ci","Tleaf","Qin","curve_identifier")]
A_Ci = A_Ci[-c(1,2),]

require(ggplot2)
library(data.table)
A_Ci = data.frame(A_Ci)
A_Ci$A = as.numeric(A_Ci$A)
A_Ci$Ci = as.numeric(A_Ci$Ci)
A_Ci$Tleaf = as.numeric(A_Ci$Tleaf)
A_Ci$Qin = as.numeric(A_Ci$Qin)

ggplot(A_Ci, aes(x=Ci, y=A,colour=curve_identifier)) +
geom_line()

write.csv(A_Ci,"ACi_AQ_data/ACi_DOY217_2022.csv",row.names = FALSE)
