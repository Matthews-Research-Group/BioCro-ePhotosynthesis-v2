rm(list=ls())
library(PhotoGEA)
library(BioCro)
source('functions.R')

A_Ci_df = read.csv("test_data.csv")

df1 = A_Ci_df[,2:5]
colnames(df1) = c("PAR","Tleaf","Ci","A")

output = get_vcmax_jmax(df1)
# print(output)
vcmax1 = output[1]

df2 = A_Ci_df[,c(2:4,6)]
colnames(df2) = c("PAR","Tleaf","Ci","A")
output = get_vcmax_jmax(df2)
# print(output)
vcmax2 = output[1]
