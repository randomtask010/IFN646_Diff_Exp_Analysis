
# Running through the edgeR example

source("Working Directory/EdgeR_Analysis/edgeR_parameterised_file.R")

#Variables

Tool <- "edgeR"
SourceFileVariable <- c("3_500_500", "3_750_250", "3_1000_0", "6_500_500", "6_750_250", "6_1000_0", "9_500_500", "9_750_250", "9_1000_0")
PValue <- 0.05

# loop

for (sample in SourceFileVariable){
  run_loop_edgeR(Tool, sample, PValue)
}
