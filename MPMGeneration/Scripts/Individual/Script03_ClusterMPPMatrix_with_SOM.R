## Specify inputs

# Specify load path of processed MPP Matrix (Output of Script02)
path_load_data = ''; # path to the file containing "_MPPMatrixProcessed" and finishing with ".csv", copy pase from Matlab Command Window, after executing Script_02_ProcessMPPMatrix.m

# Specify name of your experiment, for purposes of saving results files
name_output = '' # specify a name, with which you want to label all outputs, e.g. use same name as used in Script_01, e.g. Test

# Specify some variables for the SOM clustering
ydim_SOM = 20; # you can use anything from 10 to 80, e.g. 20
xdim_SOM = 20; # you can use anything from 10 to 80
dist_metric = 2; # distance matric, e.g. 2
num_runs = 20; # number of runs, e.g. 20

##

# Initiate packages
library(FlowSOM)

# Prepare output path

path_output = dirname(dirname(path_load_data));

# Load processed MPP Matrix
data = read.csv(path_load_data,header=FALSE);

# Convert to numeric matrix
data = data.matrix(data)

# Select protein marker columns to use for clustering
marker_cols <- 1:ncol(data)

##### RUN FLOWSOM
# create flowFrame object
data_FlowSOM <- flowCore::flowFrame(data)

# set seed for reproducibility
set.seed(1234)

# run FlowSOM (initial steps prior to meta-clustering)
out <- FlowSOM::ReadInput(data_FlowSOM, transform = FALSE, scale = FALSE)
out <- FlowSOM::BuildSOM(out, colsToUse = marker_cols, ydim = ydim_SOM, xdim = xdim_SOM, distf = dist_metric,rlen = num_runs)


##### Save ######
# make the date string
date = Sys.Date();
date = gsub("//-", "", date)

# make the name for complete SOM output and parts of it
str_save_name1 = paste0(date, "_",name_output, "_", "CompleteSOMOutput.rDS")
str_save_name2 = paste0(date, "_",name_output, "_", "SOMNodeID.csv")
str_save_name3 = paste0(date, "_",name_output, "_", "FeatureMediansPerNode.csv")
str_save_name4 = paste0(date, "_",name_output, "_", "ContributionNode.csv")
str_save_name5 = paste0(date, "_",name_output, "_", "CoordinatesNode.csv")
str_save_name6 = paste0(date, "_",name_output, "_", "SizesNode.csv")

# generate save directory
str_save_dir = paste0(path_output,date, "_ROutput")
dir.create(str_save_dir)

# Save whole result structure and parts of it in newly created directory
str_save_path = file.path(str_save_dir,str_save_name1)
saveRDS(out,str_save_path)

# get for each gene the cluster ID
str_save_path = file.path(str_save_dir,str_save_name2)
GeneClusterID = out$map$mapping
write.csv(GeneClusterID, str_save_path)

# get medians of all features per node
str_save_path = file.path(str_save_dir,str_save_name3)
FeatureMediansPerNode = out$map$medianValues
FeatureMediansPerNode[is.na(FeatureMediansPerNode)] = 0; # very rarely a node might contain only NAs They are then replaced by zeros.
write.csv(FeatureMediansPerNode, str_save_path)

# get feature importance per node
str_save_path = file.path(str_save_dir,str_save_name4)
contributionNode = out$map$codes
write.csv(contributionNode, str_save_path)

# get x and y coordinates of nodes
str_save_path = file.path(str_save_dir,str_save_name5)
coordinatesNode = out$MST$l
write.csv(coordinatesNode, str_save_path)

# get sizes of nodes
str_save_path = file.path(str_save_dir,str_save_name6)
sizesNode = out$MST$size
write.csv(sizesNode, str_save_path)