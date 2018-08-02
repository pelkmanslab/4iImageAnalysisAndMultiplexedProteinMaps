%% This script is a compilation of all 6 scripts which build an Multiplexed Protein Map. Execute the individual Scripts out of this document. 
%% (except for Script_03 for which you will need to swith to R), or execute each step from a separate .m file.
%% Do not forget to specify input arguments for each script/section. Example values are given at the end of expenations of each variable afer an e.g.
%% If you execute the Scripts in sequence from this file, we have already filled in variable names in later sections, to make your life easier.
%% Enjoy

%% Script_01_GenerateMPPMatrix
% Specify inputs
path_project = ''; % specify path to project, see supporting pdf for expected directory structure
str_sharedfoldername = ''; % specify string which is present in all folders of the individual cycles
path_cellmetadata_info = ''; % specify path to metadata of cells, from which a Multiplexed Protein Map should be built
sites_per_well = []; % number of microscopy sites acquired per well, e.g. 42
str_seg_object = ''; % specify the segmetation object, which you want to use as a mask to extract multiplexed pixel profiles, e.g. Cells
channel_dapi = []; % Specify digit which specifies DAPI signal of your 4i cycle in your filenames e.g. 1 for "Images_C03_T0001F001A01Z01C01.png"
channel_green = []; % Specify digit which specifies second signal of your 4i cycle (often green) in your filenames e.g. 2 for "Images_C03_T0001F001A01Z01C02.png"
channel_red = []; % Specify digit which specifies third signal of your 4i cycle (often green) in your filenames e.g. 3 for "Images_C03_T0001F001A01Z01C03.png"
num_parallel = []; % specify number of parallelizations you want, e.g. 8
name_output = ''; % specify a name, with which you want to label all outputs, e.g. Test
path_output = ''; % specify directory, where all outputs should be saved

% Run extract_MPPMatrix function. It will generate MMP Matrix of num_cells specified randomly selected cells and save MPP Matrix specified directory
% extract_MPPMatrix will extract DAPI signal (channel_dapi) only for the first cycle.
path_output = fullfile(path_output, 'MatlabOutput');
[mpp_matrix, mpp_metadata, cells_metadata, list_channel, path_save_raw] = extract_MPPMatrix(   path_project,...    
                                                                                                                str_sharedfoldername,...
                                                                                                                path_cellmetadata_info,...
                                                                                                                sites_per_well,...
                                                                                                                str_seg_object,...
                                                                                                                channel_dapi,...
                                                                                                                channel_green,...
                                                                                                                channel_red,...
                                                                                                                num_parallel,...
                                                                                                                name_output,...
                                                                                                                path_output);                                                                                                      
clear cells_MPPMatrix sites_per_well str_seg_object channel_dapi channel_green channel_red num_parallel

%% Script_02_ProcessMPPMatrix
% Specify inputs
path_MPPMatrix = path_save_raw; % specify path to the .mat file generated in Script_01 (called path_save_raw) which contains MPP related variables, e.g. path_save_raw
stains_to_exclude = []; % specify which 4i channels should be excluded from further analysis
val_background = []; % specify background intensity value of the camera chip  e.g. 110
val_quantile = []; % specify quantile value, which should be used to rescale intensties of each 4i channel between 0 and 100, where val_quantile will be 100, e.g. 0.98
val_threshold = []; % specify which singal value should be used to exclude MPP profile, if all 4i channel intensities lay below val_threshold, e.g. 33.333333333
name_output = name_output; % specify a name, with which you want to label all outputs, e.g. use same name as used in Script_01
path_output = path_output; % specify directory, where all outputs should be saved, e.g. use same path as used in Script_01

% Run process_MPPMatrix_MPPmetadata function. It will generate a cleaned up version of MMP Matrix and MPP metadata according to the specifications provided.
[mpp_matrix, mpp_metadata, cells_metadata, path_save_processed] = process_MPPMatrix_MPPmetadata(    path_MPPMatrix,...
                                                                                                    stains_to_exclude,...
                                                                                                    val_background,...
                                                                                                    val_quantile,...
                                                                                                    val_threshold,...
                                                                                                    name_output,...
                                                                                                    path_output);
clear stains_to_exclude val_background val_quantile val_threshold

%% Script_03_ClusterMPPMatrix_with_SOM
% This step is performed in R. Please download the most recent version of R and FlowSOM.
% Ourputs of the Script will be saved in a newly created folder in the directory path_project/ROutput.
% The outputs used to generate an MPM are the following:
% '_ContributionNode.csv'contains the Self-Oranizing Map (SOM)
% '_SOMNodeID.csv' contains an index (first column) which assignes each MPP to a SOM node
% For additional description on the outputs of FlowSOM, please consult
% S. Van Gassen et al., FlowSOM: Using self-organizing maps for visualization and interpretation of cytometry data.
% Cytometry A 87, 636–645 (2015). doi: 10.1002/cyto.a.22625; Spmid: 25573116




%% Script_04_ClusterSOMNodes_with_PG
% Specify inputs
path_som = ''; % specify path to output of Script_03_ClusterMPPMatrix_with_SOM ending with '_ContributionNode.csv'
name_output = name_output; % specify name, where all outputs should be saved, e.g. use same name as used in Script_01
path_output = path_output; % specify directory, where all outputs should be saved, e.g. use same path as used in Script_01

% Run cluster_SOM_by_PG function, select neighborhood value interactively and identify MCUs.
[mcu_labels, path_save_mculabels] = cluster_SOM_by_PG(  path_som,...
                                                        name_output,...
                                                        path_output);
                                                   
                                                    
%% Script_05_Proximity_Size_Calculation                                         
% Specify inputs                                            
path_project = path_project; % specify path to project, see supporting pdf for expected directory structure, e.g. use same path as used in Script_01
str_sharedfoldername = str_sharedfoldername; % specify string which is present in all folders of the individual cycles, e.g. use same name as used in Script_01
path_save_mculabels = path_save_mculabels; % specify path to the .mat file generated in Script_04 (called path_save_mculabels) which contains list MCU indecs, e.g. path_save_mculabels
path_id_som_nodes = ''; % specify path to output of Script_03_ClusterMPPMatrix_with_SOM ending with'_SOMNodeID.csv'
path_save_processed = path_save_processed; % specify path to output of Script_02_ProcessMPPMatrix (path_save_processed), filename contains 'MPPMatrixProcessed' and ends with .mat
name_output = name_output; % specify name, where all outputs should be saved, e.g. use same name as used in Script_01
path_output = path_output; % specify directory, where all outputs should be saved, e.g. use same path as used in Script_01

% Run calculate_FCU_proximity_size function. It will calculate MCU Spatial Proximity Scores and MCU sizes
[proximity_pooled, percarea_pooled, proximity_singlecell, percarea_singlecell, uq_mcu, path_save_proxsize] = calculate_MCU_proximity_size(  path_project,...
                                                                                                                                            str_sharedfoldername,...
                                                                                                                                            path_save_mculabels,...
                                                                                                                                            path_id_som_nodes,...
                                                                                                                                            path_save_processed,...
                                                                                                                                            name_output,...
                                                                                                                                            path_output); 
                                                                                                                    
%% Script_06_MPMGeneration                                         
% Specify inputs                                            
path_som = path_som; % specify path to output of Script_03_ClusterMPPMatrix_with_SOM ending with '_ContributionNode.csv'
path_save_proxsize = path_save_proxsize; % specify path to output of Script_05_Proximity_Size_Calculation (path_save_proxsize) containing the string "MCU_Proximity_Size"
perplexity = []; % Often 3, try what looks best for you for the tSNE projection of the Multiplexed Protein Map. e.g. 3
std_display_limit_edges = []; % Threshold for displaying the SPS edges between nodes in MPM e.g. 2.5
scaling_factor_edges = []; % Scales the thicknes of the edges in the MPM, e.g. 20
scaling_factor_nodes = []; % Scales the diameter of the nodes in the MPM, e.g. 10000
list_stains = arrayfun(@(x) sprintf('%02d',x),[1:(size(mpp_matrix,2))], 'UniformOutput',false); % this will generate as list of strings starting from Stain01, finishing with StainN, where N is the number of 4i channels measured.
% list_stains = {'Type in 4iChannel 01' 'Type in 4iChannel 02' 'Type in 4iChannel N' 'etc' ''Type in 4iChannel last''}; Alternatively, type in 4i  channels manually.

% run plot_MPM_and_MCULoadings function, to siplay MPM and MCU intensity loadings
[map_prox_coordinates, mcu_mean_int_zs] = plot_MPM_and_MCULoadings( path_som,...
                                                                    path_save_proxsize,...
                                                                    perplexity,...
                                                                    std_display_limit_edges,...
                                                                    scaling_factor_edges,...
                                                                    scaling_factor_nodes,...
                                                                    list_stains);

