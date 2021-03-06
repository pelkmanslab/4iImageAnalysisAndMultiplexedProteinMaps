%% Script_05_Proximity_Size_Calculation                                           
% Specify inputs                                            
path_project = ''; % specify path to project, see supporting pdf for expected directory structure, e.g. use same path as used in Script_01
str_sharedfoldername = ''; % specify string which is present in all folders of the individual cycles, e.g. use same name as used in Script_01
path_save_mculabels = ''; % specify path to the .mat file generated in Script_04 (called path_save_mculabels) which contains list MCU indecs, e.g. path_save_mculabels
path_id_som_nodes = ''; % specify path to output of Script_03_ClusterMPPMatrix_with_SOM ending with'_SOMNodeID.csv'
path_mpp_processed = ''; % specify path to output of Script_02_ProcessMPPMatrix (path_save_processed), filename contains 'MPPMatrixProcessed' and ends with .mat
name_output = ''; % specify name, where all outputs should be saved, e.g. use same name as used in Script_01
path_output = ''; % specify directory, where all outputs should be saved, e.g. use same path as used in Script_01

% Run calculate_FCU_proximity_size function. It will calculate MCU Spatial Proximity Scores and MCU sizes
[proximity_pooled, percarea_pooled, proximity_singlecell, percarea_singlecell, uq_mcu, path_save_proxsize] = calculate_MCU_proximity_size(  path_project,...
                                                                                                                        str_sharedfoldername,...
                                                                                                                        path_save_mculabels,...
                                                                                                                        path_id_som_nodes,...
                                                                                                                        path_mpp_processed,...
                                                                                                                        name_output,...
                                                                                                                        path_output); 