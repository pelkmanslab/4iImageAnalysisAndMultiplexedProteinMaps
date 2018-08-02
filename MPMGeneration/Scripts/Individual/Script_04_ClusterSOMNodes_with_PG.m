%% Script_04_ClusterSOMNodes_with_PG
% Specify inputs
path_som = ''; % specify path to output of Script_03_ClusterMPPMatrix_with_SOM ending with '_ContributionNode.csv'
name_output = ''; % specify name, where all outputs should be saved, e.g. use same name as used in Script_01
path_output = ''; % specify directory, where all outputs should be saved, e.g. use same path as used in Script_01

% Run cluster_SOM_by_PG function, select neighborhood value interactively and identify MCUs.
[mcu_labels, path_save_mculabels] = cluster_SOM_by_PG(  path_som,...
                                                        name_output,...
                                                        path_output);