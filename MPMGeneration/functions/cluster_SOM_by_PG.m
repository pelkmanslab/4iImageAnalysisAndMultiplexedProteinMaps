function [mcu_labels, path_save] = cluster_SOM_by_PG(path_som, name_output, path_output)
%%% Description
% This function extracts clusters the nodes generated by the
% Self-Organizing Map in R using FlowSome into Multiplexed Cell Units
% (MCUs) by Phenograph. 
% FlowSom:
% S. Van Gassen et al., FlowSOM: Using self-organizing maps for visualization and interpretation of cytometry data.
% Cytometry A 87, 636�645 (2015). doi: 10.1002/cyto.a.22625; Spmid: 25573116 
% Phenograph:
% J. H. Levine et al., Data-Driven Phenotypic Dissection of AML Reveals Progenitor-like Cells that Correlate with Prognosis. 
% Cell 162, 184�197 (2015). doi: 10.1016/ j.cell.2015.05.047; pmid: 26095251

% This code accompanies the publication
% G. Gut et al., Science 361, eaar7042 (2018). DOI: 10.1126/science.aar7042
% Author: Gabriele Gut, University of Zurich, Switzerland

%%% Inputs
% path_som = ''; % specify path to output of Script_03_ClusterMPPMatrix_with_SOM ending with '_ContributionNode.csv'
% name_output = name_output; % specify name, where all outputs should be saved, e.g. use same name as used in Script_01_GenerateMPPMatrix
% path_output = path_output; % specify directory, where all outputs should be saved, e.g. use same path as used in Script_01_GenerateMPPMatrix

%%% Outputs
% mcu_labels, a vector which contains the index to map SOM nodes to the newly identified MCUs
% path_save,  path to where mcu_labels was saved.

%% Load SOM node intensity profiles, mpp asignment to SOM nodes, node median intensity values
som_nodes = csvread(path_som,1,1);

%%%%% MCU identidication
%% Find best neighbourhood val
neigh_to_check = [2:1:15];
neighs = NaN(1,numel(neigh_to_check));
for ix = 1:numel(neigh_to_check)
    [tmp_fcu_labels,~,~] = phenograph(som_nodes, neigh_to_check(ix), 'jaccard'); % download latest version of Phenograph 
    neighs(ix) = numel(unique(tmp_fcu_labels));
end

plot(neigh_to_check,neighs)
set(gca, 'XTick',neigh_to_check, 'XTickLabels', num2strcell(neigh_to_check))
axis square
xlabel('Neighborhood value')
ylabel('Number of MCUs detected')

%% user input for neighborhood value
neigh_selected = input(sprintf('%s: Specify Neighborhood value where the number of MCUs starts to stabilise.\n Type X coordinate value in to Command Window and hit Enter.\n\n Type here: ',mfilename))
fprintf('You chose %d as a nieghborhood value. You will identify %d MCUs. \n',neigh_selected, neighs(neigh_selected))

%% Cluster SOM nodes in MCUs using Phenograph
neigh_selected = neigh_selected; % choose neighbouhood value where number of clusters stabilises
[mcu_labels] = phenograph(som_nodes, neigh_selected, 'jaccard'); % download latest version of Phenograph 

%% Save the pixel assignement to MCUs
% get ready for saving
date_time = datestr(datetime('now'));
date_time = strrep(date_time,' ','_');
date_time = strrep(date_time,':','');

name_save = sprintf('MCUAssignementSOMNodes_%s',name_output);
name_save = strcat(date_time, '_',name_save,'.mat');
path_save = fullfile(path_output,name_save);

save(path_save, 'mcu_labels', '-v7.3');
fprintf('MatLab file saved\n')
end

%% Subfunction
function c = num2strcell(n, format)
% num2strcell Convert vector of numbers to cell array of strings
% function c = num2strcell(n, format)
%
% If format is omitted, we use
% c{i} = sprintf('%d', n(i))

if nargin < 2, format = '%d'; end

N = length(n);
c = cell(1,N);
for i=1:N
  c{i} = sprintf(format, n(i));
end
end
  