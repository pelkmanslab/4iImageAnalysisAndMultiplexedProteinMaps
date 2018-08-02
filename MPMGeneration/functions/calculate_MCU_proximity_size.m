function [proximity_pooled, percarea_pooled, proximity_singlecell, percarea_singlecell, uq_mcu, path_output] = calculate_MCU_proximity_size(path_project, str_sharedfoldername, path_save_mculabels, path_somnodeid, path_mpp_processed, name_output, path_output)
%%% Description
% This function quantifies the spatial proximity score of all MCUs present
% in a cell for all cells included in the MPM analysis. It also quantifies
% the area of each MCU relative to the cell area in which the MCUs were
% present. 

% This code accompanies the publication
% G. Gut et al., Science 361, eaar7042 (2018). DOI: 10.1126/science.aar7042
% Author: Gabriele Gut, University of Zurich, Switzerland

%%% Inputs
% path_project = path_project; % specify path to project, see supporting pdf for expected directory structure, e.g. use same path as used in Script_01_GenerateMPPMatrix
% str_sharedfoldername = str_sharedfoldername; % specify string which is present in all folders of the individual cycles, e.g. use same name as used in Script_01_GenerateMPPMatrix
% path_save_mculabels = path_save_mculabels; % specify path to the .mat file generated in Script_04 (called path_save_mculabels) which contains list MCU indecs, e.g. path_save_mculabels
% path_id_som_nodes = ''; % specify path to output of Script_03_ClusterMPPMatrix_with_SOM ending with'_SOMNodeID.csv'
% path_save_processed = path_save_processed; % specify path to output of Script_02_ProcessMPPMatrix (path_save_processed), filename contains 'MPPMatrixProcessed' and ends with .mat
% name_output = name_output; % specify name, where all outputs should be saved, e.g. use same name as used in Script_01_GenerateMPPMatrix
% path_output = path_output; % specify directory, where all outputs should be saved, e.g. use same path as used in Script_01_GenerateMPPMatrix

%%% Outputs
% proximity_pooled, a 2D matrix which contains the mean spatial proximity score of each MCU with all others (pairwise cominations) over all analysed cells. Matrix is symmetric, MCU1 is is first row firs column, MCU last is last position in row and column
% percarea_pooled, a vector which contain the mean of the size of all MCUs normalized to the area of the cell they were measured in over all cells. Mean area of MCU1 is in the first row, Mean area of MCU last in the last row.
% proximity_singlecell, a 3D matrix which contains the pairwise spatial proximity score of each MCU with all others (pairwise cominations) over all analysed cells. Matrix is symmetric, MCU1 is is first row firs column, MCU last is last position in row and column. Cells are in 3rd dimension, where the first cell is in 1 and the last cell is in n, where n is the number of cells.
% percarea_singlecell,a 2D matrix which contains the size of all MCUs normalized to the area of the cell they were measured in. Rows represent MCUs, Columns represent cells, were MCU1 is in the first row and MCU last in the last, and were cell 1 is in the first column and the last cell is in column last. uq_mcu, labels of the MCUs, from 1 to n, where n is the last MCU.
% path_output, path to where all MCUproximity and size related variables are saved


%% Prepare paths
[name_cycle] = findfolderwithregexpi(path_project, str_sharedfoldername);
path_cycle = fullfile(path_project, name_cycle);
list_path = [1];
list_channel = [1];

%% Load SOM
id_nodes = csvread(path_somnodeid,1,1);
id_nodes = id_nodes(:,1);

%% Load MCU labels
mcu_labels = load(path_save_mculabels);
mcu_labels = mcu_labels.mcu_labels;

%% Load processed mpp files
mpp_input = load(path_mpp_processed);
cells_metadata = mpp_input.cells_metadata;
mpp_metadata = mpp_input.mpp_metadata;
clear mpp_input

%% Prepare matrices to save results
uq_mcu = unique(mcu_labels);
num_mcu = numel(uq_mcu);
percarea_singlecell = NaN(size(cells_metadata,1),num_mcu);
proximity_singlecell = NaN(num_mcu,num_mcu, size(cells_metadata,1));

%% loop over each cell to quantify MCU spatial proximity score profiles and size (% of cell area)
parfor cx = 1:size(cells_metadata,1)
    cell_oi = cx;
    
    %% get the cell and colorcode pixels
    [im_cell_stack, ~] = pick_segmentation_by_metadata2(cells_metadata(cx,:), 42, path_cycle(list_path), list_channel);
    
    %% show all clusters pooled
    im_pool = zeros(size(im_cell_stack{1},1),size(im_cell_stack{1},2),1);
    t = find(lin(im_cell_stack{1}));
    for ix = 1:numel(uq_mcu)
        
        label_to_check = uq_mcu(ix);
        
        % get overlaps
        node_oi = find(label_to_check == mcu_labels);
        linx_tmp_nodes = find(ismember(id_nodes(:,1),node_oi));
        linx_tmp_cell = find(mpp_metadata(:,1) == cell_oi);
        pix_phenocluster = intersect(linx_tmp_nodes,linx_tmp_cell);
        pix_phenocluster = mpp_metadata(pix_phenocluster,2);
        
        t2 = t(pix_phenocluster);
        im_pool(t2) = ix;
        
    end
    
    %% generate Adjency matrix
    tmp_adjacency = zeros(numel(uq_mcu),numel(uq_mcu));
    tmp_adjacency_corr = zeros(numel(uq_mcu),numel(uq_mcu));
    tmp_percentage = zeros(1,numel(uq_mcu));
    
    for ix = 1:numel(uq_mcu)
        bw_cluster = im_pool == uq_mcu(ix);
        bw_cluster_source = bwmorph(bw_cluster, 'clean',1);
        bw_cluster = imdilate(bw_cluster_source,strel('square',3));
        bw_cluster(bw_cluster_source) = 0;
        
        % percentage of interaction
        neigh_clust = unique(im_pool(bw_cluster));
        countElA=histc(im_pool(bw_cluster),neigh_clust);
        countElA(neigh_clust == 0) = [];
        neigh_clust(neigh_clust == 0) = [];
        countElA(neigh_clust == uq_mcu(ix)) = [];
        neigh_clust(neigh_clust == uq_mcu(ix)) = [];
        neigh_clust(neigh_clust == 0) = [];
        relFreq=countElA/sum(countElA);
        
        % orig without randomization
        tmp_adjacency(ix,neigh_clust) = relFreq;
        tmp_adjacency(neigh_clust,ix) = relFreq;
        
        %% randomized control for djacency
        rand_iteration = 100;
        tmp_adjacency_random = NaN(rand_iteration,numel(uq_mcu));
        for rx = 1:rand_iteration
            im_pool_rand = zeros(size(im_pool));
            rand_position = t(randsample(numel(t),numel(t)));
            im_pool_rand(t)  = im_pool(rand_position);
            
            rand_bw_cluster = t(randsample(sum(bw_cluster(:)),sum(bw_cluster(:))));%new
            neigh_clust_random = unique(im_pool_rand(rand_bw_cluster));% new
            countElA_random = histc(im_pool_rand(rand_bw_cluster),neigh_clust_random);% new
            
            countElA_random(neigh_clust_random == 0) = [];
            neigh_clust_random(neigh_clust_random == 0) = [];
            countElA_random(neigh_clust_random == uq_mcu(ix)) = [];
            neigh_clust_random(neigh_clust_random == uq_mcu(ix)) = [];
            relFreq_random=countElA_random/sum(countElA_random);
            
            tmp_adjacency_random(rx,neigh_clust_random) = relFreq_random;
            
        end
        tmp_adjacency_random = nanmean(tmp_adjacency_random(:,neigh_clust),1)';
        
        %% correct score for randomness
        if ~isempty(neigh_clust)
            tmp_adjacency_corr(ix,neigh_clust) = relFreq-tmp_adjacency_random;
            tmp_adjacency_corr(neigh_clust,ix) = relFreq-tmp_adjacency_random;
        end
        
        %% calculate percentage of are used
        tmp_percentage(1,ix) = sum(bw_cluster_source(:));
    end
    
    proximity_singlecell(:,:,cx) = tmp_adjacency_corr;
    percarea_singlecell(cx,:) = tmp_percentage./sum(tmp_percentage,2);
    
end

%% combine analysis here
%get means
percarea_pooled = nanmean(percarea_singlecell,1)';
percarea_pooled = percarea_pooled./sum(percarea_pooled);
proximity_pooled = nanmean(proximity_singlecell,3);
proximity_pooled(find(diag(ones(size(proximity_pooled,1),1)))) = 0;

%% save in super matrix
path_output = fullfile(path_output, sprintf('MCU_Proximity_Size_%s.mat', name_output));
save(path_output,'mcu_labels', 'cells_metadata', 'proximity_pooled', 'percarea_pooled', 'percarea_singlecell','proximity_singlecell', 'uq_mcu', '-v7.3')

fprintf('Saved MatLab file\n')
fprintf('Finished\n')
end

%% Subfunctions
function cellstrFileName = findfolderwithregexpi(strRootPath,strRegExp, boolFullPath)

if nargin<3
    boolFullPath = false;
end

% get directory content listing
cellTargetFolderList = CPdir(strRootPath)';
cellTargetFolderList(~[cellTargetFolderList.isdir]) = [];
cellTargetFolderList = {cellTargetFolderList.name}';
cellTargetFolderList([1 2]) = [];

% get rid of non-"Well" directories
matMatch = ~cellfun(@isempty,regexpi(cellTargetFolderList,strRegExp));

% if one hit is found, return string, otherwise cell-array
if size(find(matMatch))==1
    cellstrFileName = cellTargetFolderList{matMatch};
else
    cellstrFileName = cellTargetFolderList(matMatch);
end
% convert to full path if requested
if boolFullPath && ~isempty(cellstrFileName)
    if iscell(cellstrFileName)
        cellstrFileName = cellfun(@(x) fullfile(strRootPath,x),cellstrFileName,'UniformOutput',false);
    else
        cellstrFileName = fullfile(strRootPath,cellstrFileName);
    end
end
end

function cellstrFileName = findfilewithregexpi(strRootPath,strRegExp)

% get directory content listing
cellTargetFolderList = CPdir(strRootPath)';
cellTargetFolderList([cellTargetFolderList.isdir]) = [];
cellTargetFolderList = {cellTargetFolderList.name}';

% get rid of non-"Well" directories
matMatch = ~cellfun(@isempty,regexpi(cellTargetFolderList,strRegExp));

% if one hit is found, return string, otherwise cell-array
if size(find(matMatch))==1
    cellstrFileName = cellTargetFolderList{matMatch};
else
    cellstrFileName = cellTargetFolderList(matMatch);
end

end

function [im_cell_stack, im_nuclei_stack, im_seg_cells] = pick_segmentation_by_metadata2(cell_metadata, site_per_well, list_path, list_channel)

%% identify position of de cell in the dataset
row_names = {'A' 'B' 'C' 'D' 'E' 'F' 'G' 'H' 'I' 'J' 'K' 'L' 'M' 'N'};

position_row = row_names{cell_metadata(1)};
position_col = sprintf('%02d', cell_metadata(2));

site = rem(cell_metadata(3),site_per_well);
if site == 0
    position_site = sprintf('%03d',site_per_well);
else
    position_site = sprintf('%03d',site);
end

position_cellid = cell_metadata(4);

position_search = sprintf('.*%s%s.*F%s', position_row,position_col,position_site)


%% get both image and sementation
im_cell_stack = cell(numel(list_channel),1);
im_nuclei_stack = cell(numel(list_channel),1);

for ix = 1:numel(list_channel)
    fprintf('Loading %d of %d\n',ix, numel(list_channel))
    
    path_tmp = list_path{ix};
    
    %% new
    all_files = findfilewithregexpi(fullfile(path_tmp,'SEGMENTATION'),'.png');
    name_seg_cells = char(all_files(~cellfun(@isempty,regexpi(all_files,strcat(position_search,'.*','Cells')))));
    name_seg_nuclei = char(all_files(~cellfun(@isempty,regexpi(all_files,strcat(position_search,'.*','Nuclei')))));
    
    im_seg_cells = imread(fullfile(fullfile(path_tmp,'SEGMENTATION'), name_seg_cells));
    im_seg_nuclei = imread(fullfile(fullfile(path_tmp,'SEGMENTATION'), name_seg_nuclei));
    
    
    %% grab single cell
    BB = regionprops(im_seg_cells == position_cellid,'BoundingBox');
    boundingbox_details = round(BB.BoundingBox);
    
    %% store
    im_cell_stack{ix} = imcrop(im_seg_cells == position_cellid,boundingbox_details);
    im_nuclei_stack{ix} = imcrop(im_seg_nuclei == position_cellid,boundingbox_details);
    
end
end

