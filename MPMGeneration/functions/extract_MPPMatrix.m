function [mpp_matrix, mpp_metadata, cells_metadata, list_channel, path_save] = extract_MPPMatrix(path_project, str_sharedfoldername, path_cellmetadata_info, sites_per_well, str_seg_object, channel_dapi, channel_green, channel_red, num_parallel, name_output, path_output)
%%% Description
% This function extracts Multiplexed Pixel Profiles from 4i images for all
% pixel positions within the segmentation of a specified cell
% (cells_metadata). The function and its subfunctions expect a certain
% directory structure as well as some nameing conventions for images and
% segmentation files. Please read the pdf file named
% Explenations_Generation_MulriplexedProteinMaps.pdf for more information.
% 
% This code accompanies the publication
% G. Gut et al., Science 361, eaar7042 (2018). DOI: 10.1126/science.aar7042
% Author: Gabriele Gut, University of Zurich, Switzerland

%%% Inputs
% path_project = ''; % specify path to project, see supporting pdf for expected directory structure
% str_sharedfoldername = ''; % specify string which is present in all folders of the individual cycles
% path_cellmetadata_info = ''; % specify path to metadata of cells, from which a Multiplexed Protein Map should be built
% sites_per_well = []; % number of microscopy sites acquired per well, e.g. 42
% str_seg_object = ''; % specify the segmetation object, which you want to use as a mask to extract multiplexed pixel profiles, e.g. Cells
% channel_dapi = []; % Specify digit which specifies DAPI signal of your 4i cycle in your filenames e.g. 1 for "Images_C03_T0001F001A01Z01C01.png"
% channel_green = []; % Specify digit which specifies second signal of your 4i cycle (often green) in your filenames e.g. 2 for "Images_C03_T0001F001A01Z01C02.png"
% channel_red = []; % Specify digit which specifies third signal of your 4i cycle (often green) in your filenames e.g. 3 for "Images_C03_T0001F001A01Z01C03.png"
% num_parallel = []; % specify number of parallelizations you want, e.g. 8
% name_output = ''; % specify a name, with which you want to label all outputs, e.g. Test
% path_output = ''; % specify directory, where all outputs should be saved

%%% Outputs
% mpp_matrix, a 2D matrix containing the Multiplexed Pixel Profiles (MPP) extracted from the cells specified in cells_metadata. Rows correspond to MPPs, columns the 4i channels in which the pixel intensities were measured
% mpp_metadata, a 2D matrix containing the metadata information of each MPP. Rows correspond to MPPs MPPMatrix, the first column contains the cell ID, the second column the pixel ID.
% cells_metadata, a 2D matrix containing the metadata information of cells in CellProfiler 1 format. Rows correspond to cells, the first columns contains the row identifier (in numneric, e.g. A = 1, B = 2) of the well of the cell. The second columns contains the column identifier of the well of the cell. The Third columns contains the microscopy site identifier in which the cell was found (e.g. if the cell was imaged in the 5 microscopy site of the experiment, the identifier should say 5).
% list_channel, a vector which contains the channel which were extracted for each cycle.
% path_save, path to where all MPP related files have been saved.

%% Load cell meta data
cells_metadata = csvread(path_cellmetadata_info)

%% Generate paths of each cycle
[name_cycle] = findfolderwithregexpi(path_project, str_sharedfoldername); % identifies all sub folders in the project directory, which include str_sharedfoldername in their name.
path_cycle = fullfile(path_project, name_cycle); % generates paths to all subdirectories

%% Generate list_channel and list_path, two matrices with the same length, which specify for which channels the pixel intensities should be extracted.
%% In the current set up DAPI will be measured for the first listed cycle, the Red and the Green channel will be measured in all listed cycles. 
%% The sequence of channels will be as follows: Cycle_01_01 (Dapi), Cycle_01_02 (Green), Cycle_01_03 (Red), Cycle_02_02 (Green), Cycle_02_03 (Red), Cycle_03_02 (Green), Cycle_03_03 (Red), etc.  
list_channel = [channel_dapi; lin(repmat([channel_green;channel_red],1,numel(path_cycle)))];
list_path = path_cycle([1; lin(repmat([1:numel(path_cycle)],2,1))]); % 


%% close and open parpools
delete(gcp('nocreate'))
parpool('local',num_parallel)

%% Generate MPP Matrix by randomly selecting specified number of cells and extractiong MPP for each 2D pixel position within cell segmentation
[cells_MPPMatrix] = pick_MPP_from_objects(cells_metadata, sites_per_well, str_seg_object, list_path,list_channel);

%% Generate metadata information for each Multiplexed Pixel Profile
%% first column of mpp_metadata will contain the Cell ID the MPP came from
%% second column of mpp_metadata will contain the pixel ID, which is unique within one cell.
mpp_metadata = cell(numel(cells_MPPMatrix),1);
for ix = 1:numel(cells_MPPMatrix)
   
    % Generate cell id
    tmp_id = cell(numel(cells_MPPMatrix{ix}),1);
        
    % Generate pixel and cell ID
    tmp = NaN(size(cells_MPPMatrix{ix},1),2);
    tmp(:,1) = ix;
    tmp(:,2) = 1:size(tmp,1);
    tmp_id{ix} = tmp;
    
    mpp_metadata{ix} = cat(1,tmp_id{:});
    
end
mpp_metadata = cat(1,mpp_metadata{:});

%% Concatinate MPP Matrix
mpp_matrix = flatten(cells_MPPMatrix)';
mpp_matrix = cat(1,mpp_matrix{:});

%% get ready for saving
date_time = datestr(datetime('now'));
date_time = strrep(date_time,' ','_');
date_time = strrep(date_time,':','');

%% check that the folder is present, else make it
if ~isdir(path_output)
    mkdir(path_output);
end

%% prepare names and paths for saving
name_save = sprintf('MPPMatrix_%s',name_output);
name_save = strcat(date_time, '_',name_save,'.mat');
path_save = fullfile(path_output,name_save);

%% save matrices
save(path_save, 'cells_MPPMatrix', 'cells_metadata', 'mpp_matrix', 'mpp_metadata', 'list_channel', '-v7.3');
fprintf('Saved MatLab file as %s in %s \n', name_save, pwd)
fprintf('Finished\n')
end

%% Subfunctions
function [cells_pix_array] = pick_MPP_from_objects(metadata, sites_per_well, str_seg_object, list_path,list_channel)

cells_pix_array = cell(size(metadata,1),1);
cells_metadata = cell(size(metadata,1),1);
parfor ix = 1:size(metadata,1)   
    fprintf('Running %d of %d cells\n', ix, size(metadata,1))
    pix_of_cell = get_pix_info_per_multiplexed_cell(metadata(ix,:), sites_per_well, list_path, list_channel, str_seg_object);
    
    pix_problem = size(pix_of_cell,2) ~= numel(list_channel);
    if pix_problem
        while pix_problem
            fprintf('We have pix problem and entered while loop\n')
            pix_of_cell = get_pix_info_per_multiplexed_cell(metadata(ix,:), sites_per_well, list_path, list_channel, str_seg_object);
            pix_problem = size(pix_of_cell,2) ~= numel(list_channel);
        end
        [cells_pix_array{ix}] = pix_of_cell;
        fprintf('Exit while loop\n')
    else
        [cells_pix_array{ix}] = pix_of_cell;
    end
    cells_metadata{ix} = metadata(ix,:);
end
cells_metadata = cat(1,cells_metadata{:});
end

function [a] = get_pix_info_per_multiplexed_cell(cell_metadata, sites_per_well, list_path, list_channel, str_seg_object)

%% identify position of the cell in the dataset
row_names = {'A' 'B' 'C' 'D' 'E' 'F' 'G' 'H' 'I' 'J' 'K' 'L' 'M' 'N'};

position_row = row_names{cell_metadata(1)}; % identifies row
position_col = sprintf('%02d', cell_metadata(2)); % identifies column

site = rem(cell_metadata(3),sites_per_well); % identifies microscopy site
if site == 0
    position_site = sprintf('%03d',sites_per_well);
else
    position_site = sprintf('%03d',site);
end

position_cellid = cell_metadata(4); % identifies cell id

position_search = sprintf('.*%s%s.*F%s', position_row,position_col,position_site);


%% get both image and sementation
im_tiff_stack = cell(numel(list_channel),1);
im_cell_stack = cell(numel(list_channel),1);

channel_id = cellfun(@(x) sprintf('C0%d',x), mat2cell2(list_channel), 'UniformOutput', false); % change here, if channels are specified differently in your images

for ix = 1:numel(list_channel)
   
    path_tmp = list_path{ix};
    
    name_tiff = findfilewithregexpi(fullfile(path_tmp,'TIFF'),strcat(position_search, '.*', channel_id{ix})); % loads images
    im_tiff = imread_illumination_corrected(fullfile(path_tmp,'TIFF', name_tiff));
    
    name_seg_cells = findfilewithregexpi(fullfile(path_tmp,'SEGMENTATION'),strcat(position_search,'.*',str_seg_object)); % loads segmentations
    im_seg_cells = imread(fullfile(fullfile(path_tmp,'SEGMENTATION'), name_seg_cells));
    
    im_tiff_stack{ix} = im_tiff;
    im_cell_stack{ix} = im_seg_cells;
end

dims = size(im_tiff_stack{1});
im_tiff_stack_zscore = cellfun(@(x) reshape(lin(x),dims), im_tiff_stack, 'UniformOutput', false);

% make calculations
a = cellfun(@(x,y) x(y == position_cellid), im_tiff_stack_zscore, im_cell_stack, 'UniformOutput', false);
t = cellfun(@(x) numel(x), a);

if numel(unique(t)) > 1
else
    a = cat(2,a{:});
end
end

function corr_image = imread_illumination_corrected(strImage, shallCacheCorrection)
% reads images and applies illumination correction

if nargin < 2
    shallCacheCorrection = false;
end

strImageDir  = fileparts(strImage);
strBatchDir = strrep(strImageDir,'TIFF','BATCH');

iChannel = check_image_channel(strImage);
[matMeanImage, matStdImage, hasIlluminationCorrection] = getIlluminationReference(strBatchDir,iChannel,shallCacheCorrection);

if hasIlluminationCorrection == false
    error('could not find illumination correction')
else
    Image = double(imread(strImage));
    isLog = 1;
    corr_image = IllumCorrect(Image,matMeanImage,matStdImage,isLog);    
end

end
