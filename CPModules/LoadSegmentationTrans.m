function handles = LoadSegmentationTrans(handles)

% Help for the LoadSegmentationTrans module:
% LoadSegmentationTrans is heavily based on Load More Images
% Category: Object Processing

% $Revision: 1725 $

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%% new

%textVAR01 = Whats the name of the primary segmentation you want to import
%from trans?
%defaultVAR01 = Nuclei
%infotypeVAR01 = objectgroup indep
ObjectName_Primary = char(handles.Settings.VariableValues{CurrentModuleNum,1});

%textVAR02 = How do you want to call the primary trans segmentation?
%defaultVAR02 = TransNuclei
%infotypeVAR02 = objectgroup indep
TransObjectName_Primary = char(handles.Settings.VariableValues{CurrentModuleNum,2});

%textVAR03 = Whats the name of the secondary segmentation you want to import
%from trans? Ignore by keeping /
%defaultVAR03 = Cells
%infotypeVAR03 = objectgroup indep
ObjectName_Secondary = char(handles.Settings.VariableValues{CurrentModuleNum,3});

%textVAR04 = How do you want to call the secondary trans segmentation?
%defaultVAR04 = TransCells
%infotypeVAR04 = objectgroup indep
TransObjectName_Secondary = char(handles.Settings.VariableValues{CurrentModuleNum,4});

%textVAR05 = Whats the name of the tertiary segmentation you want to import
%from trans? Ignore by keeping /
%defaultVAR05 = Cytoplasm
%infotypeVAR05 = objectgroup indep
ObjectName_Tertiary = char(handles.Settings.VariableValues{CurrentModuleNum,5});

%textVAR06 = How do you want to call the tertiary trans segmentation?
%defaultVAR06 = TransCytoplasm
%infotypeVAR06 = objectgroup indep
TransObjectName_Tertiary = char(handles.Settings.VariableValues{CurrentModuleNum,6});

%textVAR07 = Which channel forms the basis of the filenames?
%infotypeVAR07 = imagegroup
OrigImageName = char(handles.Settings.VariableValues{CurrentModuleNum,7});
%inputtypeVAR07 = popupmenu

%pathnametextVAR08 = From which reference acquistion should objects be imported?
%defaultVAR08 = L:\Data\Users\Gabriele\Multiplexing\20160529_MPSimulation01\20160531_MPSimulation01_stain02
Pathname = char(handles.Settings.VariableValues{CurrentModuleNum,8});

%%%VariableRevisionNumber = 1

%%%%%%%%%%%%%%%%%%%%%%%%
%%% INITIALIZATION   %%%
%%%%%%%%%%%%%%%%%%%%%%%%

strTransPlate = Pathname;
if ~any(fileattrib(strTransPlate))
    error('Could not find reference plate');
end

if handles.Current.SetBeingAnalyzed == 1
    makeExternalBuffer(handles, strTransPlate, ObjectName_Primary, ObjectName_Secondary, ObjectName_Tertiary)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY CALCULATIONS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% Determines which cycle is being analyzed.
SetBeingAnalyzed = handles.Current.SetBeingAnalyzed;

%%% preproduce trans object names
TransObjectName = cell(1);
if hasObjectBeenDefined(ObjectName_Primary)
    TransObjectName{end+1,1} = ['Trans' ObjectName_Primary];
end
if hasObjectBeenDefined(ObjectName_Secondary)
    TransObjectName{end+1,1} = ['Trans' ObjectName_Secondary];
end
if hasObjectBeenDefined(ObjectName_Tertiary)
    TransObjectName{end+1,1} = ['Trans' ObjectName_Tertiary];
end
TransObjectName = TransObjectName(2:end,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FIRST CYCLE FILE HANDLING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% Extracting the list of files to be analyzed occurs only the first time
%%% through this module.

% get segmentation paths
[~, strTransSegmentationFolder] = getSubFoldersFromTransPlate(strTransPlate);

% check wheter you are overwriting orig images with trans images
if SetBeingAnalyzed == 1
    
    for i = 1:length(TransObjectName)
        if isfield(handles.Pipeline,TransObjectName{i})
            error(['Image processing was cancelled in the ', ModuleName, ' module because you are trying to load two sets of images with the same name (e.g. OrigBlue). The last set loaded will always overwrite the first set and make it obselete. Please remove one of these modules.']);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LOADING Segmentation and save it %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

drawnow
%Load
tmp = cell(numel(TransObjectName),1);
[strCorrespondingSegmentationPrimary_trans, couldFindSameSite_segmentation_Primary, strCorrespondingSegmentationSecondary_trans, couldFindSameSite_segmentation_Secondary, strCorrespondingSegmentationTertiary_trans, couldFindSameSite_segmentation_Tertiary] = ...
    getFileNameViaReferenceFile(handles, OrigImageName, ObjectName_Primary, ObjectName_Secondary, ObjectName_Tertiary);

if couldFindSameSite_segmentation_Primary
    fp = fullfile(strTransSegmentationFolder, strCorrespondingSegmentationPrimary_trans);
    matSegmentationImage = double(imread(fp));
    
    % save object into handles
    fieldname = ['Segmented',TransObjectName_Primary];
    handles.Pipeline.(fieldname) = matSegmentationImage;
    fieldname = ['UneditedSegmented',TransObjectName_Primary];
    handles.Pipeline.(fieldname) = matSegmentationImage;
    fieldname = ['SmallRemovedSegmented',TransObjectName_Primary];
    handles.Pipeline.(fieldname) = matSegmentationImage;
    
    column = find(strcmp(handles.Measurements.Image.ObjectCountFeatures,TransObjectName_Primary));
    if isempty(column)
        handles.Measurements.Image.ObjectCountFeatures(end+1) = {TransObjectName_Primary};
        column = length(handles.Measurements.Image.ObjectCountFeatures);
    end
    ObjCount = max(matSegmentationImage(:));
    handles.Measurements.Image.ObjectCount{handles.Current.SetBeingAnalyzed}(1,column) = ObjCount;
    RP = regionprops(matSegmentationImage,'Centroid');
    Centroid = cat(1,RP.Centroid);
    handles.Measurements.(TransObjectName{1}).LocationFeatures = {'CenterX','CenterY'};
    handles.Measurements.(TransObjectName{1}).Location(handles.Current.SetBeingAnalyzed) = {Centroid};
    
    % store image for display
    tmp{1} = matSegmentationImage;
else
    error('No corresponding image could be found')
end

if couldFindSameSite_segmentation_Secondary
    fp = fullfile(strTransSegmentationFolder, strCorrespondingSegmentationSecondary_trans);
    matSegmentationImage = double(imread(fp));
    
    % save object into handles
    fieldname = ['Segmented',TransObjectName_Secondary];
    handles.Pipeline.(fieldname) = matSegmentationImage;
    fieldname = ['UneditedSegmented',TransObjectName_Secondary];
    handles.Pipeline.(fieldname) = matSegmentationImage;
    fieldname = ['SmallRemovedSegmented',TransObjectName_Secondary];
    handles.Pipeline.(fieldname) = matSegmentationImage;
    
     column = find(strcmp(handles.Measurements.Image.ObjectCountFeatures,TransObjectName_Secondary));
    if isempty(column)
        handles.Measurements.Image.ObjectCountFeatures(end+1) = {TransObjectName_Secondary};
        column = length(handles.Measurements.Image.ObjectCountFeatures);
    end
    ObjCount = max(matSegmentationImage(:));
    handles.Measurements.Image.ObjectCount{handles.Current.SetBeingAnalyzed}(1,column) = ObjCount;
    RP = regionprops(matSegmentationImage,'Centroid');
    Centroid = cat(1,RP.Centroid);
    handles.Measurements.(TransObjectName{2}).LocationFeatures = {'CenterX','CenterY'};
    handles.Measurements.(TransObjectName{2}).Location(handles.Current.SetBeingAnalyzed) = {Centroid};
    
    % store image for display   
    tmp{2} = matSegmentationImage;
else
    error('No corresponding image could be found')
end

if couldFindSameSite_segmentation_Tertiary
    fp = fullfile(strTransSegmentationFolder, strCorrespondingSegmentationTertiary_trans);
    matSegmentationImage = double(imread(fp));
    
    fieldname = ['Segmented',TransObjectName_Tertiary];
    handles.Pipeline.(fieldname) = matSegmentationImage;
    fieldname = ['UneditedSegmented',TransObjectName_Tertiary];
    handles.Pipeline.(fieldname) = matSegmentationImage;
    fieldname = ['SmallRemovedSegmented',TransObjectName_Tertiary];
    handles.Pipeline.(fieldname) = matSegmentationImage;
    
         column = find(strcmp(handles.Measurements.Image.ObjectCountFeatures,TransObjectName_Tertiary));
    if isempty(column)
        handles.Measurements.Image.ObjectCountFeatures(end+1) = {TransObjectName_Tertiary};
        column = length(handles.Measurements.Image.ObjectCountFeatures);
    end
    ObjCount = max(matSegmentationImage(:));
    handles.Measurements.Image.ObjectCount{handles.Current.SetBeingAnalyzed}(1,column) = ObjCount;
    RP = regionprops(matSegmentationImage,'Centroid');
    Centroid = cat(1,RP.Centroid);
    handles.Measurements.(TransObjectName{3}).LocationFeatures = {'CenterX','CenterY'};
    handles.Measurements.(TransObjectName{3}).Location(handles.Current.SetBeingAnalyzed) = {Centroid};
    
    % store image for display   
    tmp{3} = matSegmentationImage;
else
    error('No corresponding image could be found')
end

%%%%%%%%%%%%%%%%%%%%%%%
%%% DISPLAY RESULTS %%%
%%%%%%%%%%%%%%%%%%%%%%%
drawnow

ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
if any(findobj == ThisModuleFigureNumber)
    if CPisHeadless == false
        %%% Activates the appropriate figure window.
        CPfigure(handles,'Image',ThisModuleFigureNumber);
        for ix = 1:numel(TransObjectName)
            subplot(1,numel(TransObjectName), ix)
            % RGB color
            ColoredLabelMatrixImage = CPlabel2rgb(handles,tmp{ix});
            CPimagesc(ColoredLabelMatrixImage,handles);
            axis image
            %imagesc(ColoredLabelMatrixImage);
            
            title(sprintf('Loaded %s segmentation , cycle # %d',TransObjectName{ix},handles.Current.SetBeingAnalyzed));
        end
    end
end

end


%% Subfunctions
% to find corresponding image in trans
function makeExternalBuffer(handles, strTransPlate, ObjectName_Primary, ObjectName_Secondary, ObjectName_Tertiary)
% note that preferentially all data would be stored in pipeline. However
% storing custom fields in handles.pipelines can lead to various errors.
% While storing data in handles.measurements is technically possible there
% would be various mistakes occuring if order of sites becomes rearranged
% in batch files. Since there will always only be one segmentation of a
% given name, there can be a buffer, that links to that name, and which is
% overwritten, if a new pipeline is made


[~, SegmentationFolder_trans] = getSubFoldersFromTransPlate(strTransPlate);

%% new
% Segmentations (primary)
filterForSegmentationsOfTrans = ['Segmented' ObjectName_Primary '\.'];
eB.strSegmentations_trans = getFilesAndDirectories(SegmentationFolder_trans, filterForSegmentationsOfTrans);
if isempty(eB.strSegmentations_trans)
    error(['Could not find Segementations for ' ObjectName_Primary]);
end

[eB.Row_segmentations_trans, eB.intColumn_segmentations_trans, eB.intImagePosition_segmentations_trans] = cellfun(@(x) MetaFromImageName(x), eB.strSegmentations_trans, 'UniformOutput', true);

% Segmentations (secondary)
if hasObjectBeenDefined(ObjectName_Secondary) == true
    filterForSegmentationsOfTrans_Secondary = ['Segmented' ObjectName_Secondary '\.'];
    eB.strSegmentationsSecondary_trans = getFilesAndDirectories(SegmentationFolder_trans, filterForSegmentationsOfTrans_Secondary);
    [eB.Row_segmentationsSecondary_trans, eB.intColumn_segmentationsSecondary_trans, eB.intImagePosition_segmentationsSecondary_trans] = ...
        cellfun(@(x) MetaFromImageName(x), eB.strSegmentationsSecondary_trans, 'UniformOutput', true);
    
    if isempty(eB.strSegmentationsSecondary_trans)
        error(['Could not find Segementations for ' ObjectName_Secondary]);
    end
end

% Segmentations (teritiary)
if hasObjectBeenDefined(ObjectName_Tertiary) == true
    filterForSegmentationsOfTrans_Tertiary = ['Segmented' ObjectName_Tertiary '\.'];
    eB.strSegmentationsTertiary_trans = getFilesAndDirectories(SegmentationFolder_trans, filterForSegmentationsOfTrans_Tertiary);
    [eB.Row_segmentationsTertiary_trans, eB.intColumn_segmentationsTertiary_trans, eB.intImagePosition_segmentationsTertiary_trans] = ...
        cellfun(@(x) MetaFromImageName(x), eB.strSegmentationsTertiary_trans, 'UniformOutput', true);
    
    if isempty(eB.strSegmentationsSecondary_trans)
        error(['Could not find Segementations for ' ObjectName_Tertiary]);
    end
end


% save name
outDir = handles.Current.DefaultOutputDirectory;
ex = ['SimpleMulti_' 'SegmentationTrans' '.mat'];
strBufferFile = fullfile(outDir, ex );


if any(fileattrib(strBufferFile))
    [~, ex] = fileparts(strBufferFile);
    fprintf([ex ' already exists. It will be overwritten with current one.\n']);
end

save(strBufferFile, 'eB');
end
function [TiffFolder_trans, SegmentationFolder_trans] = getSubFoldersFromTransPlate(strTransPlate)


if any(strfind(strTransPlate, [filesep 'TIFF']))
    error('Reference directory must refer to a plate folder, not the TIFF folder');
end

if any(strfind(strTransPlate, [filesep 'SEGMENTATION']))
    error('Reference directory must refer to a plate folder, not the SEGMENTATION folder');
end

TiffFolder_trans = fullfile(strTransPlate, 'TIFF');
if ~any(fileattrib(TiffFolder_trans))
    error('Could not find TIFF folder of other plate');
end

SegmentationFolder_trans = fullfile(strTransPlate, 'SEGMENTATION');
if ~any(fileattrib(SegmentationFolder_trans))
    error('Could not find SEGMENTATION folder of other plate');
end


end
function [strCorrespondingSegmentationPrimary_trans, couldFindSameSite_segmentation_Primary, strCorrespondingSegmentationSecondary_trans, couldFindSameSite_segmentation_Secondary, strCorrespondingSegmentationTertiary_trans, couldFindSameSite_segmentation_Tertiary] = ...
    getFileNameViaReferenceFile(handles, OrigImageName, ObjectName_Primary, ObjectName_Secondary, ObjectName_Tertiary)
% Avoid repeated, independent access to buffer file by reading segmentation
% and images after same loading event. Note that the code thus becomes more
% ugly and less readable.


SetBeingAnalyzed = handles.Current.SetBeingAnalyzed;
% size(handles.Pipeline.(['Filename' OrigImageName])) %[GG20160419] old
% SetBeingAnalyzed %[GG20160419] old
ImageName_cis = handles.Pipeline.(['Filename' OrigImageName]){SetBeingAnalyzed}; %[GG20160419] old
% ImageName_cis = handles.Pipeline.(['Filename' OrigImageName]){1}; %[GG20160419] new

outDir = handles.Current.DefaultOutputDirectory;
ex = ['SimpleMulti_' 'SegmentationTrans' '.mat'];
strBufferFile = fullfile(outDir, ex );

if ~any(fileattrib(strBufferFile))
    error('Could not find buffer file for object, that should be created during first cycle');
else
    load(strBufferFile);
end

[Row_cis, intColumn_cis, intImagePosition_cis] = MetaFromImageName(ImageName_cis);

% get corresponding segmentation from trans (primary)
correspondsToSameSite_segmentation_Primary = eB.Row_segmentations_trans == Row_cis & eB.intColumn_segmentations_trans == intColumn_cis & eB.intImagePosition_segmentations_trans == intImagePosition_cis;
if sum(correspondsToSameSite_segmentation_Primary) == 0
    couldFindSameSite_segmentation_Primary = false;
    strCorrespondingSegmentationPrimary_trans = '';
elseif sum(correspondsToSameSite_segmentation_Primary) == 1;
    couldFindSameSite_segmentation_Primary = true;
    strCorrespondingSegmentationPrimary_trans = eB.strSegmentations_trans{correspondsToSameSite_segmentation_Primary};
else
    error(['Could not unambiguously find '  ObjectName_Primary ' of other dataset. Please set more stringent filters.']);
end

% get corresponding segmentation from trans (secodnary)
correspondsToSameSite_segmentation_Secondary = eB.Row_segmentationsSecondary_trans == Row_cis & eB.intColumn_segmentationsSecondary_trans == intColumn_cis & eB.intImagePosition_segmentationsSecondary_trans == intImagePosition_cis;
if sum(correspondsToSameSite_segmentation_Secondary) == 0
    couldFindSameSite_segmentation_Secondary = false;
    strCorrespondingSegmentationSecondary_trans = '';
elseif sum(correspondsToSameSite_segmentation_Secondary) == 1;
    couldFindSameSite_segmentation_Secondary = true;
    strCorrespondingSegmentationSecondary_trans = eB.strSegmentationsSecondary_trans{correspondsToSameSite_segmentation_Secondary};
else
    error(['Could not unambiguously find '  ObjectName_Secondary ' of other dataset. Please set more stringent filters.']);
end

% get corresponding segmentation from trans (tertiary)
correspondsToSameSite_segmentation_Tertiary = eB.Row_segmentationsTertiary_trans == Row_cis & eB.intColumn_segmentationsTertiary_trans == intColumn_cis & eB.intImagePosition_segmentationsTertiary_trans == intImagePosition_cis;
if sum(correspondsToSameSite_segmentation_Tertiary) == 0
    couldFindSameSite_segmentation_Tertiary = false;
    strCorrespondingSegmentationTertiary_trans = '';
elseif sum(correspondsToSameSite_segmentation_Tertiary) == 1;
    couldFindSameSite_segmentation_Tertiary = true;
    strCorrespondingSegmentationTertiary_trans = eB.strSegmentationsTertiary_trans{correspondsToSameSite_segmentation_Tertiary};
else
    error(['Could not unambiguously find '  ObjectName_Tertiary ' of other dataset. Please set more stringent filters.']);
end

end
function [CurrFileList, CurrDirectoryList] = getFilesAndDirectories(strInputDir,strFilter)
CurrFileList = [];
CurrDirectoryList = [];
if nargin<2
    doFilterStep = false;
else
    doFilterStep = true;
end

CurrFileListImport = CPdir(strInputDir);
CurrFileListImport = struct2cell(CurrFileListImport);
f = ~cell2mat(CurrFileListImport(2,:));
if any(f)
    CurrFileList = CurrFileListImport(1,f)';
end
if any(~f)
    CurrDirectoryList = CurrFileListImport(1,~f);
end

if doFilterStep == true && ~isempty(CurrFileList)
    f = cell2mat(cellfun(@(x) ~isempty(regexp(x,strFilter,'once')), CurrFileList,'UniformOutput',false));
    CurrFileList = CurrFileList(f);
    
    f = cell2mat(cellfun(@(x) ~isempty(regexp(x,strFilter,'once')), CurrDirectoryList,'UniformOutput',false));
    CurrDirectoryList = CurrDirectoryList(f);
end

f = ismember(CurrDirectoryList,{'.';'..';'.duc'});
if any(f)
    CurrDirectoryList = CurrDirectoryList(~f);
end

end
function ObjectIsDefined = hasObjectBeenDefined(ObjectName)
ObjectIsDefined = ~isequal(ObjectName,'/'); % no secondary specified;
end
function [intRow, intColumn, intImagePosition, intTimepoint, intZstackNumber, intChannelNumber, strMicroscopeType, strWellName, intActionNumber] = MetaFromImageName(strImageName)
% METAFROMIMAGES collects metainformation about image acquistion from the
% file name. It serves as a hub for multiple different functions,
% previously developed within the lab, where features are derived from
% parsing the file names. STRIMAGENAME is the name of the file. You might
% add custom modifications to extract your metainformation of interest.
%
%%%%%%%% IMPORTANT: MAKE YOUR CUSTOM ADJUSTMENTS %%%%%%%%%%%%%%
% a) Use regular expressions (see matlab help) to obtain the information. 
% This might look something like
% strChannelMatch = regexp(strImageName, '_([^_]{3})_(T\d{4})F(\d{3})L(\d{2})A(\d{2})Z(\d{2})C(\d{2})', 'Tokens');
% b) Some of the output arguments are actually not required for Spot detection, 
% so you can set them to default value (indicated by * in
% description)
%
%   
%   Authors:
%   Nico Battich
%   Thomas Stoeger
%   Lucas Pelkmans
%
% Battich et al., 2013.
% Website: http://www.imls.uzh.ch/research/pelkmans.html
% *************************************************************************
%
%   NAME                TYPE       
%   intChannelNumber    double      
%   Number describing the color, eg. 1 for blue, 2 for green, 3 for red, 4 
%   for far red,.... (might be custom)
%
%   intZstackNumber     double      
%   Number describing the Z-plane. NaN if not suitable, otherwise: start 
%   with 1 and increment by 1 corresponding to subsequent stage positions.
%   *(1) *only if no 3D analysis
%   
%   intActionNumber     double      
%   Number describing the action. eg.: if multiple channels were acquired 
%   at the same time by parallel optics/cameras. NaN if not suitable. *(1)
%
%   intImagePosition    double
%   Number describing the acquisition site within a single well. 
%
%   strMicroscopeType   character
%   Name of microscope. Mostly used for better error messages. *('myMic')
%
%   intRow      double
%   Row of well. Row A is 1. Row B is 2. Row C is 3.....
%
%   intColumn       double
%   Column of well.
%
%   strWellName   double
%   Name of Well accoding with combination of letter (for rows) and number
%   (for columns) eg. A05 . *('A01')
%
%   int Timepoint   double *
%   Number indicating the timepoint of a time series. Starts with 1 and
%   increments by 1 for subsequent time points. *(1)


% Add security and informative error message for publication
if isempty(which('check_image_channel')) || isempty(which('check_image_channel')) || isempty(which('filterimagenamedata'))
    error('Sorry. You have to do some coding: You have to write custom functions, which obtain the metadata (such as image channel, Position in Z-stack...) from the filename (or alternative source). Place them in MetaFromImageName.m')
else
    % Call functions, which anlyse image file names for metadata
    [intChannelNumber,intZstackNumber,intActionNumber] = check_image_channel(strImageName);
    [intImagePosition,strMicroscopeType] = check_image_position(strImageName);
    [intRow, intColumn, strWellName, intTimepoint] = filterimagenamedata(strImageName);
end

end
function [intChannelNumber,intZstackNumber,intActionNumber] = check_image_channel(strImageName)

    if nargin==0
         strImageName = 'blablabl_NIKON_t0001_C21_s1_z0001_w01_ledGFP.png';
    end

    
    intChannelNumber = NaN;
    intZstackNumber = NaN;
    intActionNumber = NaN;
    
    % CW
    strNomenclature1 = regexp(strImageName,'f\d\dd\d.(tif|png)','Match');
    strNomenclature1a = regexp(strImageName,'f\dd\d.(tif|png)','Match');
    strNomenclature2 = regexp(strImageName,'_w\d\d\d.(tif|png)','Match');
    strNomenclature3 = regexp(strImageName,' - n\d{2,}.(tif|png)','Match');
    
    % NIKON
    strNomenclaturePre4 = regexp(strImageName,'NIKON.*_t\d{1,}(_z\d{1,}_|_)[A-Z]\d{2,3}_s\d{1,}_w\d{1,}[^\.]*\.(tiff?|png)','Match');

    % MD MICROEXPRESS
    strNomenclature4 = regexp(strImageName,'_\w\d\d_s\d{1,}_w\d','Match');
    strNomenclature4a = regexp(strImageName,'_\w\d\d_s\d{1,}[A-Z0-9\-]{36}\.(tif|png)','Match');
    
    % CV7K
    strNomenclature5 = regexp(strImageName, ...
        '_([^_]{3})_(T\d{4})(F\d{3})(L\d{2})(A\d{2})(Z\d{2})(C\d{2})', 'Match');

    % VisiScope in slide scan mode
    strNomenclature6 = regexp(strImageName,'_s\d{04}_r\d{02}_c\d{02}_[^_]+_C\d{02}','Match');
    
    strChannelWavelengths = {'_w460','_w530','_w595','_w685'};
    strChannelDescriptions = {'DAPI','FITC','TRITC','CY5'};

    if not(isempty(strNomenclature1))
        %%% CELLWORX
        strChannelMatch = regexp(strImageName,'\d.(tif|png)','Match');
        intChannelNumber = str2double(strrep(strrep(char(strChannelMatch{1}),'.tif',''),'.png',''))+1;
    elseif not(isempty(strNomenclature2))
        %%% CELLWORX        
        strChannelMatch = regexp(strImageName,'_w\d\d\d.(tif|png)','Match');
        intChannelNumber = find(~cellfun('isempty',strfind(strChannelWavelengths,char(strrep(strrep(strChannelMatch,'.tif',''),'.png','')))));
   elseif not(isempty(strNomenclature1a))
        %%% CELLWORX       
        strChannelMatch = regexp(strImageName,'\d.(tif|png)','Match');
        intChannelNumber = str2double(strrep(strrep(char(strChannelMatch{1}),'.tif',''),'.png',''))+1;
   elseif not(isempty(strNomenclature3))
        %%% BD-PATHWAY
        for i = 1:length(strChannelDescriptions)
            if ~isempty(findstr(upper(strChannelDescriptions{i}),upper(strImageName)))
                intChannelNumber = i;
            end
        end
    elseif not(isempty(strNomenclaturePre4))
        %%% NIKON
        %strChannelMatch = regexp(strImageName,'NIKON.*_\w\d\d_s\d{3,}_w(\d\d)','Tokens');
        strChannelMatch = regexpi(strImageName,'NIKON.*_t\d{1,}(_z\d{1,}_|_)[A-Z]\d{2,3}_s(\d{1,})_w(\d{1,})[^\.]*\.(tiff?|png)','Tokens');
        intChannelNumber = str2double(strChannelMatch{1}(3));
        strZstack = char(strChannelMatch{1}(1));
        if ~strcmp(strZstack,'_') % does filename actually contain a z-stack info?
            intZstackNumber = str2double(strZstack(3:end-1));
        else
            intZstackNumber = NaN;
        end
    elseif not(isempty(strNomenclature4))
        %%% MD
        strChannelMatch = regexp(strImageName,'_\w\d\d_s\d{1,}_w(\d)','Tokens');
        intChannelNumber = str2double(strChannelMatch{1});
    elseif not(isempty(strNomenclature4a))
        %%% MD - if there is only one channel present
        intChannelNumber = 1;
    elseif not(isempty(strNomenclature5))    
        %%% CV7K - we have a match against "Yokogawa" filenaming tail 
        strChannelMatch = regexp(strImageName, '_([^_]{3})_(T\d{4})F(\d{3})L(\d{2})A(\d{2})Z(\d{2})C(\d{2})', 'Tokens');
        intChannelNumber = str2double(strChannelMatch{1}(7));
        intZstackNumber = str2double(strChannelMatch{1}(6));
        intActionNumber = str2double(strChannelMatch{1}(5));
    else
        warning('iBRAIN:check_image_channel','unknown file name %s',strImageName)
    end
end
function [intImagePosition,strMicroscopeType] = check_image_position(strImageName)
% help for check_image_position()
% BS, 082015
% usage [intImagePosition,strMicroscopeType] = check_image_position(strImageName)
%
% possible values for strMicroscopeType are 'CW', 'BD', and 'MD'

    if nargin == 0
%         strImageName = '061224_YD_50K_KY_P1_1_1_H07f00d0.png';
    % strImageName = '090313_SV40_FAKKO_GLY_GLYGGPP_GM1_noGM1_A01_01_w460_SegmentedCells.png';
    
%    strImageName = '070314_TDS_VSV_50K_P2_1_B02f0d0_SegmentedCells.png';
        
%         strImageName = '070420_Tf_KY_P2_20x_C10_19_w460.tif'
%         strImageName = 'Dapi - n000000.tif'
%         strImageName = 'DAPI_p53SLS - n000000.png';    
%         strImageName = '040423_frank_si_VSV_B02_25_w460.tif';  
%         strImageName = 'Y:\Data\Users\Jean-Philippe\p53hitvalidation\ko_p53Pro72_plate1_triplicate1\TIFF\Well H12\Dapi_p53SLS - n000000.tif';
%         strImageName = '080611-olivia-1_A10_s1_w12E22EFEB-B167-43E0-A05F-997CCA19728A.tif'
%         strImageName = '080815-Frank-VSV-pH-dyn_B02_s1_w11BBF4034-97B9-4912-9DA5-6FBAF05BA7E4.tif'
%         strImageName = '081008_VV_rescreen_CP078-1aa_K04_8_w530.png';
%         strImageName = '2008-08-14_HPV16_batch1_CP001-1ea_A20_9_w530.tif'
%         strImageName = 'BAC-siRNA-Optimisation_C01_s3CB0B5EFE-CA88-49D1-B8B8-2115D7B91A6F.png'
         strImageName = 'RDUP20120522-CNX_K04_s6_w210E554F9-A2B4-4D39-969A-C1216D5C5893.png';        
    end

    strMicroscopeType = '';
    intImagePosition = 1;
    
    % CW
    strNomenclature1 = regexp(strImageName,'f\d\dd\d.(tif|png)','Match');
    strNomenclature1a = regexp(strImageName,'f\dd\d.(tif|png)','Match'); 
    strNomenclature1b = regexp(strImageName,'f\d\dd\d','Match');
    strNomenclature1c = regexp(strImageName,'f\dd\d','Match');
    strNomenclature2 = regexp(strImageName,'_w\d\d\d.(tif|png)','Match'); 
    strNomenclature3 = regexp(strImageName,' - n\d{2,}.(tif|png)','Match');
    
    % NIKON
    strNomenclaturePre4 = regexp(strImageName,'NIKON.*_t\d{1,}(_z\d{1,}_|_)[A-Z]\d{2,3}_s\d{1,}_w\d{1,}[^\.]*\.(tiff?|png)','Match');
    
    % MD MICROEXPRESS
    strNomenclature4 = regexp(strImageName,'_\w\d\d_s\d{1,}_w\d','Match');    
    % MD with only one channel matches this but not the previous
    strNomenclature4a = regexp(strImageName,'_\w\d\d_s\d{1,}[A-Z0-9\-]{36}','Match');
    strNomenclature4b = regexp(strImageName,'_\w\d\d_\d{1,}_w\d','Match');    
    
    % CV7K
    %%[NB] here is a fix from the first regexp below. It was not compatible
    %%with the iBrainTrackerV1.
    strNomenclature5 = regexp(strImageName, ...
        '_([^_]{3})_(T\d{4})(F\d{3})(L\d{2})(A\d{2})(Z\d{2})(C\d{2})', 'Match');
    %    strNomenclature5 = regexp(strImageName, ...
    %    '_([^_]{3})_(T\d{4})(F\d{3})(L\d{2})(A\d{2})(Z\d{2})(C\d{2})\.(tif|png)$', 'Match');

    % VisiScope in slide scan mode
    strNomenclature6 = regexp(strImageName,'_s\d{04}_r\d{02}_c\d{02}_[^_]+_C\d{02}','Match');
    
    % fallback
	strNomenclature7 = regexp(strImageName,'_s\d{1,}_','Match');
    
    if not(isempty(strNomenclature1))
        strMicroscopeType = 'CW';
        strImagePosition = regexp(strImageName,'f\d\dd\d.(tif|png)','Match');
        strImagePosition = strImagePosition{1,1}(1,1:3);
        strImagePosition = strrep(strImagePosition,'f','');
        intImagePosition = str2double(strImagePosition)+1;        
    elseif not(isempty(strNomenclature2))
        strMicroscopeType = 'CW';        
        strImagePosition = regexp(strImageName,'_\d\d_w\d\d\d.(tif|png)','Match');
        if not(isempty(strImagePosition))
            intImagePosition = str2double(strImagePosition{1,1}(1,2:3));
        else
            strImagePosition = regexp(strImageName,'_\d_w\d\d\d.(tif|png)','Match');
            intImagePosition = str2double(strImagePosition{1,1}(1,2));
        end
    elseif not(isempty(strNomenclature1a))
        strMicroscopeType = 'CW';        
        strImagePosition = regexp(strImageName,'f\dd\d.(tif|png)','Match');
        strImagePosition = strImagePosition{1,1}(1,1:2);
        strImagePosition = strrep(strImagePosition,'f','');
        intImagePosition = str2double(strImagePosition)+1;        
    elseif not(isempty(strNomenclature1b))
        strMicroscopeType = 'CW';
        strImagePosition = regexp(strImageName,'f\d\dd\d','Match');
        strImagePosition = strImagePosition{1,1}(1,1:3);
        strImagePosition = strrep(strImagePosition,'f','');
        intImagePosition = str2double(strImagePosition)+1;     
    elseif not(isempty(strNomenclature1c))
        strMicroscopeType = 'CW';
        strImagePosition = regexp(strImageName,'f\dd\d','Match');
        strImagePosition = strImagePosition{1,1}(1,2);
        intImagePosition = str2double(strImagePosition)+1;             
    elseif not(isempty(strNomenclature3))
        strMicroscopeType = 'BD';        
        strImagePosition = regexpi(strImageName,' - n(\d{2,}).(tif|png)','Tokens');
        if ~isempty(strImagePosition)
            intImagePosition = str2double(strImagePosition{1}{1})+1;
        end   
    elseif not(isempty(strNomenclaturePre4))
        strMicroscopeType = 'NIKON';
        strImagePosition = regexpi(strImageName,'NIKON.*_t\d{1,}(_z\d{1,}_|_)[A-Z]\d{2,3}_s(\d{1,})_w(\d{1,})[^\.]*\.(tiff?|png)','Tokens');
        intImagePosition = str2double(char(strImagePosition{1}(2)));
    elseif not(isempty(strNomenclature4))
        strMicroscopeType = 'MD';        
        strImagePosition = regexpi(strImageName,'_\w\d\d_s(\d{1,})_w\d','Tokens');        
        intImagePosition = str2double(strImagePosition{1});
    elseif not(isempty(strNomenclature4a))
        strMicroscopeType = 'MD';
        strImagePosition = regexpi(strImageName,'_\w\d\d_s(\d{1,})[A-Z0-9\-]{36}','Tokens');  
        intImagePosition = str2double(strImagePosition{1});        
    elseif not(isempty(strNomenclature4b))
        strMicroscopeType = 'MD';
        strImagePosition = regexpi(strImageName,'_\w\d\d_(\d{1,})_w\d','Tokens');        
        intImagePosition = str2double(strImagePosition{1});
    elseif not(isempty(strNomenclature5))
        %%% CV7K - we have a match against "Yokogawa" filenaming tail
        strMicroscopeType = 'CV7K';
        %%[NB] I change this to be compatible with tracker as before
        strImagePosition = regexp(strImageName, '_([^_]{3})_(T\d{4})F(\d{3})(L\d{2})(A\d{2})(Z\d{2})C(\d{2})', 'Tokens');
        %strImagePosition = regexp(strImageName, '_([^_]{3})_(T\d{4})F(\d{3})(L\d{2})(A\d{2})(Z\d{2})C(\d{2})\.(tif|png)$', 'Tokens');
        intImagePosition = str2double(strImagePosition{1}(3));
    elseif not(isempty(strNomenclature7)) 
        strMicroscopeType = 'MD';
        strImagePosition = regexpi(strImageName,'_s(\d{1,})_','Tokens');
        intImagePosition = str2double(strImagePosition{1});
    else         
        error('unknown file name %s',strImageName)
    end
end
function [intRow, intColumn, strWellName, intTimepoint] = filterimagenamedata(strImageName)
% see HELP of METAFROMIMAGENAME for more information

    intRow = NaN;
    intColumn = NaN;
    strWellName = NaN;
    intTimepoint  = 1;
    
    if nargin == 0

        strImageName = 'blablabl_NIKON_t0001_C21_s1_z0001_w01_ledGFP.png';
    end

    strWellName = char(strrep(regexp(strImageName,'_[A-Z]\d\d_','Match'),'_',''));
    strWellName2 = char(strrep(regexp(strImageName,'_[A-Z]\d\df','Match'),'_',''));
    strWellName3 = char(strrep(regexp(strImageName,'Well [A-Z]\d{2,}','Match'),'_',''));
    
    % MD MICROEXPRESS
    strNomenclature4 = regexp(strImageName,'_\w\d\d_s\d{1,}_w\d','Match');    
        
    
    % In case of time-points, perhaps override the well position and assign
    % different time points artificial well positions? (temp hack really)
    % match timepoint
    cellstrTimepoint = regexp(strImageName,'_[tT]_?(\d{1,})_?','Tokens');  
    if not(isempty(cellstrTimepoint))
        intTimepoint = str2double(cellstrTimepoint{1});
    else
        intTimepoint = 1;
    end
    
    % This was a timepoint hack, but disabled, was giving problems...
%     if not(isempty(strWellNameTime))
%         matFakeRows = lin(repmat(1:12,24,1));
%         matFakeColumns = lin(repmat([1:24],12,1)');
%         int_time = str2double(char(strWellNameTime{1}));
%         intRow = matFakeRows(int_time);
%         intColumn = matFakeColumns(int_time);
%         strWellName = sprintf('%s%02d',char(intRow+64),intColumn);
%         fprintf('%s: time point detected, faking 384 well position (row=%02d, col=%02d = %s) data based on time point t=%03d\n',mfilename,intRow,intColumn,strWellName,int_time)
%     end

% Note that we shold always take the last match... as user might input
% something before that would/could look like a well.

    
    if not(isempty(strWellName))
        intRow=double(strWellName(end,1))-64;
        intColumn=str2double(strWellName(end,2:3));
        strWellName=strWellName(end,1:end);

    elseif not(isempty(strWellName2))
        intRow=double(strWellName2(end,1))-64;
        intColumn=str2double(strWellName2(end,2:3));

        strWellName=strWellName2(end,1:end-1);

    elseif not(isempty(strWellName3))
        strImageData = regexp(strImageName,'Well ([A-Z])(\d{2,})','Tokens');
        intRow=double(strImageData{1}{1})-64;
        intColumn=str2double(strImageData{1}{2});
        strWellName=sprintf('%s%.02d',strImageData{1}{1},intColumn);

    elseif not(isempty(strNomenclature4))
        %%% MD
        strImageData = regexpi(strImageName,'_(\w)(\d\d)_s','Tokens');
        intRow=double(strImageData{1})-64;
        intColumn=str2double(strImageData{2});
        strWellName=[strImageData{1},strImageData{2}];
    else
        intRow = NaN;
        intColumn = NaN;
        strWellName = NaN;

        warning('filterimagenamedata: unable to get well data from image name %s',strImageName)
    end

	
end

