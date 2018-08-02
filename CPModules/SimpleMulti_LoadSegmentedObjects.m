function handles = SimpleMulti_LoadSegmentedObjects(handles)

% Help for the SimpleMulti_LoadSegmentedObjects module:
% Use SimpleMulti_LoadSegmentedObjects to load and align segmentations which were
% generated for another 4i cycle by a previous CP1 pipeline.
% Category: Other
%
%% *************************************************************************
%
%
% $Revision: 2159 $


%% [GG 20160809] Prolongued float format to allow exact background subtraction
format LONG
%% [GG 20160809] end, Note that the is a resetting of the format at the end of the module

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);


%textVAR01 = Which objects do you want to import (required)? (e.g.: Nuclei)
%defaultVAR01 = Nuclei
%infotypeVAR01 = objectgroup indep
ObjectName_Primary = char(handles.Settings.VariableValues{CurrentModuleNum,1});

%textVAR02 = Which second objects do you want to import (optional)? (e.g.: Cells). Ignore by keeping /
%defaultVAR02 = /
%infotypeVAR02 = objectgroup indep
ObjectName_Secondary = char(handles.Settings.VariableValues{CurrentModuleNum,2});

%pathnametextVAR03 = From which reference acquistion should objects be imported?
%defaultVAR03 = /BIOL/sonas/biol_uzh_pelkmans_s5/Data/Users/Thomas/151020-AribitraryPlate/
Pathname = char(handles.Settings.VariableValues{CurrentModuleNum,3});

%textVAR04 = Which images of the current pipeline should be used for aligning the acquisitions?
%infotypeVAR04 = imagegroup
OrigImageName = char(handles.Settings.VariableValues{CurrentModuleNum,4});
%inputtypeVAR04 = popupmenu

%textVAR05 = Which filter describes the images of the reference acquisition that should be used for aligning the acquisitions?
%defaultVAR05 = C01.
filterForImagesOfTrans = char(handles.Settings.VariableValues{CurrentModuleNum,5});

%textVAR06 = What is the background of the images?
%defaultVAR06 = 100
val_background = str2double(char(handles.Settings.VariableValues{CurrentModuleNum,6}));

%textVAR07 = What is the maximum shift you allow between cycles?
%defaultVAR07 = 200
max_shift_allowed = str2double(char(handles.Settings.VariableValues{CurrentModuleNum,7}));

%textVAR08 = How many times should we try good overlap?
%defaultVAR08 = 5
opt_rounds = str2double(char(handles.Settings.VariableValues{CurrentModuleNum,8}));

%textVAR09 = Do you want to correct your reference image for illumnation bias?
%defaultVAR09 = Yes
bool_illcorim = strcmpi(handles.Settings.VariableValues{CurrentModuleNum,9}, 'Yes');

%%%VariableRevisionNumber = 14



%%%%%%%%%%%%%%%%%%%%%%%%
%%% INITIALIZATION   %%%
%%%%%%%%%%%%%%%%%%%%%%%%

strTransPlate = Pathname;
if ~any(fileattrib(strTransPlate))
    error('Could not find reference plate');
end

if hasSecondaryObjectBeenDefined(ObjectName_Secondary) == false    % TS: ad hoc patch for including two objects: just do cheap computation twice rather than changing code
    ObjectName_Secondary = ObjectName_Primary;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Computation  %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%% FOR FIRST CYCLE %%%%%%%%%
if handles.Current.SetBeingAnalyzed == 1
    makeExternalBuffer(handles, strTransPlate, filterForImagesOfTrans, ObjectName_Primary, ObjectName_Secondary);
end

%%%%%% FOR ALL CYCLES %%%%%%%%%%
createShiftedSegmentation = true; % depending on results, this switch will possbily be turned (TS151021: creating a segmentation without objects).

[TiffFolder_trans, SegmentationFolder_trans] = getSubFoldersFromTransPlate(strTransPlate);
[strCorrespondingImage_trans, couldFindSameSite_image, strCorrespondingSegmentation_trans, couldFindSameSite_segmentation, strCorrespondingSegmentationSecondary_trans, couldFindSameSite_segmentation_Secondary] = ...
    getFileNameViaReferenceFile(handles, ObjectName_Primary, OrigImageName, ObjectName_Secondary);

% Get current image
Image_Cis = CPretrieveimage(handles,OrigImageName,ModuleName);

% Obtain coordinates of overalap of cis image (this acquisition) and trans image (reference acquistion)
if couldFindSameSite_image == true;
    %     readFun = @(x) double(imread(x)) ./ (2^16-1); %[GG20160419] old
    %     readFun = @(x) double(imread_illumination_corrected(x,1)) ./ (2^16-1); %[GG20160419] new
    
    if bool_illcorim %[GG20160621] added with user input 9
        readFun = @(x) double(imread_illumination_corrected(x,1)); %[GG20160419] new
    else
        readFun = @(x) double(imread(x)); %[GG20160419] new
    end
    
    Image_Trans = readFun(fullfile(TiffFolder_trans, strCorrespondingImage_trans));
      
    %% pad array image to get them to the same size, [GG20160621] added this section for the MPSimulation
    dim_cis = size(Image_Cis);
    dim_trans = size(Image_Trans);
    fprintf('Measure size of images\n')
    
    bool_diffsized = any(dim_cis ~= dim_trans);
    bool_transbigger = sum(dim_trans - dim_cis) >= 0;
    fprintf('Claculate dimensions of pad \n')
    
    if bool_diffsized
        if bool_transbigger
            dim_pad = ceil((dim_trans - dim_cis)./2);
            Image_Cis = padarray(Image_Cis, dim_pad, 0);
            Image_Cis = Image_Cis(1:dim_trans(1), 1:dim_trans(2));
            fprintf('Pad Image_Cis \n')
        else
            dim_pad = ceil((dim_cis - dim_trans)./2);
            Image_Trans = padarray(Image_Trans, dim_pad, 0);
            Image_Trans = Image_Trans(1:dim_cis(1), 1:dim_cis(2));
            fprintf('Pad Image_Trans \n')
        end
    end
%     
    %% get overlap
    fprintf('Getting overlap \n')
    [NSWE_Cis, NSWE_Trans] = getRobustAlignmentAtLowIntensity(Image_Trans, Image_Cis, max_shift_allowed, opt_rounds);
%     [NSWE_Cis, NSWE_Trans]  = getNSWEofOverlappingImageparts(Image_Cis, Image_Trans);
    
        
else
    createShiftedSegmentation = false;
end

% Obtain Segmentation from reference acquisition
if couldFindSameSite_segmentation == true;
    OrigSegmentation = double(imread(fullfile(SegmentationFolder_trans, strCorrespondingSegmentation_trans)));
    if bool_diffsized
        OrigSegmentation = padarray(OrigSegmentation, dim_pad, 0);
        if bool_transbigger
            OrigSegmentation = OrigSegmentation(1:dim_trans(1), 1:dim_trans(2));
        else
            OrigSegmentation = OrigSegmentation(1:dim_cis(1), 1:dim_cis(2));
        end
    end
    
else
    createShiftedSegmentation = false;
    OrigSegmentation = zeros(size(Image_Cis));
end

% Obtain Segmentation from secondaryObject
if couldFindSameSite_segmentation_Secondary == true;
    OrigSegmentation_Secondary = double(imread(fullfile(SegmentationFolder_trans, strCorrespondingSegmentationSecondary_trans)));
    if bool_diffsized
        OrigSegmentation_Secondary = padarray(OrigSegmentation_Secondary, dim_pad, 0);
        if bool_transbigger
            OrigSegmentation_Secondary = OrigSegmentation_Secondary(1:dim_trans(1), 1:dim_trans(2));
        else
            OrigSegmentation_Secondary = OrigSegmentation_Secondary(1:dim_cis(1), 1:dim_cis(2));
        end
    end
else
    OrigSegmentation_Secondary = zeros(size(Image_Cis));
end


if createShiftedSegmentation == true
    ShiftedSegmentation_Primary = zeros(size(Image_Cis),'double');
    ShiftedSegmentation_Secondary = zeros(size(Image_Cis),'double');
    
    % SHIFT ORIGINAL SEGMENTATION %
    ShiftedSegmentation_Primary(NSWE_Cis(1):NSWE_Cis(2), NSWE_Cis(3):NSWE_Cis(4)) = ...
        OrigSegmentation(NSWE_Trans(1):NSWE_Trans(2), NSWE_Trans(3):NSWE_Trans(4));
    
    ShiftedSegmentation_Secondary(NSWE_Cis(1):NSWE_Cis(2), NSWE_Cis(3):NSWE_Cis(4)) = ...
        OrigSegmentation_Secondary(NSWE_Trans(1):NSWE_Trans(2), NSWE_Trans(3):NSWE_Trans(4));
    
    % RELABEL  %
    
    % Since several CellProfiler modules will introduce error, if labelling
    % is not continous, make continuous labels;
    
    % SAVE RELABEL % %%[GG 20160208]
    % ensure save continuous relabeling %
    uOrigLabel_Primary  = unique(OrigSegmentation);
    uShiftedLabel_Primary = unique(ShiftedSegmentation_Primary);
    
    uOrigLabel_Secondary  = unique(OrigSegmentation_Secondary);
    uShiftedLabel_Secondary = unique(ShiftedSegmentation_Secondary);
    
    % find those nuclei/cells which were shfted out of image
    lostNucleiByShift = uOrigLabel_Primary(~ismember(uOrigLabel_Primary, uShiftedLabel_Primary));
    lostCellsByShift = uOrigLabel_Secondary(~ismember(uOrigLabel_Secondary, uShiftedLabel_Secondary));
    
    % combine both cells and nuclei shifted and correct for both
    lostObjectsByShift = unique([lostNucleiByShift; lostCellsByShift]);
    
    % reposition centroids of nuclei and cells to image border, shirink
    % them to single spot, this makes sure, that they are excluded once
    % loaded by GetRawProbwModelData2 as they are detected as bordercells
    [Rows, Columns] = size(OrigSegmentation);
    Repositioner = [Columns Rows];
    if any(lostObjectsByShift)% check: lost any nucleus
        t = regionprops(OrigSegmentation,'Centroid');
        t = struct2cell(t(lostObjectsByShift,:))';
        t = round(cell2mat(t));
        shift = NSWE_Cis([3 1])-NSWE_Trans([3 1]);
        lostNuclei_shiftedcentroids = t + repmat(shift,size(t,1),1);
        
        % check for nuclei out of range
        t = lostNuclei_shiftedcentroids(:,1) > Repositioner(1);
        lostNuclei_shiftedcentroids(t,1) = Repositioner(1);
        t = lostNuclei_shiftedcentroids(:,2) > Repositioner(2);
        lostNuclei_shiftedcentroids(t,2) = Repositioner(2);
        
        % check whether there are any objets with overlaping nuclear coordinates
        A = lostNuclei_shiftedcentroids(:,1);
        B = lostNuclei_shiftedcentroids(:,2);
        if size(unique(A),1) ~= size(A,1)
            [uRow_Pos] = unique(A);
            free_positions = 1:Repositioner(1);
            d = A;
            d(d < 1) = 1;
            d(d > Repositioner(1)) = Repositioner(1);
            free_positions(d) = [];
            % randomly assignes other Row position to cells with the same
            % Row position
            keep = cell(numel(uRow_Pos),1);
            for ix = 1:numel(uRow_Pos)
                tmp = find(A == uRow_Pos(ix));
                tmp = tmp(2:end);
                keep{ix} = tmp;
            end
            keep(cellfun(@(x) isempty(x), keep)) = [];
            vals2change = cat(1,keep{:});
            
            ix_R = randsample(size(free_positions,2), size(vals2change,1));
            lostNuclei_shiftedcentroids(vals2change,1) = free_positions(ix_R);
        end
        
        if size(unique(B),1) ~= size(B,1)
            [uCol_Pos] = unique(B);
            free_positions = 1:Repositioner(2);
            d = B;
            d(d < 1) = 1;
            d(d > Repositioner(2)) = Repositioner(2);
            free_positions(d) = [];
            % randomly assignes other Col position to cells with the same
            % col position
            keep = cell(numel(uCol_Pos),1);
            for ix = 1:numel(uCol_Pos)
                tmp = find(B == uCol_Pos(ix));
                tmp = tmp(2:end);
                keep{ix} = tmp;
            end
            keep(cellfun(@(x) isempty(x), keep)) = [];
            vals2change = cat(1,keep{:});
            
            ix_R = randsample(size(free_positions,2), size(vals2change,1));%[GG20161002] added replacement false
            lostNuclei_shiftedcentroids(vals2change,2) = free_positions(ix_R);
        end
        
        % put centroids to border od the image
        ix_X = lostNuclei_shiftedcentroids(:,1) <= 0;
        ix_Y = lostNuclei_shiftedcentroids(:,2) <= 0;
        lostNuclei_shiftedcentroids(ix_X,1) = 1;
        lostNuclei_shiftedcentroids(ix_Y,2) = 1;
        
        
        lin_lostnuclei = sub2ind(size(ShiftedSegmentation_Primary),lostNuclei_shiftedcentroids(:,2),lostNuclei_shiftedcentroids(:,1));
        ShiftedSegmentation_Primary(lin_lostnuclei) = lostObjectsByShift; % by this shifted nucleus will land within the current cell segmentation fo the object, but at the border of the image, so no relabeling needed
        lin_lostcells = sub2ind(size(ShiftedSegmentation_Secondary),lostNuclei_shiftedcentroids(:,2),lostNuclei_shiftedcentroids(:,1));
        ShiftedSegmentation_Secondary(lin_lostcells) = lostObjectsByShift; % by this shifted nucleus will land within the current cell segmentation fo the object, but at the border of the image, so no relabeling needed
    end
    
    
    %     if any(lostObjectsByShift)% check: lost any cells
    %         t = regionprops(OrigSegmentation,'Centroid');
    %         t = struct2cell(t(lostObjectsByShift,:))';
    %         t = cell2mat(t);
    %         shift = NSWE_Cis([3 1])-NSWE_Trans([3 1]);
    %         lostCells_shiftedcentroids = t + repmat(shift,size(t,1),1);
    %         lostCells_shiftedcentroids(lostCells_shiftedcentroids<0.5) = 1;
    %         lostCells_shiftedcentroids = round(lostCells_shiftedcentroids);
    %         Xpos_out_of_im = lostCells_shiftedcentroids(:,1) > size(OrigSegmentation,2);
    %         Ypos_out_of_im = lostCells_shiftedcentroids(:,2) > size(OrigSegmentation,1);
    %         if  any(Xpos_out_of_im)
    %             lostCells_shiftedcentroids(Xpos_out_of_im,1) = size(OrigSegmentation,2);
    %         end
    %         if  any(Ypos_out_of_im)
    %             lostCells_shiftedcentroids(Ypos_out_of_im,2) = size(OrigSegmentation,1);
    %         end
    %         lin_lostcells = sub2ind(size(ShiftedSegmentation_Secondary),lostCells_shiftedcentroids(:,2),lostCells_shiftedcentroids(:,1));
    %         ShiftedSegmentation_Secondary(lin_lostcells) = lostObjectsByShift; % by this shifted nucleus will land within the current cell segmentation fo the object, but at the border of the image, so no relabeling needed
    %     end
    
    % end Gabri code %
    
    uOrigLabel_Primary  = unique(OrigSegmentation);
    uShiftedLabel_Primary = unique(ShiftedSegmentation_Primary);
    
    uOrigLabel_Secondary  = unique(OrigSegmentation_Secondary);
    uShiftedLabel_Secondary = unique(ShiftedSegmentation_Secondary);
    
    
    if ~(any(uOrigLabel_Primary>0) && any(uShiftedLabel_Primary>0) && any(uOrigLabel_Secondary>0) && any(uShiftedLabel_Secondary>0))  % in case that at least one segmentation is empty
        createShiftedSegmentation = false;
    else
        uOrigLabel_Primary    = uOrigLabel_Primary(uOrigLabel_Primary~=0);   % ignore background
        uOrigLabel_Secondary  = uOrigLabel_Secondary(uOrigLabel_Secondary~=0);
        
        %[GG 20160208] Start Outcommented
        %         lFun = @(x) x(:);
        %         % check if labelling has been continous in original segmentation
        %         if (~isequal(lFun(uOrigLabel_Primary), lFun(unique(1:length(uOrigLabel_Primary))))) || (~isequal(lFun(uOrigLabel_Secondary), lFun(unique(1:length(uOrigLabel_Secondary)))))
        %             error('Labelling within original segmentation is not continuous. This will cause wrong data (but no compuational errors) in several standard modules of CellProfiler, and is therefore strongly advised not to do! Please reconsider your full analysis, instead of making your pipeline of ignore this discontinuous labelling!');
        %         else
        %             uShiftedLabel_Primary = uShiftedLabel_Primary(uShiftedLabel_Primary~=0); % ignore background
        %             uShiftedLabel_Secondary = uShiftedLabel_Secondary(uShiftedLabel_Secondary~=0);
        %
        %             commonUniqueElementsOfPrimaryAndSecondary = intersect(uShiftedLabel_Primary, uShiftedLabel_Secondary);
        %
        %
        %             elementsInShiftedLabel = length(commonUniqueElementsOfPrimaryAndSecondary);
        %
        %             uShiftedLabel_Shared = sort(commonUniqueElementsOfPrimaryAndSecondary,'ascend'); % do relabelling
        %             relabelledSegmentation_Primary = zeros(size(ShiftedSegmentation_Primary));
        %             for j=1:elementsInShiftedLabel
        %                 c = uShiftedLabel_Shared(j);
        %                 f = ShiftedSegmentation_Primary == c;
        %                 relabelledSegmentation_Primary(f) = j;
        %             end
        %
        %             relabelledSegmentation_Secondary = zeros(size(ShiftedSegmentation_Primary));
        %             for j=1:elementsInShiftedLabel
        %                 c = uShiftedLabel_Shared(j);
        %                 f = ShiftedSegmentation_Secondary == c;
        %                 relabelledSegmentation_Secondary(f) = j;
        %             end
        %
        %             correspondingOrigLabelForShifted = uShiftedLabel_Shared;
        %         end
        %[GG 20160208] END Outcommented
    end
    
end

% [GG 20160208] Start Added these two lines to make module work as before

relabelledSegmentation_Primary = ShiftedSegmentation_Primary;
relabelledSegmentation_Secondary = ShiftedSegmentation_Secondary;
commonUniqueElementsOfPrimaryAndSecondary = intersect(uShiftedLabel_Primary, uShiftedLabel_Secondary);
uShiftedLabel_Shared = sort(commonUniqueElementsOfPrimaryAndSecondary,'ascend'); % do relabelling
correspondingOrigLabelForShifted = uShiftedLabel_Shared;
% [GG 20160208] End

if createShiftedSegmentation == false  % create blank image (e.g.: if no object has been present, or some reference could not be found)
    relabelledSegmentation_Primary = zeros(size(Image_Cis),'double'); % default values
    relabelledSegmentation_Secondary = zeros(size(Image_Cis),'double');
    
    correspondingOrigLabelForShifted = 0;
    Centroid_Primary = [0 0];
    Centroid_Secondary = [0 0];
else % create further measurements, if valid objects are present
    tmp = regionprops(relabelledSegmentation_Primary,'Centroid');
    Centroid_Primary = cat(1,tmp.Centroid);
    
    tmp = regionprops(relabelledSegmentation_Secondary,'Centroid');
    Centroid_Secondary = cat(1,tmp.Centroid);
end

if numel(unique(relabelledSegmentation_Primary)) ~= numel(unique(OrigSegmentation))
    error('Obeject count between reference image and current cycle is wrong');
end
% %%%%%%%%%%%%%%
% %%% OUTPUT %%%
% %%%%%%%%%%%%%%

%% Saves the segmented image, not edited for objects along the edges or
%%% for size, to the handles structure.
fieldname = ['UneditedSegmented',ObjectName_Primary];
handles.Pipeline.(fieldname) = relabelledSegmentation_Primary;

%%% Saves the segmented image, only edited for small objects, to the
%%% handles structure.
fieldname = ['SmallRemovedSegmented',ObjectName_Primary];
handles.Pipeline.(fieldname) = relabelledSegmentation_Primary;

%%% Saves the final segmented label matrix image to the handles structure.
fieldname = ['Segmented',ObjectName_Primary];
handles.Pipeline.(fieldname) = relabelledSegmentation_Primary;

fieldname = ['UneditedSegmented',ObjectName_Secondary];
handles.Pipeline.(fieldname) = relabelledSegmentation_Secondary;

%%% Saves the segmented image, only edited for small objects, to the
%%% handles structure.
fieldname = ['SmallRemovedSegmented',ObjectName_Secondary];
handles.Pipeline.(fieldname) = relabelledSegmentation_Secondary;

%%% Saves the final segmented label matrix image to the handles structure.
fieldname = ['Segmented',ObjectName_Secondary];
handles.Pipeline.(fieldname) = relabelledSegmentation_Secondary;





%%% Saves the ObjectCount, i.e., the number of segmented objects.
%%% See comments for the Threshold saving above
if ~isfield(handles.Measurements.Image,'ObjectCountFeatures')
    handles.Measurements.Image.ObjectCountFeatures = {};
    handles.Measurements.Image.ObjectCount = {};
end

% primary
column = find(strcmp(handles.Measurements.Image.ObjectCountFeatures,ObjectName_Primary));
if isempty(column)
    handles.Measurements.Image.ObjectCountFeatures(end+1) = {ObjectName_Primary};
    column = length(handles.Measurements.Image.ObjectCountFeatures);
end
ObjCount = max(relabelledSegmentation_Primary(:));
handles.Measurements.Image.ObjectCount{handles.Current.SetBeingAnalyzed}(1,column) = ObjCount;

% secondary
column = find(strcmp(handles.Measurements.Image.ObjectCountFeatures,ObjectName_Secondary));
if isempty(column)
    handles.Measurements.Image.ObjectCountFeatures(end+1) = {ObjectName_Secondary};
    column = length(handles.Measurements.Image.ObjectCountFeatures);
end
ObjCount = max(relabelledSegmentation_Secondary(:));
handles.Measurements.Image.ObjectCount{handles.Current.SetBeingAnalyzed}(1,column) = ObjCount;


%%% Saves the location of each segmented object
% Follow CP convention for empty images (e.g.: as in IdentifySecondary
% module)
handles.Measurements.(ObjectName_Primary).LocationFeatures = {'CenterX','CenterY'};
handles.Measurements.(ObjectName_Primary).Location(handles.Current.SetBeingAnalyzed) = {Centroid_Primary};

%%% save relation to original object ID to handles
handles.Measurements.(ObjectName_Primary).UnshiftedObjectIdFeatures{handles.Current.SetBeingAnalyzed} = 'UnshiftedObjectId';
handles.Measurements.(ObjectName_Primary).UnshiftedObjectId{handles.Current.SetBeingAnalyzed} = correspondingOrigLabelForShifted;


% Follow CP convention for empty images (e.g.: as in IdentifySecondary
% module)
handles.Measurements.(ObjectName_Secondary).LocationFeatures = {'CenterX','CenterY'};
handles.Measurements.(ObjectName_Secondary).Location(handles.Current.SetBeingAnalyzed) = {Centroid_Secondary};

%%% save relation to original object ID to handles
handles.Measurements.(ObjectName_Secondary).UnshiftedObjectIdFeatures{handles.Current.SetBeingAnalyzed} = 'UnshiftedObjectId';
handles.Measurements.(ObjectName_Secondary).UnshiftedObjectId{handles.Current.SetBeingAnalyzed} = correspondingOrigLabelForShifted;


% [GG20160419] Replaces old segmentation with new freshly generated
% shifted one
% fieldname = ['Segmented',ObjectName_Primary];
% handles.Pipeline.(fieldname) = UneditedLabelMatrixImage;


%%%%%%%%%%%%%%%
%%% DISPLAY %%%
%%%%%%%%%%%%%%%


drawnow

if ~CPisHeadless()
    ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
    if any(findobj == ThisModuleFigureNumber)
        
        CPfigure(handles,'Image',ThisModuleFigureNumber);
        
        %%% A subplot of the figure window is set to display the original image.
        subplot(2,2,1);
        CPimagesc(Image_Cis,handles);
        title(['Input Image, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
        %%% A subplot of the figure window is set to display the colored label
        %%% matrix image.
        subplot(2,2,2);
        if couldFindSameSite_image == true
            im = Image_Trans;
        else
            im = zeros(size(Image_Cis));
        end
        CPimagesc(im,handles);
        title('Image from reference plate');
        
        subplot(2,2,3);
        ObjectOutlinesOnOrigImage_Primary = combineImageAndSegmentation(Image_Cis, relabelledSegmentation_Primary);
        ObjectOutlinesOnOrigImage_Secondary = combineImageAndSegmentation(Image_Cis, relabelledSegmentation_Secondary);
        ObjectOutlinesOnOrigImage = max(cat(3,ObjectOutlinesOnOrigImage_Primary, ObjectOutlinesOnOrigImage_Secondary),[],3);
        
        CPimagesc(ObjectOutlinesOnOrigImage,handles);
        title([ObjectName_Primary, ' and'  ObjectName_Secondary 'Outlines on Input Image']);
        
        
        subplot(2,2,4);
        CPimagesc(relabelledSegmentation_Secondary,handles);
        colormap('jet');
        clear ColoredLabelMatrixImage
        title(['Segmentation of ',ObjectName_Secondary]);
        
        drawnow
    end
end

%% [GG 20160809] Prolongued float format to allow exact background subtraction
format SHORT
%% [GG 20160809]end

end


function makeExternalBuffer(handles, strTransPlate, filterForImagesOfTrans, ObjectName_Primary, ObjectName_Secondary)
% note that preferentially all data would be stored in pipeline. However
% storing custom fields in handles.pipelines can lead to various errors.
% While storing data in handles.measurements is technically possible there
% would be various mistakes occuring if order of sites becomes rearranged
% in batch files. Since there will always only be one segmentation of a
% given name, there can be a buffer, that links to that name, and which is
% overwritten, if a new pipeline is made


[TiffFolder_trans, SegmentationFolder_trans] = getSubFoldersFromTransPlate(strTransPlate);



% Images
eB.strImages_trans = getFilesAndDirectories(TiffFolder_trans, filterForImagesOfTrans);
[eB.Row_image_trans, eB.intColumn_image_trans, eB.intImagePosition_image_trans] = cellfun(@(x) MetaFromImageName(x), eB.strImages_trans, 'UniformOutput', true);

% Segmentations (primary)
filterForSegmentationsOfTrans = ['Segmented' ObjectName_Primary '\.'];
eB.strSegmentations_trans = getFilesAndDirectories(SegmentationFolder_trans, filterForSegmentationsOfTrans);
if isempty(eB.strSegmentations_trans)
    error(['Could not find Segementations for ' ObjectName_Primary]);
end

[eB.Row_segmentations_trans, eB.intColumn_segmentations_trans, eB.intImagePosition_segmentations_trans] = cellfun(@(x) MetaFromImageName(x), eB.strSegmentations_trans, 'UniformOutput', true);

% Segmentations (secondary)
if hasSecondaryObjectBeenDefined(ObjectName_Secondary) == true
    filterForSegmentationsOfTrans_Secondary = ['Segmented' ObjectName_Secondary '\.'];
    eB.strSegmentationsSecondary_trans = getFilesAndDirectories(SegmentationFolder_trans, filterForSegmentationsOfTrans_Secondary);
    [eB.Row_segmentationsSecondary_trans, eB.intColumn_segmentationsSecondary_trans, eB.intImagePosition_segmentationsSecondary_trans] = ...
        cellfun(@(x) MetaFromImageName(x), eB.strSegmentationsSecondary_trans, 'UniformOutput', true);
    
    if isempty(eB.strSegmentationsSecondary_trans)
        error(['Could not find Segementations for ' ObjectName_Secondary]);
    end
end



strBufferFile = getFileNameOfBuffer(handles, ObjectName_Primary, ObjectName_Secondary);

if any(fileattrib(strBufferFile))
    [~, ex] = fileparts(strBufferFile);
    fprintf([ex ' already exists. It will be overwritten with current one.\n']);
end

save(strBufferFile, 'eB');
end

function strBufferFile = getFileNameOfBuffer(handles, ObjectName_Primary, ObjectName_Secondary)

outDir = handles.Current.DefaultOutputDirectory;
ex = ['SimpleMulti_' ObjectName_Primary ObjectName_Secondary '.mat'];
strBufferFile = fullfile(outDir, ex );

end

function SecondaryIsDefined = hasSecondaryObjectBeenDefined(ObjectName_Secondary)
SecondaryIsDefined = ~isequal(ObjectName_Secondary,'/'); % no secondary specified;

end


function [strCorrespondingImage_trans, couldFindSameSite_image, strCorrespondingSegmentation_trans, couldFindSameSite_segmentation, strCorrespondingSegmentationSecondary_trans, couldFindSameSite_segmentation_Secondary] = ...
    getFileNameViaReferenceFile(handles, ObjectName_Primary, OrigImageName, ObjectName_Secondary)
% Avoid repeated, independent access to buffer file by reading segmentation
% and images after same loading event. Note that the code thus becomes more
% ugly and less readable.


% SetBeingAnalyzed = handles.Current.SetBeingAnalyzed;
% size(handles.Pipeline.(['Filename' OrigImageName])) %[GG20160419] old
% SetBeingAnalyzed %[GG20160419] old
% ImageName_cis = handles.Pipeline.(['Filename' OrigImageName]){SetBeingAnalyzed}; %[GG20160419] old
ImageName_cis = handles.Pipeline.(['Filename' OrigImageName]){1}; %[GG20160419] new

strBufferFile = getFileNameOfBuffer(handles, ObjectName_Primary, ObjectName_Secondary);

if ~any(fileattrib(strBufferFile))
    error('Could not find buffer file for object, that should be created during first cycle');
else
    load(strBufferFile);
end

[Row_cis, intColumn_cis, intImagePosition_cis] = MetaFromImageName(ImageName_cis);

% get corresponding image from trans
correspondsToSameSite_image = eB.Row_image_trans == Row_cis & eB.intColumn_image_trans == intColumn_cis & eB.intImagePosition_image_trans == intImagePosition_cis;
if sum(correspondsToSameSite_image) == 0
    couldFindSameSite_image = false;
    strCorrespondingImage_trans = '';
elseif sum(correspondsToSameSite_image) == 1;
    couldFindSameSite_image = true;
    strCorrespondingImage_trans = eB.strImages_trans{correspondsToSameSite_image};
else
    error('Could not unambiguously find corresponding image of other dataset. Please set more stringent filters.');
end

% get corresponding segmentation from trans (primary)
correspondsToSameSite_segmentation = eB.Row_segmentations_trans == Row_cis & eB.intColumn_segmentations_trans == intColumn_cis & eB.intImagePosition_segmentations_trans == intImagePosition_cis;
if sum(correspondsToSameSite_segmentation) == 0
    couldFindSameSite_segmentation = false;
    strCorrespondingSegmentation_trans = '';
elseif sum(correspondsToSameSite_segmentation) == 1;
    couldFindSameSite_segmentation = true;
    strCorrespondingSegmentation_trans = eB.strSegmentations_trans{correspondsToSameSite_segmentation};
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

%% subfunctions required for gerRobustAlignmentAtlowIntensity
function [NSWE_Cis, NSWE_Trans] = getRobustAlignmentAtLowIntensity(Image_Trans, Image_Cis, max_shift_allowed, opt_rounds)
%% getRobustAlignmentAtLowIntensity(Image_Trans, Image_Cis, max_shift_allowed, opt_rounds)
%% Image_Cis = reference image, Image_Trans = image to be shifted, max_shift_allowed = maximal shift to be expected between images, opt_rounds = maximal rouds of shift seeking

%% rescale image
Image_Trans = imrescale(Image_Trans, 0.001, 0.99);
Image_Cis = imrescale(Image_Cis, 0.001, 0.99);

%% find shift
[NSWE_Cis, NSWE_Trans]  = getNSWEofOverlappingImageparts(Image_Cis, Image_Trans);

%% while loop to find good shift
current_round = 0;
shifts = abs(NSWE_Cis - NSWE_Trans);
while any(shifts > max_shift_allowed) && current_round <= opt_rounds
    
    % rescale
    Image_Trans = imrescale(Image_Trans, 0.001, 0.99);
    Image_Cis = imrescale(Image_Cis, 0.001, 0.99);
    
    % find overlap
    [NSWE_Cis, NSWE_Trans]  = getNSWEofOverlappingImageparts(Image_Cis, Image_Trans);
    
    % calculate absolute shifts to check for while loop
    shifts = abs(NSWE_Cis - NSWE_Trans);
    
    % update counter
    current_round = current_round + 1;
    
    sprintf('%d. round of rescaling',  current_round+1')
end
sprintf('%d x Rescaled images used for registration\n', current_round+1)
end

%%% subfunctions
function [im_res] = imrescale(im, low_quantile, top_quantile)

val_res = quantile(im(:),[low_quantile top_quantile]);

im_res = (im - val_res(1))./(val_res(2) - val_res(1));
im_res(im_res > 1) = 1;
im_res(im_res < 0) = 0;

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
function [matMeanImage, matStdImage, hasIlluminationCorrection] = getIlluminationReference(strBatchDir,iChannel,cacheInRam)
if nargin < 3
    cacheInRam = false;
end
fprintf('Within SimpleMulti function used: getIlluminationReference \n')

strPathToCurrentIllumination = fullfile(strBatchDir,...
    sprintf('Measurements_batch_illcor_channel%03d_zstack000.mat',iChannel));

matMeanImage = [];
matStdImage =[];

if ~any(fileattrib(strPathToCurrentIllumination))
    hasIlluminationCorrection = false;
    warning('matlab:bsBla','%s:  failed to load illumination correction %s',mfilename,strPathToCurrentIllumination);
else
    hasIlluminationCorrection = true;
    if cacheInRam == false
        ImportedIlluminationCorrection = load(strPathToCurrentIllumination);
        matMeanImage = double(ImportedIlluminationCorrection.stat_values.mean);
        matStdImage = double(ImportedIlluminationCorrection.stat_values.std);
    else
        ImportedIlluminationCorrection = cacheForIlluminationStats(strPathToCurrentIllumination);
        matMeanImage = ImportedIlluminationCorrection.stat_values.mean;  % conversion to double in cache to save time
        matStdImage = ImportedIlluminationCorrection.stat_values.std;
    end
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% COORDINATES OF OVERLAPPING IMAGE   %%%%%%%%%%%%%%%%%%%

function [NSWE_Cis, NSWE_Trans] = getNSWEofOverlappingImageparts(Image_Cis, Image_Trans)
% upon providing two images (Image_Cis, Image_Trans), The top row (N -orth),
% lowest row (S -outh),  left column (W -est) and right column (E -east) of
% the overlapping region are returned;

[imageHeight, imageWidth] = size(Image_Cis);

fourierTrans = fft2(Image_Trans);
fourierCis = fft2(Image_Cis);
[registrationOutput] = dftregistration(fourierTrans,fourierCis,1); %[GG20160420] orig is 1 I tried 10

SurplusRows = registrationOutput(3);
SurplusColumns = registrationOutput(4);

% Coordinates of for Rows

if SurplusRows > 0
    
    N_cis = 1;
    S_cis = imageHeight - SurplusRows;
    
    N_trans = 1 + SurplusRows;
    S_trans = imageHeight;
    
elseif SurplusRows == 0
    
    N_cis = 1;
    S_cis = imageHeight;
    N_trans = 1;
    S_trans = imageHeight;
    
elseif SurplusRows < 0
    
    N_cis = -SurplusRows + 1;
    S_cis = imageHeight;
    
    N_trans = 1;
    S_trans = imageHeight + SurplusRows;
    
else
    error('something terribly wrong');
end

% Coordinates of for Columns

if SurplusColumns > 0
    
    W_cis = 1;
    E_cis = imageWidth - SurplusColumns;
    
    W_trans = 1 + SurplusColumns;
    E_trans = imageWidth;
    
elseif SurplusColumns == 0
    
    W_cis = 1;
    E_cis = imageWidth;
    
    W_trans = 1;
    E_trans = imageWidth;
    
elseif SurplusColumns < 0
    
    W_cis = -SurplusColumns + 1;
    E_cis = imageWidth;
    
    W_trans = 1;
    E_trans = imageWidth + SurplusColumns;
    
else
    error('something terribly wrong');
end

NSWE_Cis = [N_cis, S_cis, W_cis, E_cis];
NSWE_Trans = [N_trans, S_trans, W_trans, E_trans];


end
function ObjectOutlinesOnOrigImage = combineImageAndSegmentation(OrigImage, FinalLabelMatrixImage)
% This code is copy/pasted from original IdentifySecondary


MaxFilteredImage = ordfilt2(FinalLabelMatrixImage,9,ones(3,3),'symmetric');
%%% Determines the outlines.
IntensityOutlines = FinalLabelMatrixImage - MaxFilteredImage;
%%% [PLab] hack.s ave memory.
clear MaxFilteredImage;
%%% Converts to logical.
warning off MATLAB:conversionToLogical
LogicalOutlines = logical(IntensityOutlines);
%%% [PLab] hack.s ave memory.
clear IntensityOutlines;
warning on MATLAB:conversionToLogical

ObjectOutlinesOnOrigImage = OrigImage;

qmin = quantile(ObjectOutlinesOnOrigImage(:),0.01); % [TS] Keep intention of original code by rescaling without outliers
qmax = quantile(ObjectOutlinesOnOrigImage(:),0.99);
ObjectOutlinesOnOrigImage(ObjectOutlinesOnOrigImage < qmin) = qmin;
ObjectOutlinesOnOrigImage(ObjectOutlinesOnOrigImage > qmax) = qmax;
ObjectOutlinesOnOrigImage = (ObjectOutlinesOnOrigImage-qmin) ./ (qmax-qmin);

ObjectOutlinesOnOrigImage = insertRedLine(ObjectOutlinesOnOrigImage, LogicalOutlines);

end
function RGBimage = insertRedLine(image, bwLineImage)

if size(image,3) == 1
    image = repmat(image,[1 1 3]);
end

r = image(:,:,1);
g = image(:,:,2);
b = image(:,:,3);

r(bwLineImage) = 1;
g(bwLineImage) = 0;
b(bwLineImage) = 0;

RGBimage = cat(3, r, g, b);
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



%%%%%% Illumionation correction %%%%%%%
function corr_image = imread_illumination_corrected(strImage, shallCacheCorrection)
% reads images and applies illumination correction
fprintf('Within SimpleMulti function used: imread_illumination_corrected \n')


% strImage = 'Z:\Data\Users\RNAFish\IndividualExperiments\150430_mycHprtRescanB\JoinedHprtMycRescan\TIFF\Q3_A01_T0001F109L01A03Z01C03.png';
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
function dat = cacheForIlluminationStats(strFileName)
% Initialize Persistent variables for caching
fprintf('Within SimpleMulti function used: cacheForIlluminationStats \n')
persistent CachedMeasurments;
persistent OriginalPathOfChached;

if isempty(CachedMeasurments)
    CachedMeasurments = cell(0);
end

if isempty(OriginalPathOfChached)
    OriginalPathOfChached = cell(0);
end

nStrFileName = strFileName; % ensure that each file only stored once in cache once

[isCached cachedLoc]= ismember(nStrFileName,OriginalPathOfChached);

if ~isCached   % load into cache, if absent there
    fprintf('Caching illumination correction ... ');
    cachedLoc = length(CachedMeasurments) + 1;
    
    ImportedIlluminationCorrection = load(nStrFileName);
    ImportedIlluminationCorrection.stat_values.matMeanImage = double(ImportedIlluminationCorrection.stat_values.mean);
    ImportedIlluminationCorrection.stat_values.std = double(ImportedIlluminationCorrection.stat_values.std);
    
    OriginalPathOfChached{cachedLoc} = nStrFileName;
    CachedMeasurments{cachedLoc} = ImportedIlluminationCorrection;
    
    fprintf('complete \n');
end

dat = CachedMeasurments{cachedLoc}; % retreive data
end
function CorrectedImage = IllumCorrect(Image,matMean,matStd,isLog)
%[NB] this functions does the Illumination correction by Z-scorring and
%reconstructing the original image values
%Usage:
%CorrImage = ILLUMCORRECT(IMAGE,MATMEAN,MATSTD,ISLOG). Where IMAGE is the
%original image to be corrected, MATMEAN is the perpixel mean values of the
%image, MATSTD is the perpixel standard deviation of the image and ISLOG is
%1 when MATMEAN and MATSTD are calculated form log10 transformed images and
%0 when they are same scale as IMAGE. If left undefined the default value
%is 1.
fprintf('Within SimpleMulti function used: IllumCorrect \n')
% convert input into double, which is required for calculation
Image = double(Image);
matMean = double(matMean);
matStd = double(matStd);


%Check inputs are at least 3
if nargin < 3
    error('%s: The minimum number of inputs is 3, Please check you have the correct number of inputs.',mfilename)
end

%check the fourth input
if nargin < 4
    warning('%s: No isLog value imputed. Asuming that mean and std values are in log10 scale.',mfilename)
    isLog = 1;
end


%Check that the mean std deviation images are of the same size
if ~(sum(size(matMean)==size(matStd)) == 2)
    error('%s: the Mean and Std matrices must have the same size.',mfilename)
end


%calculate the resize factor for matMean and matStd
ReFact = size(Image)./size(matMean);
if ReFact(1)~=ReFact(2)
    error('%s: the size of the input Image is not a multiple or equal to the size of the Mean and Std matrices. Please check the inputs.',mfilename)
end
ReFact = ReFact(1);

%Resize matMean and matStd
matMean = imresize(matMean,ReFact);
matStd = imresize(matStd,ReFact);


% do correction
if isLog == 1
    % Avoid -Inf values after log10 transform.
    Image(Image == 0) = 1;
    % Apply z-score normalization for each single pixel.
    CorrectedImage = (log10(Image)-matMean)./matStd;
    % Reverse z-score.
    CorrectedImage = (CorrectedImage.*mean(matStd(:)))+mean(matMean(:));
    % Reverse log10 transform that was applied to images when learning
    % mean/std statistics as well the corrected image.
    CorrectedImage = 10.^CorrectedImage;
else
    % Apply z-score normalization for each single pixel.
    CorrectedImage = (Image-matMean)./matStd;
    % Reverse z-score.
    CorrectedImage = (CorrectedImage.*mean(matStd(:)))+mean(matMean(:));
end

% Avoid negative values in the final image
CorrectedImage(CorrectedImage<0)=0;

% Fix potentially broken pixels (which are not variable)
CorrectedImage = fixNonNumericalValueInImage(CorrectedImage);

end
function matImage = fixNonNumericalValueInImage(matImage)
% replaces NaN and Inf within matIMAGE by an estimate of the local
% intensities (derived from smoothing the surrounding);
fprintf('Within SimpleMulti function used: fixNonNumericalValueInImage \n')

% Predefined Constants
minimalBlur = 30;
maxSigmaBlur = 5;

% Find bad pixels
bwBadPixels = isinf(matImage) | isnan(matImage);

if ~any(bwBadPixels(:)) % no bad pixel
    return
elseif ~any(~bwBadPixels) % all bad pixels
    fprintf('%s: no pixel has a numerical value \n',mfilename)
    return
else % some bad pixels
    
    % Estimate size of artifacts
    CoordinatesOfBounding = cell2mat(struct2cell(regionprops(bwBadPixels,'BoundingBox')));
    MaxmialObjectDiameterOfArtifact = ceil(max([max(CoordinatesOfBounding(:,2)) max(CoordinatesOfBounding(:,4))]));
    SmoothingSize = max([minimalBlur 2*MaxmialObjectDiameterOfArtifact]);
    numRows = size(matImage,1);
    numColumns = size(matImage,2);
    SmoothingSigma = min([maxSigmaBlur round(SmoothingSize./2)]);
    
    % Expand bad pixels (in a quick way by boxing): only process these
    % regions in later steps
    ExpandedObjects = false(size(matImage));
    
    for j=1:size(CoordinatesOfBounding,1)
        N = floor(CoordinatesOfBounding(j,2) - SmoothingSize);
        S = ceil(CoordinatesOfBounding(j,2)+CoordinatesOfBounding(j,4) + SmoothingSize);
        W = floor(CoordinatesOfBounding(j,1) - SmoothingSize);
        E = ceil(CoordinatesOfBounding(j,1)+CoordinatesOfBounding(j,3) + SmoothingSize);
        
        N = max([1 N]);
        S = min([numRows S]);
        W = max([1 W]);
        E = min([numColumns E]);
        
        ExpandedObjects(N:S,W:E) = true;
    end
    
    % Smoothen image (only within boxed regions)
    CoordinatesOfBounding = cell2mat(struct2cell(regionprops(ExpandedObjects,'BoundingBox')));
    SmoothenedImage = zeros(size(matImage));
    
    for j=1:size(CoordinatesOfBounding,1)
        N = floor(CoordinatesOfBounding(j,2));
        S = ceil(CoordinatesOfBounding(j,2)+CoordinatesOfBounding(j,4));
        W = floor(CoordinatesOfBounding(j,1));
        E = ceil(CoordinatesOfBounding(j,1)+CoordinatesOfBounding(j,3));
        
        N = max([1 N]);
        S = min([numRows S]);
        W = max([1 W]);
        E = min([numColumns E]);
        
        CurrCropImage = matImage(N:S,W:E);
        CurrBwBadPixels = bwBadPixels(N:S,W:E);
        
        % 1st round: Replace by Local median
        LocalIntensities = CurrCropImage(:);
        hasNoNumericalValue = isinf(LocalIntensities) | isnan(LocalIntensities);
        LocalIntensities = LocalIntensities(~hasNoNumericalValue);
        if any(LocalIntensities)
            CurrCropImage(CurrBwBadPixels) = median(LocalIntensities);
        else
            CurrCropImage(CurrBwBadPixels) = 0;
        end
        
        % 2nd round: Smooth locally
        H = fspecial('gaussian',[SmoothingSize SmoothingSize],SmoothingSize./SmoothingSigma);
        SmoothenedImage(N:S,W:E) = imfilter(CurrCropImage,H,'symmetric');
    end
end

matImage(bwBadPixels) = SmoothenedImage(bwBadPixels);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% CODE FROM MATLABCENTRAL   %%%%%%%%%%%%%%%%%%%

function [output, Greg] = dftregistration(buf1ft,buf2ft,usfac)
% function [output Greg] = dftregistration(buf1ft,buf2ft,usfac);
% Efficient subpixel image registration by crosscorrelation. This code
% gives the same precision as the FFT upsampled cross correlation in a
% small fraction of the computation time and with reduced memory
% requirements. It obtains an initial estimate of the crosscorrelation peak
% by an FFT and then refines the shift estimation by upsampling the DFT
% only in a small neighborhood of that estimate by means of a
% matrix-multiply DFT. With this procedure all the image points are used to
% compute the upsampled crosscorrelation.
% Manuel Guizar - Dec 13, 2007

% Portions of this code were taken from code written by Ann M. Kowalczyk
% and James R. Fienup.
% J.R. Fienup and A.M. Kowalczyk, "Phase retrieval for a complex-valued
% object by using a low-resolution image," J. Opt. Soc. Am. A 7, 450-458
% (1990).

% Citation for this algorithm:
% Manuel Guizar-Sicairos, Samuel T. Thurman, and James R. Fienup,
% "Efficient subpixel image registration algorithms," Opt. Lett. 33,
% 156-158 (2008).

% Inputs
% buf1ft    Fourier transform of reference image,
%           DC in (1,1)   [DO NOT FFTSHIFT]
% buf2ft    Fourier transform of image to register,
%           DC in (1,1) [DO NOT FFTSHIFT]
% usfac     Upsampling factor (integer). Images will be registered to
%           within 1/usfac of a pixel. For example usfac = 20 means the
%           images will be registered within 1/20 of a pixel. (default = 1)

% Outputs
% output =  [error,diffphase,net_row_shift,net_col_shift]
% error     Translation invariant normalized RMS error between f and g
% diffphase     Global phase difference between the two images (should be
%               zero if images are non-negative).
% net_row_shift net_col_shift   Pixel shifts between images
% Greg      (Optional) Fourier transform of registered version of buf2ft,
%           the global phase difference is compensated for.

% Default usfac to 1
if exist('usfac')~=1, usfac=1; end

% Compute error for no pixel shift
if usfac == 0,
    CCmax = sum(sum(buf1ft.*conj(buf2ft)));
    rfzero = sum(abs(buf1ft(:)).^2);
    rgzero = sum(abs(buf2ft(:)).^2);
    error = 1.0 - CCmax.*conj(CCmax)/(rgzero*rfzero);
    error = sqrt(abs(error));
    diffphase=atan2(imag(CCmax),real(CCmax));
    output=[error,diffphase];
    
    % Whole-pixel shift - Compute crosscorrelation by an IFFT and locate the
    % peak
elseif usfac == 1,
    [m,n]=size(buf1ft);
    CC = ifft2(buf1ft.*conj(buf2ft));
    [max1,loc1] = max(CC);
    [max2,loc2] = max(max1);
    rloc=loc1(loc2);
    cloc=loc2;
    CCmax=CC(rloc,cloc);
    rfzero = sum(abs(buf1ft(:)).^2)/(m*n);
    rgzero = sum(abs(buf2ft(:)).^2)/(m*n);
    error = 1.0 - CCmax.*conj(CCmax)/(rgzero(1,1)*rfzero(1,1));
    error = sqrt(abs(error));
    diffphase=atan2(imag(CCmax),real(CCmax));
    md2 = fix(m/2);
    nd2 = fix(n/2);
    if rloc > md2
        row_shift = rloc - m - 1;
    else
        row_shift = rloc - 1;
    end
    
    if cloc > nd2
        col_shift = cloc - n - 1;
    else
        col_shift = cloc - 1;
    end
    output=[error,diffphase,row_shift,col_shift];
    
    % Partial-pixel shift
else
    
    % First upsample by a factor of 2 to obtain initial estimate
    % Embed Fourier data in a 2x larger array
    [m,n]=size(buf1ft);
    mlarge=m*2;
    nlarge=n*2;
    CC=zeros(mlarge,nlarge);
    CC(m+1-fix(m/2):m+1+fix((m-1)/2),n+1-fix(n/2):n+1+fix((n-1)/2)) = ...
        fftshift(buf1ft).*conj(fftshift(buf2ft));
    
    % Compute crosscorrelation and locate the peak
    CC = ifft2(ifftshift(CC)); % Calculate cross-correlation
    [max1,loc1] = max(CC);
    [max2,loc2] = max(max1);
    rloc=loc1(loc2);cloc=loc2;
    CCmax=CC(rloc,cloc);
    
    % Obtain shift in original pixel grid from the position of the
    % crosscorrelation peak
    [m,n] = size(CC); md2 = fix(m/2); nd2 = fix(n/2);
    if rloc > md2
        row_shift = rloc - m - 1;
    else
        row_shift = rloc - 1;
    end
    if cloc > nd2
        col_shift = cloc - n - 1;
    else
        col_shift = cloc - 1;
    end
    row_shift=row_shift/2;
    col_shift=col_shift/2;
    
    % If upsampling > 2, then refine estimate with matrix multiply DFT
    if usfac > 2,
        %%% DFT computation %%%
        % Initial shift estimate in upsampled grid
        row_shift = round(row_shift*usfac)/usfac;
        col_shift = round(col_shift*usfac)/usfac;
        dftshift = fix(ceil(usfac*1.5)/2); %% Center of output array at dftshift+1
        % Matrix multiply DFT around the current shift estimate
        CC = conj(dftups(buf2ft.*conj(buf1ft),ceil(usfac*1.5),ceil(usfac*1.5),usfac,...
            dftshift-row_shift*usfac,dftshift-col_shift*usfac))/(md2*nd2*usfac^2);
        % Locate maximum and map back to original pixel grid
        [max1,loc1] = max(CC);
        [max2,loc2] = max(max1);
        rloc = loc1(loc2); cloc = loc2;
        CCmax = CC(rloc,cloc);
        rg00 = dftups(buf1ft.*conj(buf1ft),1,1,usfac)/(md2*nd2*usfac^2);
        rf00 = dftups(buf2ft.*conj(buf2ft),1,1,usfac)/(md2*nd2*usfac^2);
        rloc = rloc - dftshift - 1;
        cloc = cloc - dftshift - 1;
        row_shift = row_shift + rloc/usfac;
        col_shift = col_shift + cloc/usfac;
        
        % If upsampling = 2, no additional pixel shift refinement
    else
        rg00 = sum(sum( buf1ft.*conj(buf1ft) ))/m/n;
        rf00 = sum(sum( buf2ft.*conj(buf2ft) ))/m/n;
    end
    error = 1.0 - CCmax.*conj(CCmax)/(rg00*rf00);
    error = sqrt(abs(error));
    diffphase=atan2(imag(CCmax),real(CCmax));
    % If its only one row or column the shift along that dimension has no
    % effect. We set to zero.
    if md2 == 1,
        row_shift = 0;
    end
    if nd2 == 1,
        col_shift = 0;
    end
    output=[error,diffphase,row_shift,col_shift];
end

% Compute registered version of buf2ft
if (nargout > 1)&&(usfac > 0),
    [nr,nc]=size(buf2ft);
    Nr = ifftshift([-fix(nr/2):ceil(nr/2)-1]);
    Nc = ifftshift([-fix(nc/2):ceil(nc/2)-1]);
    [Nc,Nr] = meshgrid(Nc,Nr);
    Greg = buf2ft.*exp(i*2*pi*(-row_shift*Nr/nr-col_shift*Nc/nc));
    Greg = Greg*exp(i*diffphase);
elseif (nargout > 1)&&(usfac == 0)
    Greg = buf2ft*exp(i*diffphase);
end
return
end
function out=dftups(in,nor,noc,usfac,roff,coff)
% function out=dftups(in,nor,noc,usfac,roff,coff);
% Upsampled DFT by matrix multiplies, can compute an upsampled DFT in just
% a small region.
% usfac         Upsampling factor (default usfac = 1)
% [nor,noc]     Number of pixels in the output upsampled DFT, in
%               units of upsampled pixels (default = size(in))
% roff, coff    Row and column offsets, allow to shift the output array to
%               a region of interest on the DFT (default = 0)
% Recieves DC in upper left corner, image center must be in (1,1)
% Manuel Guizar - Dec 13, 2007
% Modified from dftus, by J.R. Fienup 7/31/06

% This code is intended to provide the same result as if the following
% operations were performed
%   - Embed the array "in" in an array that is usfac times larger in each
%     dimension. ifftshift to bring the center of the image to (1,1).
%   - Take the FFT of the larger array
%   - Extract an [nor, noc] region of the result. Starting with the
%     [roff+1 coff+1] element.

% It achieves this result by computing the DFT in the output array without
% the need to zeropad. Much faster and memory efficient than the
% zero-padded FFT approach if [nor noc] are much smaller than [nr*usfac nc*usfac]

[nr,nc]=size(in);
% Set defaults
if exist('roff')~=1, roff=0; end %#ok<*EXIST>
if exist('coff')~=1, coff=0; end
if exist('usfac')~=1, usfac=1; end
if exist('noc')~=1, noc=nc; end
if exist('nor')~=1, nor=nr; end
% Compute kernels and obtain DFT by matrix products
kernc=exp((-i*2*pi/(nc*usfac))*( ifftshift([0:nc-1]).' - floor(nc/2) )*( [0:noc-1] - coff ));
kernr=exp((-i*2*pi/(nr*usfac))*( [0:nor-1].' - roff )*( ifftshift([0:nr-1]) - floor(nr/2)  ));
out=kernr*in*kernc;
return

end


