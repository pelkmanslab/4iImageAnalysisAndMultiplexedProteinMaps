function handles = LoadImagesTrans(handles)

% Help for the Load Images Trans module:
% LoadImagesTrans is heavily based on LoadMoreImages
% Category: File Processing

% $Revision: 1725 $

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = How do you want to load these files?
%choiceVAR01 = Text-Exact match
%choiceVAR01 = Text-Regular expressions
%choiceVAR01 = Order
LoadChoice = char(handles.Settings.VariableValues{CurrentModuleNum,1});
%inputtypeVAR01 = popupmenu

if strcmp(LoadChoice(1),'T')
    ExactOrRegExp = LoadChoice(6);
end

%textVAR02 = Type the text that one type of image has in common (for TEXT options), or their position in each group (for ORDER option):
%defaultVAR02 = C01.
TextToFind{1} = char(handles.Settings.VariableValues{CurrentModuleNum,2});

%textVAR03 = What do you want to call these images within CellProfiler?
%defaultVAR03 = TransOrigBlue
%infotypeVAR03 = imagegroup indep
ImageTransName{1} = char(handles.Settings.VariableValues{CurrentModuleNum,3});

%textVAR04 = To which images does the image in trans correspond?
%defaultVAR04 = OrigBlue
%infotypeVAR04 = imagegroup indep
ImageCisName{1} = char(handles.Settings.VariableValues{CurrentModuleNum,4});

%textVAR05 = Type the text that one type of image has in common (for TEXT options), or their position in each group (for ORDER option):
%defaultVAR05 = C02.
TextToFind{2} = char(handles.Settings.VariableValues{CurrentModuleNum,5});

%textVAR06 = What do you want to call these images within CellProfiler?
%defaultVAR06 = TransOrigGreen
%infotypeVAR06 = imagegroup indep
ImageTransName{2} = char(handles.Settings.VariableValues{CurrentModuleNum,6});

%textVAR07 = To which images does the image in trans correspond?
%defaultVAR07 = OrigGreen
%infotypeVAR07 = imagegroup indep
ImageCisName{2} = char(handles.Settings.VariableValues{CurrentModuleNum,7});

%textVAR08 = Type the text that one type of image has in common (for TEXT options), or their position in each group (for ORDER option):
%defaultVAR08 = C03.
TextToFind{3} = char(handles.Settings.VariableValues{CurrentModuleNum,8});

%textVAR09 = What do you want to call these images within CellProfiler?
%defaultVAR09 = TransOrigRed
%infotypeVAR09 = imagegroup indep
ImageTransName{3} = char(handles.Settings.VariableValues{CurrentModuleNum,9});

%textVAR10 = To which images does the image in trans correspond?
%defaultVAR10 = OrigRed
%infotypeVAR10 = imagegroup indep
ImageCisName{3} = char(handles.Settings.VariableValues{CurrentModuleNum,10});

%textVAR11 = Do you also want to correct images for illumination bias ?
%choiceVAR11 = No
%choiceVAR11 = Yes
DoIllCorrIm = char(handles.Settings.VariableValues{CurrentModuleNum,11});
%inputtypeVAR11 = popupmenu

%pathnametextVAR12 = From which reference acquistion should objects be imported?
%defaultVAR12 = L:\Data\Users\Gabriele\Multiplexing\20160529_MPSimulation01\20160531_MPSimulation01_stain02
Pathname = char(handles.Settings.VariableValues{CurrentModuleNum,12});


%%%VariableRevisionNumber = 1

%%%%%%%%%%%%%%%%%%%%%%%%
%%% INITIALIZATION   %%%
%%%%%%%%%%%%%%%%%%%%%%%%

strTransPlate = Pathname;
if ~any(fileattrib(strTransPlate))
    error('Could not find reference plate');
end

if handles.Current.SetBeingAnalyzed == 1
    makeExternalBuffer(handles, strTransPlate, TextToFind);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY CALCULATIONS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% Determines which cycle is being analyzed.
SetBeingAnalyzed = handles.Current.SetBeingAnalyzed;

%%% Remove slashes entries with N/A or no filename from the input,
%%% i.e., store only valid entries
tmp1 = {};
tmp2 = {};
for n = 1:numel(ImageTransName)
    if ~strcmp(TextToFind{n}, '/') && ~strcmp(ImageTransName{n}, '/')
        tmp1{end+1} = TextToFind{n};
        tmp2{end+1} = ImageTransName{n};
    end
end
TextToFind = tmp1;
ImageTransName = tmp2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FIRST CYCLE FILE HANDLING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% Extracting the list of files to be analyzed occurs only the first time
%%% through this module.

% get image paths
[Pathname] = getSubFoldersFromTransPlate(strTransPlate);

% check wheter you are overwriting orig images with trans images
if SetBeingAnalyzed == 1
    
    for i = 1:length(ImageTransName)
        if isfield(handles.Pipeline,ImageTransName{i})
            error(['Image processing was cancelled in the ', ModuleName, ' module because you are trying to load two sets of images with the same name (e.g. OrigBlue). The last set loaded will always overwrite the first set and make it obselete. Please remove one of these modules.']);
        end
    end   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LOADING IMAGES EACH TIME %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%% for all cycles, new
FileNames = cell(numel(ImageTransName),1);

for ix = 1:numel(ImageTransName)
[strCorrespondingImage_trans, couldFindSameSite_image] = getFileNameViaReferenceFile(handles, ImageCisName{ix}, Pathname);

if couldFindSameSite_image
    fp = fullfile(Pathname, strCorrespondingImage_trans);
    
    if DoIllCorrIm
        LoadedImage = imread_illumination_corrected(fp,0);
    else
        LoadedImage = imread(fp);
    end
    
    fieldname = ['Filename', ImageTransName{ix}];
    handles.Pipeline.(fieldname){SetBeingAnalyzed} = strCorrespondingImage_trans;
%     handles.Pipeline.(ImageTransName{ix}){SetBeingAnalyzed} = LoadedImage; % this line was used when combine seg in pipeline 2. 
    if range(LoadedImage(:)) > 1 % this setting is used when using for correlation measurements
        LoadedImage = LoadedImage./double(uint16(inf));
    end
    handles.Pipeline.(ImageTransName{ix}) = LoadedImage; % this setting is used when using for correlation measurements
    FileNames{ix} = strCorrespondingImage_trans;
else
    error('No corresponding image could be found')
end

end

%%%%%%%%%%%%%%%%%%%%%%%
%%% DISPLAY RESULTS %%%
%%%%%%%%%%%%%%%%%%%%%%%
drawnow

ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
if any(findobj == ThisModuleFigureNumber);
    if handles.Current.SetBeingAnalyzed == handles.Current.StartingImageSet
        CPresizefigure('','NarrowText',ThisModuleFigureNumber)
    end
    for n = 1:length(ImageTransName)
        %%% Activates the appropriate figure window.
        currentfig = CPfigure(handles,'Text',ThisModuleFigureNumber);
        if iscell(ImageTransName)
            TextString = [ImageTransName{n},': ',FileNames{n}];
        else
            TextString = [ImageTransName,': ',FileNames];
        end
        uicontrol(currentfig,'style','text','units','normalized','fontsize',handles.Preferences.FontSize,'HorizontalAlignment','left','string',TextString,'position',[.05 .85-(n-1)*.15 .95 .1],'BackgroundColor',[.7 .7 .9])
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SAVE DATA TO HANDLES %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% NOTE: The structure for filenames and pathnames will be a cell array of cell arrays

%%% First, fix feature names and the pathname
PathNames = cell(1,length(ImageTransName));
FileNamesText = cell(1,length(ImageTransName));
PathNamesText = cell(1,length(ImageTransName));
for n = 1:length(ImageTransName)
    PathNames{n} = Pathname;
    FileNamesText{n} = [ImageTransName{n}];
    PathNamesText{n} = [ImageTransName{n}];
end

%%% Since there may be several load/save modules in the pipeline which all
%%% write to the handles.Measurements.Image.FileName field, we store
%%% filenames in an "appending" style. Here we check if any of the modules
%%% above the current module in the pipeline has written to
%%% handles.Measurements.Image.Filenames. Then we should append the current
%%% filenames and path names to the already written ones. If this is the
%%% first module to put anything into the handles.Measurements.Image
%%% structure, then this section is skipped and the FileNamesText fields
%%% are created with their initial entry coming from this module.

if  isfield(handles,'Measurements') && isfield(handles.Measurements,'Image') &&...
        isfield(handles.Measurements.Image,'FileNamesTrans') && length(handles.Measurements.Image.FileNames) == SetBeingAnalyzed
    % Get existing file/path names. Returns a cell array of names
    ExistingFileNamesText = handles.Measurements.Image.FileNamesText;
    ExistingFileNames     = handles.Measurements.Image.FileNames{SetBeingAnalyzed};
    ExistingPathNamesText = handles.Measurements.Image.PathNamesText;
    ExistingPathNames     = handles.Measurements.Image.PathNames{SetBeingAnalyzed};
    % Append current file names to existing file names
    FileNamesText = cat(2,ExistingFileNamesText,FileNamesText);
    FileNames     = cat(2,ExistingFileNames,FileNames);
    PathNamesText = cat(2,ExistingPathNamesText,PathNamesText);
    PathNames     = cat(2,ExistingPathNames,PathNames);
end

%%% Write to the handles.Measurements.Image structure
handles.Measurements.Image.TransFileNamesText                   = FileNamesText;
handles.Measurements.Image.TransFileNames(SetBeingAnalyzed)         = {FileNames};
handles.Measurements.Image.TransPathNamesText                   = PathNamesText;
handles.Measurements.Image.TransPathNames(SetBeingAnalyzed)         = {PathNames};
end

%% Subfunctions
% to find corresponding image in trans
function makeExternalBuffer(handles, strTransPlate, filterForImagesOfTrans)
% note that preferentially all data would be stored in pipeline. However
% storing custom fields in handles.pipelines can lead to various errors.
% While storing data in handles.measurements is technically possible there
% would be various mistakes occuring if order of sites becomes rearranged
% in batch files. Since there will always only be one segmentation of a
% given name, there can be a buffer, that links to that name, and which is
% overwritten, if a new pipeline is made


[TiffFolder_trans] = getSubFoldersFromTransPlate(strTransPlate);

[a, ~, ~] = fileparts(TiffFolder_trans);
[~, b, ~] = fileparts(a);

% Images
eB.strImages_trans = getFilesAndDirectories(TiffFolder_trans, filterForImagesOfTrans);
logX = cellfun(@(x) isempty(x), strfind(eB.strImages_trans, '.png'));
eB.strImages_trans(logX) = [];

[eB.Row_image_trans, eB.intColumn_image_trans, eB.intImagePosition_image_trans, ~, ~, eB.intChannelName_image_trans] = cellfun(@(x) MetaFromImageName(x), eB.strImages_trans, 'UniformOutput', true);


% save name
outDir = handles.Current.DefaultOutputDirectory;
ex = ['SimpleMulti_' 'ImagesTrans_' b '.mat'];
strBufferFile = fullfile(outDir, ex );


if any(fileattrib(strBufferFile))
    [~, ex] = fileparts(strBufferFile);
    fprintf([ex ' already exists. It will be overwritten with current one.\n']);
end

save(strBufferFile, 'eB');
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
function [TiffFolder_trans] = getSubFoldersFromTransPlate(strTransPlate)


if any(strfind(strTransPlate, [filesep 'TIFF']))
    error('Reference directory must refer to a plate folder, not the TIFF folder');
end

% if any(strfind(strTransPlate, [filesep 'SEGMENTATION']))
%     error('Reference directory must refer to a plate folder, not the SEGMENTATION folder');
% end

TiffFolder_trans = fullfile(strTransPlate, 'TIFF');
if ~any(fileattrib(TiffFolder_trans))
    error('Could not find TIFF folder of other plate');
end

% SegmentationFolder_trans = fullfile(strTransPlate, 'SEGMENTATION');
% if ~any(fileattrib(SegmentationFolder_trans))
%     error('Could not find SEGMENTATION folder of other plate');
% end


end
function [strCorrespondingImage_trans, couldFindSameSite_image] = getFileNameViaReferenceFile(handles, OrigImageName, TiffFolder_trans)
% Avoid repeated, independent access to buffer file by reading segmentation
% and images after same loading event. Note that the code thus becomes more
% ugly and less readable.


SetBeingAnalyzed = handles.Current.SetBeingAnalyzed;
% size(handles.Pipeline.(['Filename' OrigImageName])) %[GG20160419] old
% SetBeingAnalyzed %[GG20160419] old
% ImageName_cis = handles.Pipeline.(['Filename' OrigImageName]){SetBeingAnalyzed}; %[GG20160419] old
ImageName_cis = handles.Pipeline.(['Filename' OrigImageName]){SetBeingAnalyzed}; %[GG20160419] new

outDir = handles.Current.DefaultOutputDirectory;

[a, ~, ~] = fileparts(TiffFolder_trans);
[~, b, ~] = fileparts(a);

ex = ['SimpleMulti_' 'ImagesTrans_' b '.mat'];
strBufferFile = fullfile(outDir, ex );

if ~any(fileattrib(strBufferFile))
    error('Could not find buffer file for object, that should be created during first cycle');
else
    load(strBufferFile);
end

[Row_cis, intColumn_cis, intImagePosition_cis, ~, ~, intChannelNumber_cis] = MetaFromImageName(ImageName_cis);

% get corresponding image from trans
correspondsToSameSite_image = eB.intChannelName_image_trans ==  intChannelNumber_cis & eB.Row_image_trans == Row_cis & eB.intColumn_image_trans == intColumn_cis & eB.intImagePosition_image_trans == intImagePosition_cis;
if sum(correspondsToSameSite_image) == 0
    couldFindSameSite_image = false;
    strCorrespondingImage_trans = '';
elseif sum(correspondsToSameSite_image) == 1;
    couldFindSameSite_image = true;
    strCorrespondingImage_trans = eB.strImages_trans{correspondsToSameSite_image};
else
    error('Could not unambiguously find corresponding image of other dataset. Please set more stringent filters.');
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
function corr_image = imread_illumination_corrected(strImage, shallCacheCorrection)
% reads images and applies illumination correction

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
function [matMeanImage matStdImage hasIlluminationCorrection] = getIlluminationReference(strBatchDir,iChannel,cacheInRam)
if nargin < 3
    cacheInRam = false;
end

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
function dat = cacheForIlluminationStats(strFileName)
% Initialize Persistent variables for caching
persistent CachedMeasurments;
persistent OriginalPathOfChached;

if isempty(CachedMeasurments)
    CachedMeasurments = cell(0);
end

if isempty(OriginalPathOfChached)
    OriginalPathOfChached = cell(0);
end

nStrFileName = strFileName; % npc to ensure that each file only stored once in cache once

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

