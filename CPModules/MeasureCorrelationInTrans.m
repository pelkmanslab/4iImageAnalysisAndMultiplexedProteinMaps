function handles = MeasureCorrelationInTrans(handles)

% Help for the Measure Correlation in Trans  module:
% Category: Measurement
%
% SHORT DESCRIPTION:
% Measures the correlation between intensities in different images (e.g.
% different color channels) coming from different 4i cycles on a pixel by pixel basis, within identified
% objects or across an entire image.
% This Module is heavily based on the module MeasureCorrelationInTrans
% written by Anne E. Carpenter and collegues, developped at the Whitehead Institure for Biomedical Research, Copyright 2003,2004,2005.
% *************************************************************************
%
% $Revision: 4557 $

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = Choose Cis image types to measure correlation:
%choiceVAR01 = Do not use
%infotypeVAR01 = imagegroup
ImageName{1} = char(handles.Settings.VariableValues{CurrentModuleNum,1});
%inputtypeVAR01 = popupmenu

%textVAR02 = Choose Cis image types to measure correlation:
%choiceVAR02 = Do not use
%infotypeVAR02 = imagegroup
ImageName{2} = char(handles.Settings.VariableValues{CurrentModuleNum,2});
%inputtypeVAR02 = popupmenu

%textVAR03 = Choose Trans image types to measure correlation:
%choiceVAR03 = Do not use
%infotypeVAR03 = imagegroup
TransImageName{1} = char(handles.Settings.VariableValues{CurrentModuleNum,3});
%inputtypeVAR03 = popupmenu

%textVAR04 = Choose Trans image types to measure correlation:
%choiceVAR04 = Do not use
%infotypeVAR04 = imagegroup
TransImageName{2} = char(handles.Settings.VariableValues{CurrentModuleNum,4});
%inputtypeVAR04 = popupmenu

%textVAR05 = Choose objects within which to measure the correlations (Choosing Image will measure correlations across the entire images)
%choiceVAR05 = Do not use
%choiceVAR05 = Image
%infotypeVAR05 = objectgroup
ObjectName{1} = char(handles.Settings.VariableValues{CurrentModuleNum,5});
%inputtypeVAR05 = popupmenu

%textVAR06 =
%choiceVAR06 = Do not use
%choiceVAR06 = Image
%infotypeVAR06 = objectgroup
ObjectName{2} = char(handles.Settings.VariableValues{CurrentModuleNum,6});
%inputtypeVAR06 = popupmenu

%textVAR07 =
%choiceVAR07 = Do not use
%choiceVAR07 = Image
%infotypeVAR07 = objectgroup
ObjectName{3} = char(handles.Settings.VariableValues{CurrentModuleNum,7});
%inputtypeVAR07 = popupmenu

%textVAR08 = Choose trans objects within which to measure the correlations (Choosing Image will measure correlations across the entire images)
%choiceVAR08 = Do not use
%choiceVAR08 = Image
%infotypeVAR08 = objectgroup
TransObjectName{1} = char(handles.Settings.VariableValues{CurrentModuleNum,8});
%inputtypeVAR08 = popupmenu

%textVAR09 =
%choiceVAR09 = Do not use
%choiceVAR09 = Image
%infotypeVAR09 = objectgroup
TransObjectName{2} = char(handles.Settings.VariableValues{CurrentModuleNum,9});
%inputtypeVAR09 = popupmenu

%textVAR10 =
%choiceVAR10 = Do not use
%choiceVAR10 = Image
%infotypeVAR10 = objectgroup
TransObjectName{3} = char(handles.Settings.VariableValues{CurrentModuleNum,10});
%inputtypeVAR10 = popupmenu

%textVAR11 = Additional save name, in case several measurements of
%different cycles are done
%defaultVAR11 = None
AdditionalSaveName = char(handles.Settings.VariableValues{CurrentModuleNum,11});


%textVAR12 = Specify filter size for image processing:
%defaultVAR12 = 10
filter_size = str2double(char(handles.Settings.VariableValues{CurrentModuleNum,12}));
%%%VariableRevisionNumber = 3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY CALCULATIONS & FILE HANDLING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% check whether same number of trans and normals objects
if numel(TransObjectName) ~= numel(ObjectName)
    error('Different number of trans and normal objects');
end

drawnow

%%% Get the images
ImageCount = 0;
for ImageNbr = 1:2
    if ~strcmp(ImageName{ImageNbr},'Do not use')
        ImageCount = ImageCount + 1;
        try
            %%% Checks whether image has been loaded.
            Image{ImageCount} = CPretrieveimage(handles,ImageName{ImageNbr},ModuleName,'MustBeGray','DontCheckScale'); %#ok Ignore MLint
            tmpImageName{ImageCount} = ImageName{ImageNbr}; %#ok Ignore MLint
            
            TransImage{ImageCount} = CPretrieveimage(handles,TransImageName{ImageNbr},ModuleName,'MustBeGray','DontCheckScale'); %#ok Ignore MLint
            tmpTransImageName{ImageCount} = TransImageName{ImageNbr}; %#ok Ignore MLint
            
        catch 
            error(['Image processing was canceled in the ', ModuleName, ' module because there was a problem loading the image you called ', ImageName{ImageNbr}, '.'])
        end
    end
end
%%% Get rid of 'Do not use' in the ImageName cell array so we don't have to care about them later.
ImageName = tmpImageName;   
TransImageName = tmpTransImageName;   


%%% Check so that at least two images have been entered
if numel(ImageName) < 1 || numel(TransImageName) < 1
    error(['Image processing was canceled in the ', ModuleName, ' module because at least two image types must be chosen.'])
end

%%% Get the masks of segmented objects
ObjectNameCount = 0;
for ObjectNameNbr = 1:3
    if ~strcmp(ObjectName{ObjectNameNbr},'Do not use')
        ObjectNameCount = ObjectNameCount + 1;
        tmpObjectName{ObjectNameCount} = ObjectName{ObjectNameNbr}; %#ok Ignore MLint
        tmpTransObjectName{ObjectNameCount} = TransObjectName{ObjectNameNbr}; %#ok Ignore MLint
        if ~strcmp(ObjectName{ObjectNameNbr},'Image')
            %%% Retrieves the label matrix image that contains the
            %%% segmented objects which will be used as a mask.
            LabelMatrixImage{ObjectNameCount} = CPretrieveimage(handles,['Segmented', ObjectName{ObjectNameNbr}],ModuleName,'MustBeGray','DontCheckScale'); %#ok Ignore MLint
            LabelMatrixImageTrans{ObjectNameCount} = CPretrieveimage(handles,['Segmented', TransObjectName{ObjectNameNbr}],ModuleName,'MustBeGray','DontCheckScale'); %#ok Ignore MLint
        else
            LabelMatrixImage{ObjectNameCount} = ones(size(Image{1}));        % Use mask of ones to indicate that the correlation should be calcualted for the entire image
            LabelMatrixImageTrans{ObjectNameCount} = ones(size(Image{1}));        % Use mask of ones to indicate that the correlation should be calcualted for the entire image
        end
    end
end
%%% Get rid of 'Do not use' in the ObjectName cell array so we don't have to care about them later.
ObjectName = tmpObjectName;
TransObjectName = tmpTransObjectName;

%%% Check so that at least one object type have been entered
if numel(ObjectName) < 1 || numel(TransObjectName) < 1
    error(['At least one object type must be entered in the ',ModuleName,' module.'])
end

%%%%%%%%%%%%%%%%%%%%%%
%%% IMAGE ANALYSIS %%%
%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% Produce feature names for all pairwise image combinations
CorrelationFeatures = {};
Correlation2DFeatures = {};

for i = 1:numel(ImageName)
    for j = 1:numel(TransImageName)
        CorrelationFeatures{end+1} = [ImageName{i},'_',TransImageName{j}];
        Correlation2DFeatures{end+1} = [ImageName{i},'_',TransImageName{j}];
    end
end

%%% For each object type and for each segmented object, calculate the correlation between all combinations of images
for ObjectNameNbr = 1:ObjectNameCount
    
    
    if any(size(Image{i}) ~= size(LabelMatrixImage{ObjectNameNbr})) || any(size(TransImage{i}) ~= size(LabelMatrixImageTrans{ObjectNameNbr}))
        error(['Image processing was canceled in the ', ModuleName, ' module. The size of the image you want to measure is not the same as the size of the image from which the ',ObjectName{ObjectNameNbr},' objects were identified.'])
    end
    
    if max(LabelMatrixImage{ObjectNameNbr}(:)) ~= max(LabelMatrixImageTrans{ObjectNameNbr}(:))
        error(['Image processing was canceled in the ', ModuleName, ' module. The number of ', ObjectName{ObjectNameNbr}, ' in the cis and the trans image do not correspond.'])
    end
    
    %%% Calculate the correlation in all objects for all pairwise image combinations
    NbrOfObjects = max(LabelMatrixImage{ObjectNameNbr}(:));          % Get number of segmented objects
    Correlation = zeros(NbrOfObjects,length(CorrelationFeatures));   % Pre-allocate memory'
    Correlation2 = zeros(NbrOfObjects,length(CorrelationFeatures));   % Pre-allocate memory'
    uq_objects = unique(LabelMatrixImage{ObjectNameNbr}(:));         % detect object ID
    uq_objects(1) = [];                                              % delete first ID as it is 0 == background.
    for ObjectNbr = 1:NbrOfObjects                                   % Loop over objects
        FeatureNbr = 1;                                              % Easiest way to keep track of the feature number, i.e. which combination of images
        for i = 1:numel(ImageName)                                     % Loop over all combinations of images
            for j = 1:numel(TransImageName)
                index_cis = LabelMatrixImage{ObjectNameNbr} == uq_objects(ObjectNbr);   % Get the indexes for the this object number
                index_trans = LabelMatrixImageTrans{ObjectNameNbr} == uq_objects(ObjectNbr);
                if sum(index_cis(:)) ~= sum(index_trans(:))
                    CorrelationForCurrentObject = NaN;
                    CorrelationForCurrentObject2 = NaN;
                elseif isempty(index_cis) || numel(index_cis) == 1 || isempty(index_trans) || numel(index_trans) == 1 % If the object does not exist in the label matrix or is only one pixel, the correlation calculation will fail, so we assign the correlation to be NaN.
                    CorrelationForCurrentObject = NaN;
                    CorrelationForCurrentObject2 = NaN;
                else
                    f = ones(filter_size, filter_size);
                    im_cis = imfilter(Image{i}(index_cis),f);
                    im_trans = imfilter(TransImage{j}(index_trans), f);
                    CorrelationForCurrentObject = corr(im_cis, im_trans);            % Get the values for these indexes in the images and calculate the correlation
                    CorrelationForCurrentObject2 = corr2(im_cis, im_trans);            % Get the values for these indexes in the images and calculate the correlation
                end
                Correlation(ObjectNbr,FeatureNbr) = CorrelationForCurrentObject; % Store the correlation
                Correlation2(ObjectNbr,FeatureNbr) = CorrelationForCurrentObject2; % Store the correlation
                FeatureNbr = FeatureNbr + 1;
            end
        end
    end
    %%% Store the correlation and slope measurements
        handles.Measurements.(ObjectName{ObjectNameNbr}).(['Correlation_', AdditionalSaveName, 'Features']) = CorrelationFeatures;
        handles.Measurements.(ObjectName{ObjectNameNbr}).(['Correlation_', AdditionalSaveName])(handles.Current.SetBeingAnalyzed) = {Correlation};

        handles.Measurements.(ObjectName{ObjectNameNbr}).(['Correlation2D_', AdditionalSaveName, 'Features']) = CorrelationFeatures;
        handles.Measurements.(ObjectName{ObjectNameNbr}).(['Correlation2D_', AdditionalSaveName])(handles.Current.SetBeingAnalyzed) = {Correlation2};
end


%%%%%%%%%%%%%%%%%%%%%%%
%%% DISPLAY RESULTS %%%
%%%%%%%%%%%%%%%%%%%%%%%
drawnow

ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
if any(findobj == ThisModuleFigureNumber)
    if CPisHeadless == false
        %%% Activates the appropriate figure window.
        for i = 1:numel(ObjectNameCount)
            CPfigure(handles,'Image',ThisModuleFigureNumber);
            
            subplot(2, 2, 1)
            a = CPlabel2rgb(handles,LabelMatrixImage{i});
            CPimagesc(a,handles);
            axis image
            title('Cis Segmentation')
            
            subplot(2, 2, 2)
            a = CPlabel2rgb(handles,LabelMatrixImageTrans{i});
            CPimagesc(a,handles);
            axis image
            title('Trans Segmentation')
            
            subplot(2, 2, 3)
            CPimagesc(Image{i},handles);
            colormap('JET')
            axis image
            title('Cis Image')
            
            subplot(2, 2, 4)
            CPimagesc(TransImage{i},handles);
            colormap('JET')
            axis image
            title('Trans Image')
        end
    end
end

end
