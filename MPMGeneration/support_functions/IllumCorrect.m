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

%% subfunciton
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

