function handles = AlignTransImages(handles)

% Help for the AlignTransImages module:
% Category: Image Processing
%
% SHORT DESCRIPTION:
% Aligns images from different 4i cycles using Fast Fourier Transform.
% *************************************************************************
%
% $Revision: 4558 $

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = What did you call the reference image on which you want to align all other images?
%infotypeVAR01 = imagegroup
%inputtypeVAR01 = popupmenu
RefCycleImageName = char(handles.Settings.VariableValues{CurrentModuleNum,1});

%textVAR02 = How do you call ref. image for first alignement?
%infotypeVAR02 = imagegroup
%inputtypeVAR02 = popupmenu
RefImageName{1} = char(handles.Settings.VariableValues{CurrentModuleNum,2});

%textVAR03 = How do you call image for first alignment
%infotypeVAR03 = imagegroup
%inputtypeVAR03 = popupmenu
ShiftImageName{1} = char(handles.Settings.VariableValues{CurrentModuleNum,3});

%textVAR04 = What do you want to call the first aligned image?
%infotypeVAR04 = imagegroup indep
%defaultVAR04 = /
AlignedImageName{1} = char(handles.Settings.VariableValues{CurrentModuleNum,4});

%textVAR05 = How do you call ref. image for second alignement?
%infotypeVAR05 = imagegroup
%inputtypeVAR05 = popupmenu
RefImageName{2} = char(handles.Settings.VariableValues{CurrentModuleNum,5});

%textVAR06 = How do you call image for second alignment
%infotypeVAR06 = imagegroup
%inputtypeVAR06 = popupmenu
ShiftImageName{2} = char(handles.Settings.VariableValues{CurrentModuleNum,6});

%textVAR07 = What do you want to call the second aligned image?
%infotypeVAR07 = imagegroup indep
%defaultVAR07 = /
AlignedImageName{2} = char(handles.Settings.VariableValues{CurrentModuleNum,7});

%textVAR08 = How do you call ref. image for third alignement?
%infotypeVAR08 = imagegroup
%inputtypeVAR08 = popupmenu
RefImageName{3} = char(handles.Settings.VariableValues{CurrentModuleNum,8});

%textVAR09 = How do you call image for third alignment
%infotypeVAR09 = imagegroup
%inputtypeVAR09 = popupmenu
ShiftImageName{3} = char(handles.Settings.VariableValues{CurrentModuleNum,9});

%textVAR10 = What do you want to call the third aligned image?
%infotypeVAR10 = imagegroup indep
%defaultVAR10 = /
AlignedImageName{3} = char(handles.Settings.VariableValues{CurrentModuleNum,10});

%textVAR11 = What is the maximal pixel shift you expect?
%defaultVAR11 = 80
max_shift = str2double(handles.Settings.VariableValues{CurrentModuleNum,11});

%textVAR12 = How many time should we try to align the images?
%defaultVAR12 = 5
max_iterations = str2double(handles.Settings.VariableValues{CurrentModuleNum,12});

%%%VariableRevisionNumber = 2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY CALCULATIONS & FILE HANDLING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Determines which cycle is being analyzed.
SetBeingAnalyzed = handles.Current.SetBeingAnalyzed;

drawnow

%% check how many images we have to align
linx_cycles_to_align = find(cellfun(@(x) isempty(x), strfind(AlignedImageName, '/')));
num_cycles = numel(linx_cycles_to_align);

%% retriev images
RefImage = cell(num_cycles,1);
ShiftImage = cell(num_cycles,1);

for ix = 1:num_cycles
    t = linx_cycles_to_align(ix);
    RefImage{ix} = handles.Pipeline.(RefImageName{t}){SetBeingAnalyzed}; 
    ShiftImage{ix} = handles.Pipeline.(ShiftImageName{t}){SetBeingAnalyzed};
%     RefImage{ix} = CPretrieveimage(handles,RefImageName{t},ModuleName,'DontCheckColor','CheckScale');
%     ShiftImage{ix} = CPretrieveimage(handles,ShiftImageName{t},ModuleName,'DontCheckColor','CheckScale');
end

%% Align all ref images to RefCycle
RefCycleImage = CPretrieveimage(handles,RefCycleImageName,ModuleName,'DontCheckColor','CheckScale');

NSWE_Cis = NaN(num_cycles,4);
NSWE_Trans = NaN(num_cycles,4);

for ix = 1:num_cycles
    [NSWE_Cis(ix,:), NSWE_Trans(ix,:)] = getRobustAlignmentAtLowIntensity(RefCycleImage, RefImage{ix}, max_shift, max_iterations);
end


%% Shift images
shiftstack = zeros(size(RefCycleImage,1),size(RefCycleImage,2),num_cycles);
for ix = 1:size(NSWE_Trans,1) % loop over cycles
    shiftstack(NSWE_Trans(ix,1):NSWE_Trans(ix,2), NSWE_Trans(ix,3):NSWE_Trans(ix,4),ix) = ...
        ShiftImage{ix}(NSWE_Cis(ix,1):NSWE_Cis(ix,2), NSWE_Cis(ix,3):NSWE_Cis(ix,4));
end

%% do some name clean up before saving
AlignedImageName = AlignedImageName(linx_cycles_to_align);
FileNames = AlignedImageName;

%% Check whether RefCycleImage is between 0 and 1
if range(RefCycleImage(:)) < 1
    for ix = 1:size(shiftstack,3) % loop over cycles
        shiftstack(:,:,ix) = shiftstack(:,:,ix)./double(uint16(inf)); % assuming we work with 16bit images
    end
end
    
%% save imges
for ix = 1:num_cycles
    
    fieldname = ['Filename', AlignedImageName{ix}];
    handles.Pipeline.(fieldname){SetBeingAnalyzed} = AlignedImageName{ix};
    handles.Pipeline.(FileNames{ix}) = shiftstack(:,:,ix);
%     FileNames{ix} = strCorrespondingImage_trans;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SAVE DATA TO HANDLES %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% NOTE: The structure for filenames and pathnames will be a cell array of cell arrays

%%% First, fix feature names and the pathname
% PathNames = cell(1,length(AlignedImageName));
FileNamesText = cell(1,length(AlignedImageName));
% PathNamesText = cell(1,length(AlignedImageName));
for n = 1:length(AlignedImageName)
%     PathNames{n} = Pathname;
    FileNamesText{n} = [AlignedImageName{n}];
%     PathNamesText{n} = [AlignedImageName{n}];
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
        isfield(handles.Measurements.Image,'TransFileNamesText') && length(handles.Measurements.Image.FileNames) == SetBeingAnalyzed
    % Get existing file/path names. Returns a cell array of names
    ExistingFileNamesText = handles.Measurements.Image.FileNamesText;
    ExistingFileNames     = handles.Measurements.Image.FileNames{SetBeingAnalyzed};
    ExistingPathNamesText = handles.Measurements.Image.PathNamesText;
    ExistingPathNames     = handles.Measurements.Image.PathNames{SetBeingAnalyzed};
    % Append current file names to existing file names
    FileNamesText = cat(2,ExistingFileNamesText,FileNamesText);
    FileNames     = cat(2,ExistingFileNames,FileNames);
%     PathNamesText = cat(2,ExistingPathNamesText,PathNamesText);
%     PathNames     = cat(2,ExistingPathNames,PathNames);
end

%%% Write to the handles.Measurements.Image structure
handles.Measurements.Image.TransFileNamesText                   = FileNamesText;
handles.Measurements.Image.TransFileNames(SetBeingAnalyzed)         = {FileNames};
% handles.Measurements.Image.TransPathNamesText                   = PathNamesText;
% handles.Measurements.Image.TransPathNames(SetBeingAnalyzed)         = {PathNames};
end
%%%%%%%%%%%%%%%%%%%%
%%% SUBFUNCTIONS %%%
%%%%%%%%%%%%%%%%%%%%

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
function [NSWE_Cis, NSWE_Trans] = getNSWEofOverlappingImageparts(Image_Cis, Image_Trans)
% upon providing two images (Image_Cis, Image_Trans), The top row (N -orth),
% lowest row (S -outh),  left column (W -est) and right column (E -east) of
% the overlapping region are returned;

[imageHeight, imageWidth] = size(Image_Cis);

fourierTrans = fft2(Image_Trans);
fourierCis = fft2(Image_Cis);
[registrationOutput] = dftregistration(fourierTrans,fourierCis,1);

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
function [output Greg] = dftregistration(buf1ft,buf2ft,usfac)
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
if exist('roff')~=1, roff=0; end
if exist('coff')~=1, coff=0; end
if exist('usfac')~=1, usfac=1; end
if exist('noc')~=1, noc=nc; end
if exist('nor')~=1, nor=nr; end
% Compute kernels and obtain DFT by matrix products
kernc=exp((-i*2*pi/(nc*usfac))*( ifftshift([0:nc-1]).' - floor(nc/2) )*( [0:noc-1] - coff ));
kernr=exp((-i*2*pi/(nr*usfac))*( [0:nor-1].' - roff )*( ifftshift([0:nr-1]) - floor(nr/2)  ));
out=kernr*in*kernc;
end
