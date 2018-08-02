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
