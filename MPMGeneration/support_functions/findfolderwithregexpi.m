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