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