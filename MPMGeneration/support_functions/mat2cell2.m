function cellOutput = mat2cell2(matInput)

if isnumeric(matInput)
    cellOutput = arrayfun(@(x) {x}, matInput,'UniformOutput',true);
end
    