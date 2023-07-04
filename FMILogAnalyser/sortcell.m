function cell = sortcell(Cell)
cell = cellfun(@(x)sort(x),Cell,'UniformOutput',false);
end