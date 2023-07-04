function cell = uniquecell(Cell)

strcell = cellfun(@(x) num2str(x(:)'),Cell,'UniformOutput',false);
[~,index] = unique(strcell);
cell = Cell(index,1) ;
end