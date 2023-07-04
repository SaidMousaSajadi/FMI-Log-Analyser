function Cell = Cluster2Cell(Clusters)
  for i = 1:size(Clusters,1)
    BUD = [] ;
    for j = 1:size(Clusters{i,1},1)
      if j < size(Clusters{i,1},1)
        Tag = ['[' num2str(Clusters{i,1}{j,:}) '],'] ;
        BUD = strcat(BUD, Tag) ;
      else
        Tag = ['[' num2str(Clusters{i,1}{j,:}) ']'] ;
        BUD = strcat(BUD, Tag) ;
      end
    end
    CStr{i,1} = BUD ;
    eval(['Struct.C' num2str(i) ' = BUD ;'])
  end
  Cell = Struct ;
end
