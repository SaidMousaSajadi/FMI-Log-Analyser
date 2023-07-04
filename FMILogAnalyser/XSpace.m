function [Caper , Hash] = XSpace(Fuji,MaximumCluster)

  for io = 1:size(Fuji,1)
    ko = MaximumCluster-Fuji{io,5} ;
    Caper{io,1} = Fuji{io,6}' ;
    for tiny = 2:ko
      List = Fuji(find(cell2mat(Fuji(:,3)) == tiny),1).' ;
      B = Fuji(io,6);
      Indo = cell2mat(cellfun(@(x) all(ismember(x,B{1})),List(1,:),'UniformOutput',false)) ;
      List = List.' ;
      List(Indo==0,:) = [] ;
      Caper{io,tiny} = List ;
    end
    Caper{io,end} = Fuji{io,6} ;
  end
  for io = 1:size(Caper,1)
    Hash(io,:) = cell2mat(cellfun(@(x) size(x,1),Caper(io,:),'UniformOutput',false)) ;
  end
end
