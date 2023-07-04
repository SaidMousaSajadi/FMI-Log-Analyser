function L = CheckCell(CC,Vector)
  for pop = 1:size(CC,1)
    Po = CC(pop,:) ;
    Wow = [] ;
    for kok = 1:length(Po)
      Wow = [Wow , reshape(cell2mat(Po{kok}),1,[])] ;
    end
    if ~isempty(setdiff(Vector , unique(Wow)))
      L(pop,1) = true ;
    else
      L(pop,1) = false ;
    end
  end
end
