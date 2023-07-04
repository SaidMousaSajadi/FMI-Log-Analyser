function Spart = MergeCell(CC)
  Spart = [] ;
  for pop = 1:size(CC,1)
    Po = CC(pop,:) ;
    Dow = [] ;
    for kok = 1:length(Po)
      Dow = [Dow ; Po{kok}] ; % Merge
    end
    Spart{pop,1} = Dow;
  end
end
