function Clusters = WhatACoolYourPhase(X,G,Kindex,cycle,Quality,Method)
  % Need Functions: AllPath,sortcell,uniquecell,MergeData,fitfunction,ClusterFinder,MergeClusters,sortcellelements,uniquecells,NumberDistributionStates,combvec_oct
  GI = logical(G) ;

  % Isolate Edge
  ISO = find(sum(GI)==0) ;

  % Collect All Connection States
  Fuji = [] ;
  for j = 1:max(Kindex)
      AccpetableConnection = [j ; find(GI(:,j))] ;
      for k = 2:length(AccpetableConnection)
          PathsOfGI = AllPath(GI,j,AccpetableConnection(k)) ;
          for i = 1:length(PathsOfGI)
              if ismember(PathsOfGI{i},AccpetableConnection)
                  Fuji = [Fuji ; PathsOfGI(i)] ;
              end
          end
      end
  end
  % Sort and Uniquelization Inner Connections
  try
      Fuji = sortcell(Fuji) ;
      Fuji = uniquecell(Fuji) ;
  catch
  end

  for i = 1:size(Fuji,1)
      x = MergeData(X,Kindex,Fuji{i,1}) ;
      Fuji{i,2} = fitfunction(x(:,1),x(:,2),20*ones(size(x(:,1))),cycle,Quality,Method) ;
  end

  if ~isempty(Fuji)
    % Remove Worst Cases
    Fuji(cell2mat(cellfun(@(x) x==0,Fuji(:,2),'UniformOutput',false)),:) = [] ;
  else
    CoClus = [] ;
    for dior = 1:numel(ISO)
      CoClus = [CoClus ; {ISO(dior)}] ;
    end
    Clusters = {CoClus} ;
  end

  if ~isempty(Fuji)
    % Connection Sizes
    Fuji = [Fuji cellfun(@length,Fuji(:,1),'UniformOutput',false)] ;

    % Add Isolated
    Fuji = [Fuji cellfun(@(x) sort([x ISO]) ,Fuji(:,1),'UniformOutput',false)] ;

    % Connection Sizes With Isolated
    Fuji = [Fuji cellfun(@length,Fuji(:,4),'UniformOutput',false)] ;

    % Sort Cell by Connection Sizes
    Fuji = sortrows(Fuji,2) ;
    Fuji = sortrows(Fuji,5) ;
    Fuji = flipud(Fuji) ;

    % New Clusters must be include:
    Fuji = [Fuji cellfun(@(x) setdiff(1:max(Kindex),x),Fuji(:,4),'UniformOutput',false)] ;


    %% Collect all possible states for clusters
    tic
    Shuji = ClusterFinder(Fuji,max(Kindex)) ;
    toc

    %% Sorting all states and selecting unique states
    Clusters = MergeClusters(Fuji,Shuji,ISO) ;
    Clusters = sortcellelements(Clusters) ;
    tic
    Clusters = uniquecells(Clusters) ;
    toc
  else
    CoClus = [] ;
    for dior = 1:numel(ISO)
      CoClus = [CoClus ; {ISO(dior)}] ;
    end
    Clusters = {CoClus} ;
  end

end
