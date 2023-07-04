function Shuji = ClusterFinder(Fuji,NumberOfClust)
  % Need Functions: NumberDistributionStates,combvec_oct,XSpace,Double2Cell,CheckCell,MergeCell


##  NumberOfClust = max(Kindex) ;

  % Creating x space [Flexible Fuji]
  [Caper , Hash] = XSpace(Fuji,NumberOfClust) ;

  for io = 1:size(Fuji,1) % Fuji

    ko = NumberOfClust - Fuji{io,5} ;
    Ate = NumberDistributionStates(size(Caper,2),ko) ; % Base , Target
    Holo = [] ;
    for nano = 1:size(Ate,1) % Ate
      clear CC Spart
      if all(Hash(io,:) >= Ate(nano,:))

        if nano == size(Ate,1) & size(Ate,1) > 1 % Better this to general Form
          Spart = Fuji{io,6} ;

        else
          for tiko = 1:size(Hash,2)
            if Ate(nano,tiko) ~=0
              eval(['A' num2str(tiko) ' = nchoosek(1:Hash(io,' num2str(tiko) '),Ate(nano,' num2str(tiko) ')) ;'])
            else
              eval(['A' num2str(tiko) ' = [] ;'])
            end
          end
          Stri = [] ;
          for tiko = 1:size(Hash,2)
            eval(['Check = A' num2str(tiko) ' ;'])
            if size(Check,1)>=1 && ~isempty(Check)
              Stri = [Stri , '1:' num2str(size(Check,1))] ;
            else
              Stri = [Stri , '0'] ;
            end
            if tiko == size(Hash,2)
            else
              Stri = [Stri , ','] ;
            end
          end
          eval(['Ansi = transpose(combvec_oct(' Stri ')) ;'])
          for kao = 1:size(Ansi,1)
            for tiko = 1:size(Hash,2)
              eval(['Check = A' num2str(tiko) ' ;'])
              if all(logical(Ansi(1,:)) == logical(Ate(nano,:))) && Ansi(kao,tiko) ~=0
                CC{kao,tiko} = Caper{io,tiko}(Check(Ansi(kao,tiko),:),1) ;
              end
            end
          end
          CC(:,1) = Double2Cell(CC(:,1)) ;
          Mindex = CheckCell(CC,Fuji{io,6}) ;
          CC(Mindex,:) = [] ;
          Spart = MergeCell(CC) ;
        end
        Holo = [Holo ; Spart] ;
      end
    end
    Shuji{io,1} = Holo ;

  end
end



##function [Caper , Hash] = XSpace(Fuji,MaximumCluster)
##
##  for io = 1:size(Fuji,1)
##    ko = MaximumCluster-Fuji{io,5} ;
##    Caper{io,1} = Fuji{io,6}' ;
##    for tiny = 2:ko
##      List = Fuji(find(cell2mat(Fuji(:,3)) == tiny),1).' ;
##      B = Fuji(io,6);
##      Indo = cell2mat(cellfun(@(x) all(ismember(x,B{1})),List(1,:),'UniformOutput',false)) ;
##      List = List.' ;
##      List(Indo==0,:) = [] ;
##      Caper{io,tiny} = List ;
##    end
##    Caper{io,end} = Fuji{io,6} ;
##  end
##  for io = 1:size(Caper,1)
##    Hash(io,:) = cell2mat(cellfun(@(x) size(x,1),Caper(io,:),'UniformOutput',false)) ;
##  end
##end

##function CC = Double2Cell(CC)
##  for emo = 1:size(CC,1)
##    CC{emo,1} = num2cell(CC{emo}) ;
##  end
##end

##function L = CheckCell(CC,Vector)
##  for pop = 1:size(CC,1)
##    Po = CC(pop,:) ;
##    Wow = [] ;
##    for kok = 1:length(Po)
##      Wow = [Wow , reshape(cell2mat(Po{kok}),1,[])] ;
##    end
##    if ~isempty(setdiff(Vector , unique(Wow)))
##      L(pop,1) = true ;
##    else
##      L(pop,1) = false ;
##    end
##  end
##end

##function Spart = MergeCell(CC)
##  Spart = [] ;
##  for pop = 1:size(CC,1)
##    Po = CC(pop,:) ;
##    Dow = [] ;
##    for kok = 1:length(Po)
##      Dow = [Dow ; Po{kok}] ; % Merge
##    end
##    Spart{pop,1} = Dow;
##  end
##end
