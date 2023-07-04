function Full = uniquecells(Clustersi)
    Indexi = cell2mat(cellfun(@length,Clustersi,'UniformOutput',false)) ;
    Index = unique(Indexi) ;
    Full = [] ;
    for jo = 1:numel(Index)
        FastClus = Clustersi(Indexi==Index(jo),:) ;
        Dup = [] ;
        for qo = 1:size(FastClus,1)
            if qo == 1
                Cup = FastClus(qo,1) ;
                Dup = [Dup ; Cup] ;
            else
                Cup = FastClus(qo,1) ;
                L = Similar(Dup,FastClus{qo,1}) ;
                if any(L)
                else
                    Dup = [Dup ; Cup] ;
                end
            end
        end
        Full = [Full ; Dup] ;
    end
end
%

function L = Similar(Dup,Cup)
    for i = 1:size(Dup,1)
        New = Dup{i,1} ;
        for j = 1:size(Cup,1)
            for k = 1:size(New,1)
                ll(j,k) = isequal(New{k,:},Cup{j,:}) ;
            end % New
        end % Cup
        L(i,1) = all(any(ll,2)) ;
    end % Dup
end 


