function Clustersi = sortcellelements(Clustersi)
    for tech = 1:size(Clustersi,1)
        Q = Clustersi{tech,1} ;
        Q(:,2) = cellfun(@length,Q,'UniformOutput',false) ;
        Q = sortrows(Q,2) ;
        Clustersi{tech,1} = Q(:,1) ;
    end
end