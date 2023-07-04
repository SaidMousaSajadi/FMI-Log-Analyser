function Clusters = MergeClusters(Fuji,Shuji,ISO)
  Clusters = [] ;
  for io = 1:size(Fuji,1) % Fuji
    BB = Shuji{io,1} ;
      for ki = 1:size(BB,1)
        DD = BB{ki,1} ;
        Prod = [DD ; isoadder(ISO) ; Fuji{io,1}] ;
        Clusters = [Clusters ; {Prod}] ;
      end
  end
end
function Pardon = isoadder(ISOs)
    Pardon = {} ;
    if isempty(ISOs)
    else
        for i = 1:length(ISOs)
            Pardon = [Pardon ; ISOs(i)] ;
        end
    end
end
