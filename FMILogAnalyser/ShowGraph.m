function ShowGraph(Ax,xs,ys,Kindex,G,ColorsInd)
  pkg load matgeom
  axes(Ax)

  % Clusters
  XS = [] ; YS = [] ;
  for i = 1:max(Kindex)
    % Colo = rand(1,3) ;
    % plot(Ax,xs(Kindex==i),ys(Kindex==i),'*','Color',Colo)
    XS = [XS mean(xs(Kindex==i))] ;
    YS = [YS mean(ys(Kindex==i))] ;
  end

  % Graph
  [r , c] = find(G~=0) ;
  [GN , GE] = drawGraph([XS' , YS'], [r , c]) ;

  % Graph Properties
  set(GN, 'markerfacecolor', ColorsInd(1), 'markeredgecolor', 'none','markersize',2) ; set(GE, 'color', ColorsInd(1),'linewidth',0.1);

  % Text
  for bsxi = 1:size(XS,2)
    text(XS(bsxi),YS(bsxi)+10,num2str(bsxi),"fontsize",8,"fontweight",'bold','color',ColorsInd(2),"backgroundcolor",ColorsInd(3),"edgecolor",ColorsInd(4))
  end

end
