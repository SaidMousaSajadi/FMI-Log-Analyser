close all ; clear all ; clc ;

graphics_toolkit qt % or fltk(PLC/old) , qt(apps) , gnuplot(sites)
pkg load image % images
pkg load optim
pkg load io

##pkg load statistics % pdist
##pkg load ltfat % normalize

warning off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Processing(STR,obj)
  switch lower(STR)
    case 'e'
      set(obj,'string',' Error Occurred.')
      set(obj,'backgroundcolor',[1 0.785 0.785])
    case 'p'
      set(obj,'string',' In Processing!')
      set(obj,'backgroundcolor',[1 0.863 0.666])
    case 'r'
      set(obj,'string',' Ready')
      set(obj,'backgroundcolor',[0.785 1 0.785])
    otherwise
      set(obj,'string',STR)
      set(obj,'backgroundcolor',[0.785 1 0.785])
  end
end


function h = EdgeReduction(h)
  BW = h.BW ;
  if h.FlagTransform == 1 % radon

    Theta = -90:0.5:90;
    [R_Space , Off_Set_From_Centre] = radon(BW,Theta);
    R_Sort = sort(unique(R_Space),"descend") ; R_Sort(R_Sort<1) = [] ;
    [row_peak,col_peak] = find(ismember(R_Space,R_Sort(1:round(numel(R_Sort)/2))));
    Off_Set_From_Centre_Peak = Off_Set_From_Centre(row_peak);
    Theta_Peak = Theta(col_peak);
    centerX = ceil(size(BW,2)/2);
    centerY = ceil(size(BW,1)/2);
    R_Coord = [] ;
    for i = 1:length(Theta_Peak)
      [x1,y1] = pol2cart(Theta_Peak(i),Off_Set_From_Centre_Peak(i));
      try
        if BW(round(centerY+y1),round(centerX-x1)) == 1
          P_Coord = [centerX-x1 , centerY+y1] ;
          R_Coord = [R_Coord ; P_Coord] ;
          % plot([centerX-x1],[centerY+y1],"*b",'LineWidth',2) ;
        elseif BW(round(centerY-y1),round(centerX+x1)) == 1
          P_Coord = [centerX+x1 , centerY-y1] ;
          R_Coord = [R_Coord ; P_Coord] ;
          % plot([centerX+x1],[centerY-y1],"*b",'LineWidth',2) ;
        end
      catch
      end
    end
    R_Coord = round(R_Coord) ;
    L = bwlabel(BW) ;
    L0 = zeros(size(L)) ;
    for i = 1:max(L(:))
      Li = L ;
      Li(Li~=i) = 0 ;
      [c , r] = ind2sub(size(Li),find(Li==i)) ;
      if ~isempty(intersect([r,c],R_Coord,'rows'))
        L0 = L0 + Li ;
      end
    end
    BW = logical(L0) ;

  else % hough

    Theta = -90:0.5:90;
    [H_Space , ~ , Off_Set_From_Centre] = hough(BW,'Theta',Theta);
    H_Sort = sort(unique(H_Space),"descend") ; H_Sort(H_Sort<0.1) = [] ;
    [row_peak,col_peak] = find(ismember(H_Space,H_Sort(1:round(numel(H_Sort)/2))));
    Off_Set_From_Centre_Peak = Off_Set_From_Centre(row_peak);
    Theta_Peak = Theta(col_peak);
    centerX = ceil(size(BW,2)/2);
    centerY = ceil(size(BW,1)/2);
    H_Coord = [] ;
    for i = 1:length(Theta_Peak)
      [x1,y1] = pol2cart(Theta_Peak(i),Off_Set_From_Centre_Peak(i));
      try
        if BW(round(centerY+y1),round(centerX-x1)) == 1
          P_Coord = [centerX-x1 , centerY+y1] ;
          H_Coord = [H_Coord ; P_Coord] ;
          % plot([centerX-x1],[centerY+y1],"*b",'LineWidth',2) ;
        elseif BW(round(centerY-y1),round(centerX+x1)) == 1
          P_Coord = [centerX+x1 , centerY-y1] ;
          H_Coord = [H_Coord ; P_Coord] ;
          % plot([centerX+x1],[centerY-y1],"*b",'LineWidth',2) ;
        end
      catch
      end
    end
    % Additional Points
    Peaks = houghpeaks(H_Space,round(numel(H_Sort)/2));
    Lines = houghlines(BW,Theta,Off_Set_From_Centre,Peaks,'FillGap',2,'MinLength',2);
    for k = 1:length(Lines)
      P_Coord = [Lines(k).point1 ; Lines(k).point2];
      H_Coord = [H_Coord ; P_Coord] ;
    end
    H_Coord = round(H_Coord) ;
    L = bwlabel(BW) ;
    L0 = zeros(size(L)) ;
    for i = 1:max(L(:))
      Li = L ;
      Li(Li~=i) = 0 ;
      [c , r] = ind2sub(size(Li),find(Li==i)) ;
      if ~isempty(intersect([r,c],H_Coord,'rows'))
        L0 = L0 + Li ;
      end
    end
    BW = logical(L0) ;

  end
  if h.FlagInvert
    h.ONeal = imshow(~BW,'parent',h.Ax2) ;
  else
    h.ONeal = imshow(BW,'parent',h.Ax2) ;
  end
  axis(h.Ax2,'on','image');
  set(h.Ax2,'linewidth',2) ;
  set(h.Ax2,'xcolor',[0.49 0.99 0]) ;
  set(h.Ax2,'ycolor',[0.49 0.99 0]) ;
  set(h.Ax2,'ytick',[]) ;
  set(h.Ax2,'xtick',[]) ;
  h.BW = BW ;
  guidata(gcf,h) % update handles


end

function Set_Fracture_Label(Ax,IM,Cluster,X,Kindex,cycle,Quality,Method,FlagFracture,FlagLabel,FracturesColor,FlagColorScale)
  if FlagColorScale
    IM_G = double(rgb2gray(IM))/255 ; % GrayScale
    imshow(IM_G,'Parent',Ax) ;
  else
    imshow(IM,'Parent',Ax) ;
  end
  hold(Ax,'on') ;
  if FlagFracture
    for j = 1:length(Cluster)
      x = MergeData(X,Kindex,Cluster{j,:}) ;
      wt = [20*ones(size(x,1),1)] ;
      [~ , Coeffs] = fitfunction(x(:,1),x(:,2),wt,cycle,Quality,Method) ;
      C = [1:1:cycle] ;
      switch lower(Method)
        case 'lsq'
          R = Sinusoidal(Coeffs,C) ;
        case 'nlin'
          B = 2*pi/cycle ;
          Sinusoidal = @(b, x) (b(1)*sin(B*(x+b(2))) + b(3));
          R = Sinusoidal(Coeffs,C) ;
      end
      plot(Ax,C,R,'Color',FracturesColor,'LineWidth',2.5)
    end
  end
  if FlagLabel
    for j = 1:length(Cluster)
      x = MergeData(X,Kindex,Cluster{j,:}) ;
      wt = [20*ones(size(x,1),1)] ;
      [~ , Coeffs] = fitfunction(x(:,1),x(:,2),wt,cycle,Quality,Method) ;
      C = [1:1:cycle] ;
      switch lower(Method)
        case 'lsq'
          R = Sinusoidal(Coeffs,C) ;
        case 'nlin'
          B = 2*pi/cycle ;
          Sinusoidal = @(b, x) (b(1)*sin(B*(x+b(2))) + b(3));
          R = Sinusoidal(Coeffs,C) ;
      end
      h = max(R) - min(R) ;
      d = cycle - 1 ;
      Dip(j) = rad2deg(pi/2) - atand (h/d) ;
      plot(Ax,C,R,'Color',FracturesColor,'LineWidth',0.01)
      text(cycle/2,R(round(cycle/2)),num2str(Dip(j),'%2.2f'),"fontsize",10,"fontweight",'bold',"backgroundcolor",'w',"edgecolor",'red')
    end
  end
  hold(Ax,'off') ;
end


function h = ProcessOnBase(h)
  if h.FlagColorScale
    imshow(h.IM_G,'parent',h.Ax1) ;
  else
    imshow(h.IM,'parent',h.Ax1) ;
  end
  axis(h.Ax1,'on','image');
  set(h.Ax1,'linewidth',2) ;
  set(h.Ax1,'xcolor',[0.49 0.99 0]) ;
  set(h.Ax1,'ycolor',[0.49 0.99 0]) ;
  set(h.Ax1,'ytick',[]) ;
  set(h.Ax1,'xtick',[]) ;
end

function h = ProcessOnEdge(h)
  Crude = double(rgb2gray(h.IM))/255 ;
  % Maybe Append: 1)CNN , 2)Wavelet Decomposition
  if h.FlagDenoise == 1
    f = fspecial("Gaussian", [1 , 1] , 2);
    Grayo = imfilter(Crude , f , "replicate");
  elseif h.FlagDenoise == 2
    Grayo = medfilt2(Crude);
  elseif h.FlagDenoise == 3
    Grayo = wiener2(Crude,[2 2]);
  else
    Grayo = Crude ;
  end

  Sigma = get(h.SigmaBar,'value')*100 ;
  Thr = [get(h.LowBar,'value') , get(h.HighBar,'value')] ;
  BW = edge(Grayo,'Canny',Thr,Sigma) ;
  h.BW = BW ; 
  h.FlagCluster = 0 ;
  if h.FlagInvert
    h.ONeal = imshow(~BW,'parent',h.Ax2) ;
  else
    h.ONeal = imshow(BW,'parent',h.Ax2) ;
  end
  
  axis(h.Ax2,'on','image');
  set(h.Ax2,'linewidth',2) ;
  set(h.Ax2,'xcolor',[0.49 0.99 0]) ;
  set(h.Ax2,'ycolor',[0.49 0.99 0]) ;
  set(h.Ax2,'ytick',[]) ;
  set(h.Ax2,'xtick',[]) ;
  guidata(gcf,h) % update handles

end

function h = ProcessOnExport(h)
  if h.FlagColorScale
    imshow(h.IM_G,'parent',h.Ax3) ;
  else
    imshow(h.IM,'parent',h.Ax3) ;
  end
  cycle = size(h.BW,2) ;
  if h.FlagCluster & (get(h.FracturesCheck,'value') | get(h.LabelsCheck,'value'))
    Set_Fracture_Label(h.Ax3,h.IM,h.Cluster,h.X,h.Kindex,cycle,h.Quality,h.Method,get(h.FracturesCheck,'value'),get(h.LabelsCheck,'value'),h.FracturesColor,h.FlagColorScale)
  end
  axis(h.Ax3,'on','image');
  set(h.Ax3,'linewidth',2) ;
  set(h.Ax3,'xcolor',[0.49 0.99 0]) ;
  set(h.Ax3,'ycolor',[0.49 0.99 0]) ;
  set(h.Ax3,'ytick',[]) ;
  set(h.Ax3,'xtick',[]) ;
end

function h = ProcessOn3D(h)
  if h.FlagCaliper
    [X,Y,Z] = CylinderWith(h.IM_G_Interp,'Caliper',h.CaliperData) ;
  else
    [X,Y,Z] = CylinderWith(h.IM_G_Interp,'',[]) ;
  end
  surf(h.Ax4,X,Y,Z,h.IM_G_Interp,'edgecolor','none','FaceColor','texturemap') ; %
  hold(h.Ax4,'on') ;
  colormap(h.Ax4,h.Map)
  ##colormap(Ax3D,'gray')
  
  patch(h.Ax4,X(1,:),Y(1,:),Z(1,:),'k','facecolor',[0.145 0.145 0.145],'edgecolor','none')
  patch(h.Ax4,X(end,:),Y(end,:),Z(end,:),'k','facecolor',[0.145 0.145 0.145],'edgecolor','none')
  

  if h.FlagCluster & get(h.FracturesCheck,'value')
    ShowFractures3D(h)
  end

  set(h.Ax4,'zdir','reverse')
  view(h.Ax4,[35 35])
  axis(h.Ax4,[-1.5, 1.5, -1.5, 1.5, 0, 1])
  set(h.Ax4,'DataAspectRatio',[1 1 1/8])
  set(h.Ax4,'color','none')
  set(h.Ax4,'ytick',[]) ;
  set(h.Ax4,'xtick',[]) ;
  set(h.Ax4,'ztick',[]) ;
  set(h.Ax4,'ycolor','none') ;
  set(h.Ax4,'xcolor','none') ;
  set(h.Ax4,'zcolor','none') ;
  hold(h.Ax4,'off') ;
end

function [Gra,X,Kindex,ONeal] = GraphSimulator(BW,Ax,cycle,Quality,Method,GraphColor,FlagInvert)
  %% With Machine Vision (Clustering)
  %% My Method
  
  [ys xs] = ind2sub(size(BW),find(BW==1)) ; % Coordinate of Edges
  X = [xs ys] ; % Merge Coordinate
  Kindex = dbscan_oct(X,10,10); % Clusters data with density-based algorithm
  [min(Kindex) max(Kindex)] ; % 0 is noisey data

  Gra = zeros(max(Kindex),max(Kindex)) ; % Graph Matrix
  
  % Selection of couples of clusters, Connecting edges in pairs
  for i = 1:max(Kindex)
    for j = i+1:max(Kindex)
      x = MergeData(X,Kindex,[i,j]) ;
      wt = [20*ones(size(xs(Kindex==i))) ; 1*ones(size(xs(Kindex==j)))] ;
      % To Show
  ##    imshow(BW,'Parent',Ax2) ;
  ##    hold on
  ##    plot(xt,yt,'*','Color',rand(1,3))
      % Save in Matrix Gra
      Gra(j,i) = fitfunction(x(:,1),x(:,2),wt,cycle,Quality,Method) ;
    end
  end

  % Make top triangle matrix
  Gra = Gra+Gra.'-diag([diag(Gra)]) ;

  %% Show Graph
  if FlagInvert
    ONeal = imshow(~BW,'parent',Ax) ;
  else
    ONeal = imshow(BW,'parent',Ax) ;
  end
  hold(Ax,'on') ;
  try
    ShowGraph(Ax,xs,ys,Kindex,Gra,GraphColor)
  catch
  end
  axis(Ax,'on','image');
  set(Ax,'linewidth',2) ;
  set(Ax,'xcolor',[0.49 0.99 0]) ;
  set(Ax,'ycolor',[0.49 0.99 0]) ;
  set(Ax,'ytick',[]) ;
  set(Ax,'xtick',[]) ;
  hold(Ax,'off') ;
end

function R = GeneInterp(Gray)
  [x, y] = meshgrid(1:size(Gray,2), 1:size(Gray,1));
  [X,Y] = meshgrid(linspace(1,size(Gray,2),round(size(Gray,2)/3)+1) , linspace(1,size(Gray,1),round(size(Gray,1)/3))); % interpolate /3
  R = interp2(x,y,Gray,X,Y,'cubic');
end

function [X,Y,Z] = CylinderWith(Gray,Method,CaliperData)
  switch lower(Method)
    case {'caliper'}
      Depth = [CaliperData(1,1),CaliperData(end,1)] ;
      CaliperData(1:end,1) = linspace(1,0,size(CaliperData,1)) ;
      xi = linspace(1,0,size(Gray,1)) ;
      yi = pchip(CaliperData(:,1),CaliperData(:,2),xi) ;
      [X, Y, Z] = cylinder(yi./mean(yi) , size(Gray,2)-1);
    otherwise
      Depth = [] ;
      [X, Y, Z] = cylinder(ones(size(Gray,1)) , size(Gray,2)-1);
  end
end

function [xx,yy,zz] = GeneFracture(Gray,XYData)
  thetafunc = @(x) (359/size(Gray,2))*(x-1) ;
  theta = thetafunc(XYData(:,1)) ; %deg, will be change to rad
  [Xsc Ysc Zsc] = pol2cart(deg2rad(theta),ones(size(theta)),XYData(:,2)) ; %change polar to cartesian
  modelfun = @(b, x) b(1) + b(2).*x(1,:) + b(3).*x(2,:) ;
  Wsc = ones(size(Xsc')) ;
  beta0 = [0; 0; 0];
  [beta_w] = nlinfit([Xsc' ; Ysc'], (Zsc./size(Gray,1)).', modelfun, beta0, [], "weights", Wsc) ;
  [xx yy] = meshgrid([-1.5:0.1:1.5] , [-1.5:0.1:1.5]) ;
  zz = beta_w(1) + beta_w(2).*xx + beta_w(3).*yy ;
end

function ShowFractures3D(h)
  hold(h.Ax4,'on') ;
  Cluster = h.Cluster ;
  for j = 1:length(Cluster)
    x = MergeData(h.X,h.Kindex,Cluster{j,:}) ;
    cycle = size(h.IM_G,2) ; % cycle to fitting Sinusoidal function
    wt = [1*ones(size(x,1),1)] ;
    [~ , Coeffs] = fitfunction(x(:,1),x(:,2),wt,cycle,h.Quality,h.Method) ;
    C = [1:1:cycle] ;
    switch lower(h.Method)
      case 'lsq'
        R = Sinusoidal(Coeffs,C) ;
      case 'nlin'
        B = 2*pi/cycle ;
        Sinusoidal = @(b, x) (b(1)*sin(B*(x+b(2))) + b(3));
        R = Sinusoidal(Coeffs,C) ;
    end
    [xx,yy,zz] = GeneFracture(h.IM_G_Interp,[C.',R.']./3) ; % interpolate /3
    surf(h.Ax4,xx,yy,zz,'edgecolor','none','facecolor',h.FracturesColor,'facealpha',0.75)
  end
  hold(h.Ax4,'off') ;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

root.Quality = 0.8 ; % Quality Control
root.Method = 'nlin' ; % Method for Regression
root.FracturesColor = 'blue' ;
root.GraphColor = "ykwr" ;



root.FlagIM = 0 ;
root.FlagDenoise = 0 ;
root.FlagTransform = 0 ;
root.FlagGraph = 0 ;
root.FlagSimu = 0 ;
root.FlagCluster = 0 ;
root.Flag3D = 0 ;
root.FlagCaliper = 0 ;
root.FlagInvert = 0 ; 
root.FlagColorScale = 0 ;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Update_UI(obj,init = false)
  h = guidata(obj) ; % get handles
  switch (gcbo)

  case {h.Fig}
    % get(h.Fig,'Position')
  case {h.Open}
    Processing('p',h.Process)
    [h.FileName, h.FilePath, h.FileIndex] = uigetfile({"*.jpeg;*.jpg;*.tiff;*.tif;*.png", "Supported Picture Formats";"*.png", "Portable Network Graphics";"*.jpeg ; *.jpg", "Joint Photographic Experts Group";"*.tif ; *.tiff", "Tagged Image File Format"},"Import A Dynamic FMI Image") ; 
    if (h.FileIndex) ~= 0
      NameSpl = strsplit(h.FileName,".") ;
      switch lower(NameSpl{1,end})
        case {'png','jpg','jpeg','tiff','tif'} % Better Work imfinfo
          [h.IM , ~]= imread([h.FilePath h.FileName]) ;
          h.Gray = colormap(gray(16)) ;
          [~ , h.CMap] = rgb2ind(h.IM) ;
          if h.FlagColorScale
            h.Map = h.Gray ;
          else
            h.Map = h.CMap ;
          end
          h.IM_G = double(rgb2gray(h.IM))/255 ; % GrayScale
          h.IM_G_Interp = GeneInterp(h.IM_G) ; % Interpolation
          h = ProcessOnBase(h) ;
          h = ProcessOnEdge(h) ;
          h = ProcessOnExport(h) ;
          linkaxes([h.Ax1 , h.Ax2 , h.Ax3]) ;
          h.Axis = axis(h.Ax1) ;
          h.FlagIM = 1 ;
          h.FlagCaliper = 0 ;
          h.FlagGraph = 0 ;
          h.FlagSimu = 0 ;
          h.FlagCluster = 0 ;
          if h.Flag3D
            h = ProcessOn3D(h) ;
          end
          guidata(gcf,h) % update handles
      end
    end
    Processing('r',h.Process)

  case {h.Home}
    Processing('p',h.Process)
    pan off
    rotate3d off
    axis(h.Ax1,h.Axis) ; axis(h.Ax2,h.Axis) ; axis(h.Ax3,h.Axis)
    view(h.Ax1,[0 90]) ; view(h.Ax2,[0 90]) ; view(h.Ax3,[0 90])
    linkaxes([h.Ax1 , h.Ax2 , h.Ax3]) ;
    Processing('r',h.Process)

  case {h.ZoomOn}
    Processing('p',h.Process)
    zoom on
    Processing('r',h.Process)

  case {h.ZoomOff}
    Processing('p',h.Process)
    zoom off
    Processing('r',h.Process)

  case {h.ZoomIn}
    Processing('p',h.Process)
    zoom(1.001)
    zoom off
    Processing('r',h.Process)

  case {h.ZoomOut}
    Processing('p',h.Process)
    zoom(0.999)
    zoom off
    Processing('r',h.Process)

  case {h.Pan}
    Processing('p',h.Process)
    pan on
    Processing('r',h.Process)

  case {h.ArrowTop_T}
    Processing('p',h.Process)
    axis(h.Ax1,axis(h.Ax1) + [0 0 -1 0]) ;
    Processing('r',h.Process)

  case {h.ArrowTop_D}
    Processing('p',h.Process)
    axis(h.Ax1,axis(h.Ax1) + [0 0 +1 0]) ;
    Processing('r',h.Process)

  case {h.ArrowDown_T}
    Processing('p',h.Process)
    axis(h.Ax1,axis(h.Ax1) + [0 0 0 -1]) ;
    Processing('r',h.Process)

  case {h.ArrowDown_D}
    Processing('p',h.Process)
    axis(h.Ax1,axis(h.Ax1) + [0 0 0 +1]) ;
    Processing('r',h.Process)

  case {h.ArrowLeft_R}
    Processing('p',h.Process)
    axis(h.Ax1,axis(h.Ax1) + [+1 0 0 0]) ;
    Processing('r',h.Process)

  case {h.ArrowLeft_L}
    Processing('p',h.Process)
    axis(h.Ax1,axis(h.Ax1) + [-1 0 0 0]) ;
    Processing('r',h.Process)

  case {h.ArrowRight_R}
    Processing('p',h.Process)
    axis(h.Ax1,axis(h.Ax1) + [0 +1 0 0]) ;
    Processing('r',h.Process)

  case {h.ArrowRight_L}
    Processing('p',h.Process)
    axis(h.Ax1,axis(h.Ax1) + [0 -1 0 0]) ;
    Processing('r',h.Process)

  case {h.RFC}
    pan off
    try
      Processing('p',h.Process)
      Dim = axis(h.Ax1) ;
      h.IM = h.IM(ceil(Dim(3))+1:fix(Dim(4))-1,ceil(Dim(1))+1:fix(Dim(2))-1,:) ;
      h.IM_G = double(rgb2gray(h.IM))/255 ; % GrayScale
      h.IM_G_Interp = GeneInterp(h.IM_G) ; % Interpolation
      [~ , h.CMap] = rgb2ind(h.IM) ;
      guidata(gcf,h) % update handles
      h = ProcessOnBase(h) ;
      h = ProcessOnEdge(h) ;
      h = ProcessOnExport(h) ;
      linkaxes([h.Ax1 , h.Ax2 , h.Ax3]) ;
      h.Axis = axis(h.Ax1) ;
      axis(h.Axis) ;
      if h.Flag3D
        h = ProcessOn3D(h) ;
      end
      guidata(gcf,h) % update handles
      Processing('r',h.Process)
    catch
      Processing('e',h.Process)
    end

  case {h.Radio0}
    Processing('p',h.Process)
    set (h.Radio0, "value", 1);
    set (h.Radio1, "value", 0);
    set (h.Radio2, "value", 0);
    set (h.Radio3, "value", 0);
    h.FlagDenoise = 0 ;
    h.FlagGraph = 0 ;
    if h.FlagIM
      h = ProcessOnEdge(h) ;
      h = ProcessOnExport(h) ;
      if h.Flag3D
        h = ProcessOn3D(h) ;
      end
      colormap(h.Ax4,h.Map)
    end
    guidata(gcf,h) % update handles
    Processing('r',h.Process)

  case {h.Radio1}
    Processing('p',h.Process)
    set (h.Radio0, "value", 0);
    set (h.Radio1, "value", 1);
    set (h.Radio2, "value", 0);
    set (h.Radio3, "value", 0);
    h.FlagDenoise = 1 ;
    h.FlagGraph = 0 ;
    if h.FlagIM
      h = ProcessOnEdge(h) ;
      h = ProcessOnExport(h) ;
      if h.Flag3D
        h = ProcessOn3D(h) ;
      end
      colormap(h.Ax4,h.Map)
    end
    guidata(gcf,h) % update handles
    Processing('r',h.Process)

  case {h.Radio2}
    Processing('p',h.Process)
    set (h.Radio0, "value", 0);
    set (h.Radio1, "value", 0);
    set (h.Radio2, "value", 1);
    set (h.Radio3, "value", 0);
    h.FlagDenoise = 2 ;
    h.FlagGraph = 0 ;
    if h.FlagIM
      h = ProcessOnEdge(h) ;
      h = ProcessOnExport(h) ;
      if h.Flag3D
        h = ProcessOn3D(h) ;
      end
      colormap(h.Ax4,h.Map)
    end
    guidata(gcf,h) % update handles
    Processing('r',h.Process)

  case {h.Radio3}
    Processing('p',h.Process)
    set (h.Radio0, "value", 0);
    set (h.Radio1, "value", 0);
    set (h.Radio2, "value", 0);
    set (h.Radio3, "value", 1);
    h.FlagDenoise = 3 ;
    h.FlagGraph = 0 ;
    if h.FlagIM
      h = ProcessOnEdge(h) ;
      h = ProcessOnExport(h) ;
      if h.Flag3D
        h = ProcessOn3D(h) ;
      end
      colormap(h.Ax4,h.Map)
    end
    guidata(gcf,h) % update handles
    Processing('r',h.Process)

  case {h.SigmaBar}
    Processing('p',h.Process)
    set(h.SigmaText,'string',['Sigma: ' num2str(get(h.SigmaBar,'value')*100,'%3.2f')]) ;
    h.FlagGraph = 0 ;
    h.FlagCluster = 0 ;
    if h.FlagIM
      h = ProcessOnEdge(h) ;
      h = ProcessOnExport(h) ;
      if h.Flag3D
        h = ProcessOn3D(h) ;
      end
      colormap(h.Ax4,h.Map)
    end
    guidata(gcf,h) % update handles
    Processing('r',h.Process)

  case {h.LowBar}
    Processing('p',h.Process)
    if get(h.HighBar,'value') <= get(h.LowBar,'value')
      set(h.LowBar,'value',get(h.HighBar,'value')-0.01) ;
      set(h.HighBar,'value',get(h.LowBar,'value')+0.01) ;
    end
    set(h.LowText,'string',['Low Thr.: ' num2str(get(h.LowBar,'value'),'%.2f')]) ;
    h.FlagGraph = 0 ;
    h.FlagCluster = 0 ;
    if h.FlagIM
      h = ProcessOnEdge(h) ;
      h = ProcessOnExport(h) ;
      if h.Flag3D
        h = ProcessOn3D(h) ;
      end
      colormap(h.Ax4,h.Map)
    end
    guidata(gcf,h) % update handles
    Processing('r',h.Process)

  case {h.HighBar}
    Processing('p',h.Process)
    if get(h.HighBar,'value') <= get(h.LowBar,'value')
      set(h.HighBar,'value',get(h.LowBar,'value')+0.01) ;
      set(h.LowBar,'value',get(h.HighBar,'value')-0.01) ;
    end
    set(h.HighText,'string',['High Thr.: ' num2str(get(h.HighBar,'value'),'%.2f')]) ;
    h.FlagGraph = 0 ;
    h.FlagCluster = 0 ;
    if h.FlagIM
      h = ProcessOnEdge(h) ;
      h = ProcessOnExport(h) ;
      if h.Flag3D
        h = ProcessOn3D(h) ;
      end
      colormap(h.Ax4,h.Map)
    end
    guidata(gcf,h) % update handles
    Processing('r',h.Process)
    
  case {h.RadioHough}
    Processing('p',h.Process)
    set (h.RadioHough, "value", 1);
    set (h.RadioRadon, "value", 0);
    h.FlagTransform = 0 ;
    h.FlagGraph = 0 ;
    h.FlagCluster = 0 ;
    guidata(gcf,h) % update handles
    Processing('r',h.Process)

  case {h.RadioRadon}
    Processing('p',h.Process)
    set (h.RadioHough, "value", 0);
    set (h.RadioRadon, "value", 1);
    h.FlagTransform = 1 ;
    h.FlagGraph = 0 ;
    h.FlagCluster = 0 ;
    guidata(gcf,h) % update handles
    Processing('r',h.Process)

  case {h.EdgeButton}
    Processing('p',h.Process)
    h.FlagGraph = 0 ;
    h.FlagCluster = 0 ;
    if h.FlagIM
      h = ProcessOnEdge(h) ;
      h = ProcessOnExport(h) ;
      if h.Flag3D
        h = ProcessOn3D(h) ;
      end
      colormap(h.Ax4,h.Map)
    end
    guidata(gcf,h) % update handles
    Processing('r',h.Process)
  
  case {h.EdgeReductionButton}
    Processing('p',h.Process)
    h.FlagGraph = 0 ;
    h.FlagCluster = 0 ;
    if h.FlagIM
      h = EdgeReduction(h) ;
      h = ProcessOnExport(h) ;
      if h.Flag3D
        h = ProcessOn3D(h) ;
      end
      colormap(h.Ax4,h.Map)
    end
    guidata(gcf,h) % update handles
    Processing('r',h.Process)

  case {h.GraphMap}
    Processing('p',h.Process)
    Quality = h.Quality ; % Fitting Quality
    Method = h.Method ; % Fitting Method
    h.FlagSimu = 0 ;
    h.FlagCluster = 0 ;
    if h.FlagIM
      cycle = size(h.BW,2) ; % cycle to fitting Sinusoidal function
      [h.Gra,h.X,h.Kindex,h.ONeal] = GraphSimulator(h.BW,h.Ax2,cycle,Quality,Method,h.GraphColor,h.FlagInvert) ;
      h = ProcessOnExport(h) ;
      if h.Flag3D
        h = ProcessOn3D(h) ;
      end
      h.FlagGraph = 1 ;
    end
    guidata(gcf,h) % update handles
    Processing('r',h.Process)
    
  case {h.Simulation}
    Processing('p',h.Process)
    Quality = h.Quality ; % Fitting Quality
    Method = h.Method ; % Fitting Method
    h.FlagCluster = 0 ;
    if h.FlagGraph & h.FlagIM
      cycle = size(h.BW,2) ; % cycle to fitting Sinusoidal function
      h.Clusters = WhatACoolYourPhase(h.X,h.Gra,h.Kindex,cycle,Quality,Method) ;
      h = ProcessOnExport(h) ;
      if h.Flag3D
        h = ProcessOn3D(h) ;
      end
      h.FlagSimu = 1 ;
      Processing('r',h.Process)
    elseif ~h.FlagGraph & h.FlagIM
      cycle = size(h.BW,2) ; % cycle to fitting Sinusoidal function
      [h.Gra,h.X,h.Kindex,h.ONeal] = GraphSimulator(h.BW,h.Ax2,cycle,Quality,Method,h.GraphColor,h.FlagInvert) ;
      h.Clusters = WhatACoolYourPhase(h.X,h.Gra,h.Kindex,cycle,Quality,Method) ;
      h = ProcessOnExport(h) ;
      if h.Flag3D
        h = ProcessOn3D(h) ;
      end
      h.FlagSimu = 1 ;
      Processing('r',h.Process)

    end
    guidata(gcf,h) % update handles

  case {h.SaveListOfClusters}
    Processing('p',h.Process)
    if h.FlagSimu
      [FileName, FilePath, FileIndex] = uiputfile({"*.json", "JavaScript Object Notation"},"Save Clusters","Clusters.json") ;
      if (FileIndex) ~= 0
        Cell = Cluster2Cell(h.Clusters) ;
        SturctToJSON(Cell,[FilePath ,  FileName])
      end
    end
    Processing('r',h.Process)

  case {h.ManualCluster}
    Processing('p',h.Process)
    if h.FlagGraph
      InputCluster ;
      uiwait ;
      Middelar = getappdata(0,'Clusters') ; % Get String of Clusters from User
      if ~isempty(Middelar)
        h.Cluster = stringarray2cell(Middelar) ;
        h.FlagCluster = 1 ;
        guidata(gcf,h) % update handles
      end
      h = ProcessOnExport(h) ;
      if h.Flag3D
        h = ProcessOn3D(h) ;
      end
    end
    guidata(gcf,h) % update handles
    Processing('r',h.Process)
    
    
  case {h.AICluster}
    Processing('p',h.Process)
    if h.FlagSimu
      disp('It does not work in this version.')
      Processing("not this version",h.Process)
      h.FlagCluster = 1 ;
    end
    guidata(gcf,h) % update handles
    Processing('r',h.Process)
    
  case {h.View3D}
    Processing('p',h.Process)
    if h.Flag3D
      set(h.P4,'visible','off')
      set(h.Ax4,'visible','off')
      cla(h.Ax4)
      set(h.Home3D,'visible','off')
      set(h.RoteOn,'visible','off')
      set(h.RoteOff,'visible','off')
      set(h.ImportCaliper,'visible','off')
      h.Flag3D = 0 ;
    else
      set(h.P4,'visible','on')
      set(h.Ax4,'visible','on')
      set(h.Home3D,'visible','on')
      set(h.RoteOn,'visible','on')
      set(h.RoteOff,'visible','on')
      set(h.ImportCaliper,'visible','on')
      h.Flag3D = 1 ;
      guidata(gcf,h) % update handles
      if h.FlagIM
        h = ProcessOn3D(h) ;
      end
    end
    guidata(gcf,h) % update handles
    Processing('r',h.Process)
  
  case {h.Home3D}
    rotate3d off
    view(h.Ax4,[35 35])
    axis(h.Ax4,[-1.5, 1.5, -1.5, 1.5, 0, 1])

  case {h.RoteOn}
    rotate3d on

  case {h.RoteOff}
    rotate3d off

  case {h.ImportCaliper}
    Processing('p',h.Process)
    if h.FlagIM
      [h.FileName, h.FilePath, h.FileIndex] = uigetfile({"*.csv", "Comma-Separated Values"},"Import A Caliper Data") ; 
      if (h.FileIndex) ~= 0
        NameSpl = strsplit(h.FileName,".") ;
        switch lower(NameSpl{1,end})
          case {'csv'}
            h.CaliperData = dlmread([h.FilePath h.FileName],'',0,0:1) ;
        end
        h.FlagCaliper = 1 ;
        h = ProcessOn3D(h) ;
      end
    end
    guidata(gcf,h) % update handles
    Processing('r',h.Process)

  case {h.FracturesCheck}
    Processing('p',h.Process)
    if h.FlagCluster 
      h = ProcessOnExport(h) ;
      if h.Flag3D
        h = ProcessOn3D(h) ;
      end
    end
    
    guidata(gcf,h) % update handles
    Processing('r',h.Process)

  case {h.LabelsCheck}
    Processing('p',h.Process)
    if h.FlagCluster 
      h = ProcessOnExport(h) ;
    end
    guidata(gcf,h) % update handles
    Processing('r',h.Process)
	
  case {h.SaveAx1}
    if h.FlagIM
      Processing('p',h.Process) 
      [FileName, FilePath, FileIndex] = uiputfile({"*.jpeg;*.jpg;*.tiff;*.tif;*.png", "Supported Picture Formats";"*.png", "Portable Network Graphics";"*.jpeg ; *.jpg", "Joint Photographic Experts Group";"*.tif ; *.tiff", "Tagged Image File Format"},"Save Original","Original") ;
      if (FileIndex) ~= 0
        Fr = getframe(h.Ax1);
        IMAG = frame2im(Fr);
        imwrite(IMAG, [FilePath ,  FileName])
      end
      Processing('r',h.Process)
    end

  case {h.SaveAx2}
    if h.FlagIM
      Processing('p',h.Process) 
      [FileName, FilePath, FileIndex] = uiputfile({"*.jpeg;*.jpg;*.tiff;*.tif;*.png", "Supported Picture Formats";"*.png", "Portable Network Graphics";"*.jpeg ; *.jpg", "Joint Photographic Experts Group";"*.tif ; *.tiff", "Tagged Image File Format"},"Save Edges and Graph","Edges") ;
      if (FileIndex) ~= 0
        Fr = getframe(h.Ax2);
        IMAG = frame2im(Fr);
        imwrite(IMAG, [FilePath ,  FileName])
      end
      Processing('r',h.Process)
    end
  
  case {h.SaveAx3}
    if h.FlagIM
      Processing('p',h.Process) 
      [FileName, FilePath, FileIndex] = uiputfile({"*.jpeg;*.jpg;*.tiff;*.tif;*.png", "Supported Picture Formats";"*.png", "Portable Network Graphics";"*.jpeg ; *.jpg", "Joint Photographic Experts Group";"*.tif ; *.tiff", "Tagged Image File Format"},"Save Fractures","Fractures") ;
      if (FileIndex) ~= 0
        Fr = getframe(h.Ax3);
        IMAG = frame2im(Fr);
        imwrite(IMAG, [FilePath ,  FileName])
      end
      Processing('r',h.Process)
    end
  
  case {h.SaveAx4}
    if h.FlagIM & h.Flag3D
      Processing('p',h.Process) 
      [FileName, FilePath, FileIndex] = uiputfile({"*.jpeg;*.jpg;*.tiff;*.tif;*.png", "Supported Picture Formats";"*.png", "Portable Network Graphics";"*.jpeg ; *.jpg", "Joint Photographic Experts Group";"*.tif ; *.tiff", "Tagged Image File Format"},"Save 3D Model","3D") ;
      if (FileIndex) ~= 0
        Fr = getframe(h.Ax4);
        IMAG = frame2im(Fr);
        imwrite(IMAG, [FilePath ,  FileName])
      end
      Processing('r',h.Process)
    end


  case {h.ColorScale}
    Processing('p',h.Process)
    if h.FlagColorScale
      h.FlagColorScale = 0 ;
      set(h.ColorScale,'label','Color Scale To Gray') ;
      if h.FlagIM
        h.Map = h.CMap ;
      end
    else
      h.FlagColorScale = 1 ;
      set(h.ColorScale,'label','Color Scale To Hot') ;
      if h.FlagIM
        h.Map = h.Gray ;
      end
    end
    if h.FlagIM
      h = ProcessOnBase(h) ;
      h = ProcessOnExport(h) ;
      if h.Flag3D
        h = ProcessOn3D(h) ;
      end
      
    end
    guidata(gcf,h) % update handles
    Processing('r',h.Process)
  
  case {h.EdgeColorScale}
    Processing('p',h.Process)
    if h.FlagInvert
      h.FlagInvert = 0 ;
    else
      h.FlagInvert = 1 ;
    end
    guidata(gcf,h) % update handles
    if h.FlagIM
      h.ONeal = imshow(~get(h.ONeal,'cdata'),'parent',h.Ax2) ;
      axis(h.Ax2,'on','image');
      set(h.Ax2,'linewidth',2) ;
      set(h.Ax2,'xcolor',[0.49 0.99 0]) ;
      set(h.Ax2,'ycolor',[0.49 0.99 0]) ;
      set(h.Ax2,'ytick',[]) ;
      set(h.Ax2,'xtick',[]) ;
      colormap(h.Ax4,h.Map)
      if h.FlagGraph
        Quality = h.Quality ; % Fitting Quality
        Method = h.Method ; % Fitting Method
        cycle = size(h.BW,2) ; % cycle to fitting Sinusoidal function
        [h.Gra,h.X,h.Kindex,h.ONeal] = GraphSimulator(h.BW,h.Ax2,cycle,Quality,Method,h.GraphColor,h.FlagInvert) ;
      end
    end

    guidata(gcf,h) % update handles
    Processing('r',h.Process)
	
  
  case {h.GraphColorScale}
    Processing('p',h.Process)
    switch h.GraphColor
      case {'ykwr'}
        h.GraphColor = 'bwkr' ;
      case {'bwkr'}
        h.GraphColor = 'ykwr' ;
    end
    guidata(gcf,h) % update handles
    if h.FlagIM
      if h.FlagGraph
        Quality = h.Quality ; % Fitting Quality
        Method = h.Method ; % Fitting Method
        cycle = size(h.BW,2) ; % cycle to fitting Sinusoidal function
        [h.Gra,h.X,h.Kindex,h.ONeal] = GraphSimulator(h.BW,h.Ax2,cycle,Quality,Method,h.GraphColor,h.FlagInvert) ;
      end
    end
    guidata(gcf,h) % update handles
    Processing('r',h.Process)



  case {h.FracturehColorScale}
    Processing('p',h.Process)
    switch lower(h.FracturesColor)
      case {'blue'}
        h.FracturesColor = 'red' ;
        set(h.FracturehColorScale,'label','Change Fractures Color to Blue') ;
      case {'red'}
        h.FracturesColor = 'blue' ;
        set(h.FracturehColorScale,'label','Change Fractures Color to Red') ;
    end
    guidata(gcf,h) % update handles
    if h.FlagCluster ;
      h = ProcessOnExport(h) ;
      if h.Flag3D
        h = ProcessOn3D(h) ;
      end
    end
    guidata(gcf,h) % update handles
    Processing('r',h.Process)
    
	
  end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

root.Fig = figure("toolbar", "none",'uicontextmenu',[],'menubar','none','name',"FMI Log Analyser",'NumberTitle','off','units','normalized',"Position", [2.1962e-03 6.2500e-02 9.9414e-01 8.6979e-01],"CloseRequestFcn",'exit','resizefcn',@Update_UI) ;

root.F = uimenu("label", "&File", "accelerator", "f");
root.E = uimenu("label", "&Edit", "accelerator", "e");
root.V = uimenu("label", "&View", "accelerator", "v");
root.P = uimenu("label", "&Processing", "accelerator", "p");
root.H = uimenu("label", "&Help", "accelerator", "h");

# Subs
root.Open = uimenu(root.F, "label", "&Open Image", "accelerator", "O", "callback", @Update_UI);
root.Exit = uimenu(root.F, "label", "E&xit", "accelerator", "X","callback", 'close(root.Fig) ; exit');


root.SaveAx1 = uimenu(root.E, "label", "Save Original Image", "accelerator", "", "callback", @Update_UI);
root.SaveAx2 = uimenu(root.E, "label", "Save Edges Image", "accelerator", "", "callback", @Update_UI);
root.SaveAx3 = uimenu(root.E, "label", "Save Fractures Image", "accelerator", "", "callback", @Update_UI);
root.SaveAx4 = uimenu(root.E, "label", "Save 3D Image", "accelerator", "", "callback", @Update_UI);


root.ColorScale = uimenu(root.V, "label", "Color Scale To Gray", "accelerator", "", "callback", @Update_UI);
root.EdgeColorScale = uimenu(root.V, "label", "Invert Edges Color Scale", "accelerator", "", "callback", @Update_UI);
root.GraphColorScale = uimenu(root.V, "label", "Change Graphs Color", "accelerator", "", "callback", @Update_UI);
root.FracturehColorScale = uimenu(root.V, "label", "Change Fractures Color to Red", "accelerator", "", "callback", @Update_UI);
root.View3D = uimenu(root.V, "label", "3D View", "accelerator", "", "callback", @Update_UI);

root.EdgeButton = uimenu(root.P, "label", "Edge &Detection(Canny Automated)", "accelerator", "D", "callback", @Update_UI);
root.EdgeReductionButton = uimenu(root.P, "label", "Edge &Reduction(Hough/Radon)", "accelerator", "R", "callback", @Update_UI);
root.EdgeSimulation = uimenu(root.P, "label", "Fractures Simulation", "accelerator", "");
root.GraphMap = uimenu(root.EdgeSimulation, "label", "Edges Connection Graph", "accelerator", "", "callback", @Update_UI);
root.Simulation = uimenu(root.EdgeSimulation, "label", "Simulating Fractures", "accelerator", "", "callback", @Update_UI);
root.SaveListOfClusters = uimenu(root.P, "label", "Save List of Fractures(JSON File)", "accelerator", "", "callback", @Update_UI);
root.ManualCluster = uimenu(root.P, "label", "Manual Selection of Fractures", "accelerator", "", "callback", @Update_UI);
root.AICluster = uimenu(root.P, "label", "Fractures Selection with AI(Not work in this version)", "accelerator", "", "callback", @Update_UI);





# In Window
% Panel 1
root.P1 = uipanel(root.Fig,'units','normalized','Position',[0.01 0.72 0.199 0.27],'visible','on','backgroundcolor',get(root.Fig,'Color'),'title','Original Options:');

root.Home = uicontrol(root.P1,"string", "🏠",'units','normalized','Position',[0.085 0.8 0.15 0.15],'fontsize',12,'callback',@Update_UI); %
root.ZoomOn = uicontrol(root.P1,"string", "🔎",'units','normalized','Position',[0.245 0.8 0.15 0.15],'fontsize',12,'callback',@Update_UI); %
root.ZoomOff = uicontrol(root.P1,"string", "🔎⛔",'units','normalized','Position',[0.405 0.8 0.15 0.15],'fontsize',12,'callback',@Update_UI); %
root.ZoomIn = uicontrol(root.P1,"string", "🔎+",'units','normalized','Position',[0.565 0.8 0.15 0.15],'fontsize',12,'callback',@Update_UI); %
root.ZoomOut = uicontrol(root.P1,"string", "🔎-",'units','normalized','Position',[0.725 0.8 0.15 0.15],'fontsize',12,'callback',@Update_UI); %
root.Pan = uicontrol(root.P1,"string", "✋",'units','normalized','Position',[0.085 0.335 0.15 0.15],'fontsize',12,'callback',@Update_UI); %
root.ArrowTop_T = uicontrol(root.P1,"string", "⏫",'units','normalized','Position',[0.535 0.635 0.1 0.15],'fontsize',12,'callback',@Update_UI); %
root.ArrowTop_D = uicontrol(root.P1,"string", "⏬",'units','normalized','Position',[0.535 0.485 0.1 0.15],'fontsize',12,'callback',@Update_UI); %
root.ArrowDown_T = uicontrol(root.P1,"string", "⏫",'units','normalized','Position',[0.535 0.185 0.1 0.15],'fontsize',12,'callback',@Update_UI); %
root.ArrowDown_D = uicontrol(root.P1,"string", "⏬",'units','normalized','Position',[0.535 0.035 0.1 0.15],'fontsize',12,'callback',@Update_UI); %
root.ArrowLeft_R = uicontrol(root.P1,"string", "⏩",'units','normalized','Position',[0.435 0.335 0.1 0.15],'fontsize',12,'callback',@Update_UI); %
root.ArrowLeft_L = uicontrol(root.P1,"string", "⏪",'units','normalized','Position',[0.335 0.335 0.1 0.15],'fontsize',12,'callback',@Update_UI); %
root.ArrowRight_L = uicontrol(root.P1,"string", "⏪",'units','normalized','Position',[0.635 0.335 0.1 0.15],'fontsize',12,'callback',@Update_UI); %
root.ArrowRight_R = uicontrol(root.P1,"string", "⏩",'units','normalized','Position',[0.735 0.335 0.1 0.15],'fontsize',12,'callback',@Update_UI); %
root.RFC = uicontrol(root.P1,"string", "📌",'units','normalized','Position',[0.085 0.035 0.15 0.15],'fontsize',12,'callback',@Update_UI); %
# ⬅️ 🔙 ⬆️ 🔝 ➔ 🔜 ⬇️  👓 💻 🔂 🔁 🍳 ☣︎ 📌

% Axes 1
root.Ax1 = axes(root.Fig,'units','normalized',"Position", [0.011 0.0325 0.197 0.68],'box','on','xtick',[],'ytick',[]);

% Panel 2
root.P2 = uipanel(root.Fig,'units','normalized','Position',[0.27 0.72 0.199 0.27],'visible','on','backgroundcolor',get(root.Fig,'Color'),'title','Edges Options:');
root.TextDenoise = uicontrol(root.P2,'style','text','units','normalized','position',[0.035 0.85 0.30 0.1],'string','Denoising:','backgroundcolor',get(root.Fig,'Color'),'horizontalalignment','left','fontsize',8);
root.Radio0 = uicontrol (root.P2, "style", "radiobutton", "string","No","units","normalized","position", [0.235 0.84 0.317 0.1],'fontsize',7,'backgroundcolor',get(root.Fig,"Color"),'value' , 1,'callback',@Update_UI);
root.Radio1 = uicontrol (root.P2, "style", "radiobutton", "string","Gaussian","units","normalized","position", [0.385 0.84 0.317 0.1],'fontsize',7,'backgroundcolor',get(root.Fig,"Color"),'value' , 0,'callback',@Update_UI);
root.Radio2 = uicontrol (root.P2, "style", "radiobutton", "string","Median","units","normalized","position", [0.610 0.84 0.317 0.1],'fontsize',7,'backgroundcolor',get(root.Fig,"Color"),'value' , 0,'callback',@Update_UI);
root.Radio3 = uicontrol (root.P2, "style", "radiobutton", "string","Wiener","units","normalized","position", [0.810 0.84 0.317 0.1],'fontsize',7,'backgroundcolor',get(root.Fig,"Color"),'value' , 0,'callback',@Update_UI);
root.SigmaBar = uicontrol(root.P2,'style','slider','units','normalized','position',[0.315 0.72 0.67 0.1],'sliderstep',[0.0001 0.01],'min',0.0001,'max',1,'Value',0.02,'callback',@Update_UI);
root.SigmaText = uicontrol(root.P2,'style','text','units','normalized','position',[0.035 0.72 0.27 0.1],'string',['Sigma: ' num2str(get(root.SigmaBar,'Value')*100,'%3.2f')],'backgroundcolor',get(root.Fig,"Color"),'horizontalalignment','left','fontsize',8);
root.LowBar = uicontrol(root.P2,'style','slider','units','normalized','position',[0.315 0.59 0.67 0.1],'sliderstep',[0.01 0.1],'min',0,'max',1,'Value',0.3,'callback',@Update_UI);
root.LowText = uicontrol(root.P2,'style','text','units','normalized','position',[0.035 0.59 0.27 0.1],'string',['Low Thr.: ' num2str(get(root.LowBar,'Value'),'%.2f')],'backgroundcolor',get(root.Fig,"Color"),'horizontalalignment','left','fontsize',8);
root.HighBar = uicontrol(root.P2,'style','slider','units','normalized','position',[0.315 0.46 0.67 0.1],'sliderstep',[0.01 0.1],'min',0,'max',1,'Value',0.7,'callback',@Update_UI);
root.HighText = uicontrol(root.P2,'style','text','units','normalized','position',[0.035 0.46 0.27 0.1],'string',['High Thr.: ' num2str(get(root.HighBar,'Value'),'%.2f')],'backgroundcolor',get(root.Fig,"Color"),'horizontalalignment','left','fontsize',8);
root.TextTransform = uicontrol(root.P2,'style','text','units','normalized','position',[0.035 0.33 0.30 0.1],'string','Transform:','backgroundcolor',get(root.Fig,'Color'),'horizontalalignment','left','fontsize',8);
root.RadioHough = uicontrol (root.P2, "style", "radiobutton", "string","Hough","units","normalized","position", [0.245 0.32 0.317 0.1],'fontsize',7,'backgroundcolor',get(root.Fig,"Color"),'value' , 1,'callback',@Update_UI);
root.RadioRadon = uicontrol (root.P2, "style", "radiobutton", "string","Radon","units","normalized","position", [0.435 0.32 0.317 0.1],'fontsize',7,'backgroundcolor',get(root.Fig,"Color"),'value' , 0,'callback',@Update_UI);


% Axes 2
root.Ax2 = axes(root.Fig,'units','normalized',"position", [0.271 0.0325 0.197 0.68],'box','on','xtick',[],'ytick',[],'colormap',colormap('gray'));


% Panel 3
root.P3 = uipanel(root.Fig,'units','normalized','Position',[0.53 0.72 0.199 0.27],'visible','on','backgroundcolor',get(root.Fig,'Color'),'title','Fractures Options:');
root.FracturesCheck = uicontrol(root.P3,'style','checkbox','units','normalized','position',[0.035 0.85 0.39 0.1],'string','Show Fractures','backgroundcolor',get(root.Fig,'Color'),'fontsize',9,'value',0,'callback',@Update_UI);
root.LabelsCheck = uicontrol(root.P3,'style','checkbox','units','normalized','position',[0.035 0.72 0.39 0.1],'string','Show Labels','backgroundcolor',get(root.Fig,'Color'),'fontsize',9,'value',0,'callback',@Update_UI);

% Axes 3
root.Ax3 = axes(root.Fig,'units','normalized',"position", [0.531 0.0325 0.197 0.68],'box','on','xtick',[],'ytick',[]);

% Panel 4
root.P4 = uipanel(root.Fig,'units','normalized','Position',[0.79 0.72 0.199 0.27],'visible','off','backgroundcolor',get(root.Fig,'Color'),'title','3D Model:');
root.Home3D = uicontrol(root.P4,"string", "🏠",'units','normalized','Position',[0.085 0.8 0.15 0.15],'visible','off','fontsize',12,'callback',@Update_UI); %
root.RoteOn = uicontrol(root.P4,"string", "🔁",'units','normalized','Position',[0.245 0.8 0.15 0.15],'visible','off','fontsize',12,'callback',@Update_UI); %
root.RoteOff = uicontrol(root.P4,"string", "🔁⛔",'units','normalized','Position',[0.405 0.8 0.15 0.15],'visible','off','fontsize',12,'callback',@Update_UI); %
root.ImportCaliper = uicontrol(root.P4,"string", "Import Caliper",'units','normalized','Position',[0.085 0.635 0.31 0.15],'visible','off','fontsize',8,'callback',@Update_UI); %

% Axes 4
root.Ax4 = axes(root.Fig,'units','normalized',"position", [0.791 0.0325 0.197 0.68],'visible','off','box','on','xtick',[],'ytick',[]);


uimenu(root.H, "label", "&Documentation", "accelerator", "","callback", "system(['start ./Manual.pdf']) ;"); %
uimenu(root.H, "label", "About Me", "accelerator", "","callback", "web('https://www.linkedin.com/in/seyed-mousa-sajadi-8284b1124/','-new')"); %


root.Process = uicontrol(root.Fig,'style','text','units','normalized','position',[0.0 0.0 0.05 0.025],'string',' Ready','backgroundcolor',[0.785 1 0.785],"horizontalalignment",'left','fontsize',7);


guidata(gcf,root) ;
Update_UI(gcf,true)
pause

% cd .\Octave\FMILogApplication
% cls ; octave Main.m
