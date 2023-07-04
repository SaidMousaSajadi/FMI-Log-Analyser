function xData = MergeData(XData,IndicesMap,Vector)
% Example:
% XData = rand(100,2) ;
% IndicesMap = kmeans(XData,7) ;
% Vector = [1,3,5,6] ;
% 
% xData = MergeData(XData,IndicesMap,Vector) ;
    xData = [] ;
    for i = 1:length(Vector)
        xData = [xData ; XData(IndicesMap==Vector(i),:)] ;
    end
end