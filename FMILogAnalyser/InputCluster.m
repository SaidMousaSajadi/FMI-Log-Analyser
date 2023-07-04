function InputCluster()
  function GetData(obj,init = false)
    h = guidata(obj) ;
    switch (gcbo)
      case {h.OK}
        setappdata(0,'Clusters',get(h.Clusters,'string')) ;
        close(h.F)
    endswitch
  end

  R.F = figure("toolbar", "none",'uicontextmenu',[],'menubar','none','name',"Manual Selection of Fractures",'NumberTitle','off','resize','off','units','normalized',"Position", [0.3 0.3 0.3 0.3]) ;
  R.Text1 = uicontrol(R.F,'style','text','units','normalized','string','Clusters:','position',[0.01 0.80 0.25 0.085],'callback','','backgroundcolor',get(R.F,"Color")) ;
  R.Clusters = uicontrol(R.F,'style','edit','units','normalized','string','[1,3,6,9],[2,4,7],[5],[8]','position',[0.21 0.80 0.69 0.085],'callback','') ;
  R.OK = uicontrol(R.F,'style','pushbutton','units','normalized','string','Confirm Clusters','position',[0.35 0.1 0.35 0.085],'callback',@GetData) ;

  guidata (R.F, R) ;
  GetData(R.F,true)
end
