function Kenzo = stringarray2cell(Anso)
  Kenzo = strsplit(Anso,'],[') ;
  if isrow(Kenzo)
    Kenzo = Kenzo.' ;
  else
  end
  %

  if size(Kenzo,1) == 1
    Kenzo{1,1} = str2num(Kenzo{1,1}) ;
  else
    for i = 1:size(Kenzo,1)
      if i == 1
        Kenzo{i,1} = [Kenzo{i,1} ']'] ;
        Kenzo{i,1} = str2num(Kenzo{i,1}) ;
      elseif i == size(Kenzo,1)
        Kenzo{i,1} = ['[' Kenzo{i,1}] ;
        Kenzo{i,1} = str2num(Kenzo{i,1}) ;
      else
        Kenzo{i,1} = ['[' Kenzo{i,1} ']'] ;
        Kenzo{i,1} = str2num(Kenzo{i,1}) ;
      end
    end
  end
end
