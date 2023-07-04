function CC = Double2Cell(CC)
  for emo = 1:size(CC,1)
    CC{emo,1} = num2cell(CC{emo}) ;
  end
end
