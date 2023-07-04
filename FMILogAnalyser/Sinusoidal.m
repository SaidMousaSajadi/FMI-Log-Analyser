function ycal = Sinusoidal(x,xdata)
  a = x(1) ;
  b = x(2) ;
  c = x(3) ;
  d = x(4) ;
  ycal = a*sin(b*(xdata+c)) + d ;
end
