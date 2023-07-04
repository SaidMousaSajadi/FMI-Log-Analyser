function [Gof_R,Coeffs] = fitfunction(xt,yt,wt,cycle,Quality,Method)
  % Need Sinusoidal
  B = 2*pi/cycle ;
  C = [1:1:cycle] ;
  switch lower(Method)
    case 'lsq'
      options = optimset('Algorithm','lsqcurvefit','Display','off','MaxIter',10000,'TolFun',1e-10,'FinDiffRelStep',1e-10) ; % 'MaxFunctionEvaluations',5000,'FiniteDifferenceStepSize',1e-16,'OptimalityTolerance',1e-16,'StepTolerance',1e-16 ,'OptimalityTolerance',1e-16,'MaxFunctionEvaluations',100000
      Guess = [0 B 0 0] ;
      Coeffs = lsqcurvefit(@Sinusoidal,Guess,xt,yt,[-inf B -inf -inf],[inf B inf inf],options) ;
    case 'nlin'
      Sinusoidal = @(b, x) (b(1)*sin(B*(x+b(2))) + b(3));
      Coeff = [0; 0; 0];
      [Coeffs] = nlinfit(xt, yt, Sinusoidal, Coeff, [], "weights", wt) ;
  end

  % Efficiency
  R = Sinusoidal(Coeffs,C) ;
  RL = Sinusoidal(Coeffs,xt)-20 ; % Down Gap
  RH = Sinusoidal(Coeffs,xt)+20; % Top Gap
  Effi = nnz(yt<RH & yt>RL)/length(yt) ;

  % R-Squared
  Y = Sinusoidal(Coeffs,xt) ;
  rsquare = 1-(sum((yt - Y).^2) ./ sum((yt - mean(yt)).^2)) ;

  Qual= Quality ; %
  Efficiency = Quality ; %

  if Effi >= Efficiency & rsquare >= Qual
    Gof_R = (rsquare+Effi)/2 ;
  else
    Gof_R = 0 ;
  end
end
