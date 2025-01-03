function [pt] = Erasure_ModerateDeviation(rho,R,n,pu)
   % INPUTS:
% rho   = Value of the snr in linear scale (no dB!)
% R     = Values of Rate at which the bound will be computed
% n     = Blocklength
% pu    = Targeted undetected error probability
%
% OUTPUT:
% pt     = Total error probability

   Zmin = -20; Zmax = 20;
   fC = @(z) exp(-z.^2/2)/sqrt(2*pi) .* (log(2) - log(1+exp(-2*rho-2*sqrt(rho)*z)));
   C = integral(fC, Zmin, Zmax);
   fV = @(z) exp(-z.^2/2)/sqrt(2*pi) .* (log(2) - log(1+exp(-2*rho-2*sqrt(rho)*z))).^2;
   V = integral(fV, Zmin, Zmax) - (C)^2;
   t_vec = linspace(0.01,0.499,100) ; 

   logM = R.*n;   
   for t_looper = 1 : length(t_vec)
      t = t_vec(t_looper); 
      a = (logM - n.*C)./(n.^(1-t)) ; 
      b_vec = linspace(0.001,a,100) ; 
      for b_looper = 1: length(b_vec)
         b = b_vec(b_looper) ; 
         pu_temp(b_looper) = exp(-n.^(1-t).*b) ; 
      end
      b_sel = interp1(pu_temp,b_vec,pu,'linear','extrap') ; 
      pt_save(t_looper) = exp(-n.^(1-2.*t).*(a-b_sel).^2 ./ (2.*V)) ; 
   end  
   pt = min(pt_save) ; 
   

end