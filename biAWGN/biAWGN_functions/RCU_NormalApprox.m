function [epsilon] = RCU_NormalApprox(rho,rate,n)
   % INPUTS:
% rho   = Value of the snr in linear scale (no dB!)
% rate = Values of Rate at which the bound will be computed
% n     = Blocklength
%
% OUTPUT:
% epsilon     = Error probability

   Zmin = -20; Zmax = 20;
   fC = @(z) exp(-z.^2/2)/sqrt(2*pi) .* (log(2) - log(1+exp(-2*rho-2*sqrt(rho)*z)));
   C = integral(fC, Zmin, Zmax);
   fV = @(z) exp(-z.^2/2)/sqrt(2*pi) .* (log(2) - log(1+exp(-2*rho-2*sqrt(rho)*z))).^2;
   V = integral(fV, Zmin, Zmax) - (C)^2;
   for ii = 1 : length(rate)
      R = rate(ii) ; 
      epsilon(ii) = qfunc( (C-R) ./ (sqrt(V/n)) ) ; 
      % epsilon(ii) = qfunc( (C-R+ log(n)./(2.*n)) ./ (sqrt(V/n)) ) ; 
   end
end