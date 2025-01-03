function [h_hat,theta_hat] = BlockPhase_ChEst(snr,np,theta)
   % INPUTS:
% snr     = Value of the snr in linear scale (no dB!)
% np      = number of pilot symbols per block
% theta   = phase rotation of a block
%
% OUTPUT:
% h_hat     = estimated phase rotation in polar coordinates
% theta_hat = estimated phase rotation (in radian scale)
   
   [h_real, h_imag]= pol2cart(theta,1); 
   h = h_real + 1i.*h_imag ; 
   x = sqrt(snr) .* ones(1,np); 
   y = x.*h + randcn(length(theta),np) ; 
   
%    [y_theta,y_rho] = cart2pol(real(y),imag(y)) ; 
%    theta_hat = mean(y_theta,2) ; 
%    h_hat = pol2cart(theta_hat,1)
   
   h_hat = y*x' ./ norm(x).^2 ; 
   h_hat = h_hat./ abs(h_hat); % Normalize so that the length will be always 1. 
   [theta_hat,~] = cart2pol(real(h_hat),imag(h_hat)) ; 
end