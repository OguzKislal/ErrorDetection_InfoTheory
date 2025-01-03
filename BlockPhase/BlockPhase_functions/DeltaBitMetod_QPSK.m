function [epsilon] = DeltaBitMetod_QPSK(snr,h,h_hat,nc,R,N)
% INPUTS:
% snr   = Value of the snr in linear scale (no dB!)
% h     = Channel gain
% h_hat = Estimated channel gain
% nc    = number of data symbols per block
% R     = Values of Rate at which the bound will be computed
% N     = Number of samples used for CGF evaluation
%
% OUTPUT:
% epsilon     = Error probability obtained as a result of the computation of the bound


   QPSK_const = sqrt(snr).*[1+0i, 0+1i, -1+0i, 0-1i ].' ; 
   qLog_decod_func = @(x,y,h_hat) - abs(y-h_hat.*x).^2; 
   fi = @(y,h_hat,i) qLog_decod_func(QPSK_const(i),y,h_hat) ;
   exp_fiTau = @(y,h_hat,tau,i) exp(qLog_decod_func(QPSK_const(i),y,h_hat).*tau) ; 
   infDens_MGF = @(y,h_hat,tau) 1/4.*(exp_fiTau(y,h_hat,tau,1) + exp_fiTau(y,h_hat,tau,2) + ...
                  exp_fiTau(y,h_hat,tau,3) +  exp_fiTau(y,h_hat,tau,4)) ;
   infDens_CGF = @(y,h_hat,tau) log(infDens_MGF(y,h_hat,tau)) ; 
   infDens_CGF_div1_nom = @(y,h_hat,tau) fi(y,h_hat,1) .* exp_fiTau(y,h_hat,tau,1) + fi(y,h_hat,2) .* exp_fiTau(y,h_hat,tau,2)...
                                + fi(y,h_hat,3) .* exp_fiTau(y,h_hat,tau,3) + fi(y,h_hat,4) .* exp_fiTau(y,h_hat,tau,4); 
   infDens_CGF_div1_denom = @(y,h_hat,tau) exp_fiTau(y,h_hat,tau,1) + exp_fiTau(y,h_hat,tau,2) + exp_fiTau(y,h_hat,tau,3) + exp_fiTau(y,h_hat,tau,4) ;
   infDens_CGF_div1 = @(y,h_hat,tau) infDens_CGF_div1_nom(y,h_hat,tau) ./ infDens_CGF_div1_denom(y,h_hat,tau) ; 
   
   infDens_CGF_div2_term1_nom =  @(y,h_hat,tau) fi(y,h_hat,1).^2 .* exp_fiTau(y,h_hat,tau,1) + fi(y,h_hat,2).^2 .* exp_fiTau(y,h_hat,tau,2)...
                                + fi(y,h_hat,3).^2 .* exp_fiTau(y,h_hat,tau,3) + fi(y,h_hat,4).^2 .* exp_fiTau(y,h_hat,tau,4);
   infDens_CGF_div2_term1_denom =  @(y,h_hat,tau) infDens_CGF_div1_denom(y,h_hat,tau) ; 
   infDens_CGF_div2_term1 = @(y,h_hat,tau) infDens_CGF_div2_term1_nom(y,h_hat,tau) ./ infDens_CGF_div2_term1_denom(y,h_hat,tau) ; 
   infDens_CGF_div2_term2 = @(y,h_hat,tau) infDens_CGF_div1(y,h_hat,tau).^2 ; 
   infDens_CGF_div2 = @(y,h_hat,tau) infDens_CGF_div2_term1(y,h_hat,tau) - infDens_CGF_div2_term2(y,h_hat,tau) ; 
   
   for data_cnt = 1 : N 
      x = QPSK_const(randi([1,4],nc,1)) ; 
      y = h.*x + randcn(nc,1)   ;  
      tau_vec = [-2, 20]  ; 
      delta(data_cnt) = sum(qLog_decod_func(x,y,h_hat),1) ; 
      for tau_cnt = 1 : 20 
         tauSel = mean(tau_vec) ; 
         CGF_div1 = sum(infDens_CGF_div1(y,h_hat,tauSel),1) ; 
         diff_temp = (CGF_div1 - delta(data_cnt));  
         if(abs(diff_temp) <= 1e-3)
            break ; 
         elseif(CGF_div1 < delta(data_cnt))
            tau_vec = [tauSel, tau_vec(2)] ; 
         else
            tau_vec = [tau_vec(1), tauSel] ; 
         end
      end
      CGF_D0 = sum(infDens_CGF(y,h_hat,tauSel),1) ; 
      CGF_D1 = sum(infDens_CGF_div1(y,h_hat,tauSel),1) ; 
      CGF_D2 = sum(infDens_CGF_div2(y,h_hat,tauSel),1) ; 
      tauShift_term = 0 ; 
      if(CGF_D2 <0) % Rarely, CGF_D2 can be ~ -1e-16 which is basically 0
%          CGF_D2 
%          disp('CGF_D2 is minus, auto-correcting!')
         CGF_D2 = 0 ; 
      end
      if(tauSel > 0)
         SP_approx(1,data_cnt) = exp(CGF_D0-tauSel.*CGF_D1).*exp((tauSel.^2/2).*CGF_D2)...
                                                   .*qfunc(tauSel.*sqrt(CGF_D2)) ; 
      else
         SP_approx(1,data_cnt) = 1-exp(CGF_D0-tauSel.*CGF_D1).*exp((tauSel.^2/2).*CGF_D2).*qfunc(-tauSel.*sqrt(CGF_D2)) ; 
      end
      
   end
   epsilon = mean(min(1,(2.^(nc.*R) -1).*SP_approx),2) ; 
end