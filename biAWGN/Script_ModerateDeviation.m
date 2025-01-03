clc; clear all; 
close all; 


addpath('biAWGN_functions')

R = 1/2; 
if(R == 1/2)
   k_vec = [24:8:140] ; 
   % k_vec = 60 ; 
else
   k_vec = [18:6:104] ; 
end
n_vec = k_vec./R ; 
snrdB_vec_init = [-1,10] ; 
UER_target = 1e-5 ; 
FER_target = 1e-3 ; 
SNR_SEARCH = 1; 

for ii = 1 :length(k_vec)
   tic
   k = k_vec(ii) ; 
   n = k./R ;  
   snrdB_vec = snrdB_vec_init ;  
   while SNR_SEARCH == 1
      snrdB = mean(snrdB_vec) ; 
      snr = 10^(snrdB/10); 
      eps_FER = Erasure_ModerateDeviation(snr,R.*log(2),n,UER_target) ; 
      if(eps_FER < FER_target)
         snrdB_vec = [snrdB_vec(1),snrdB] ; 
      else
         snrdB_vec = [snrdB,snrdB_vec(2)] ;
      end
      if(snrdB_vec(2) - snrdB_vec(1) <= 0.04)
         eps_FER_save(ii) = eps_FER ; 
         SNR_save(ii) = snrdB ; 
         break ; 
      end
   end
end
SNR_out_lin = 10.^(SNR_save/10); 
SNR_save = log10(SNR_out_lin.*(1/R)./2)*10 
TextMat = [n_vec(:), SNR_save(:)]; 
text_name = ['ModerateDeviation_biAWGN_R05','.txt']; 
writematrix(TextMat,text_name,'Delimiter',' ') ;
