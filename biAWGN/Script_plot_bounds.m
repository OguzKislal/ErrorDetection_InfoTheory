clc; clear all; 
close all; 

addpath('biAWGN_functions')

R = 1/2; % rate (in BPCU)
if(R == 1/2)
   k_vec = [24:8:125] ; % number of information bits
else
   k_vec = [16:6:104] ; 
end
snrdB_vec_init = [-1,10] ; 

METHOD_SEL = 3; % 1: Theorem 2 (CRC bound), 2: Forney's bound, 3: Theorem 3 (Threshold RCU bound), 4: Normal Approx. 

UER_target = 1e-5 ; 
FER_target = 1e-3 ; 
n_vec = round(k_vec./R); % blocklength n = k./R. 
N = round(20*(1./UER_target)); % Number of samples to be used in MC simulations (Only relevant for METHOD_SEL=3)

for rr = 1:length(n_vec) 
   tic
   n = n_vec(rr); 
   k = k_vec(rr) ; 
   snrdB_vec = snrdB_vec_init ;  
   % snrdB_vec = [-2,7]; 
   while 1
      snrdB = mean(snrdB_vec) ; 
      snr = 10^(snrdB/10);
      if(METHOD_SEL == 1)
         %% outer-code bound 
         disp('CRC Bit Sim Start')
         delta = 0:1:min(k_vec(rr)-1,20); % Number of CRC bits
         Rvec_delta = log(2) * (k + delta)/n; % rate after CRC bits are added. 
         eps_FER = rcu_saddle_biawgn(snr,Rvec_delta,n);
         eps_UER = 2.^(-delta).*eps_FER; 
         xAxis = delta ; 
      elseif(METHOD_SEL == 2)
         disp('Forney Sim Start')
         %% Forney's bound 
         T_vec = 0:0.02:1; % Forney's bound parameter
         epsFER_Forney = nan(length(T_vec),1);
         epsUER_Forney = nan(length(T_vec),1);
         for ii = 1:length(T_vec)
            [eps_FER(ii),eps_UER(ii)] = ForneyAchievError_biAWGN(n,R,snrdB,T_vec(ii));  
         end 
         xAxis = T_vec ;
      elseif(METHOD_SEL == 3)
         %% UER Thr. with RCU
         disp('UER RCU Sim Start')
         T_vec = linspace(-1.5,1.5,32) ; % Threshold (\lambda in the paper)
         T_vec = [-5:-2, T_vec, 2] ;  
         s_vec = [0.6 : 0.2 : 2.4] ; 
         % If you want a faster simulation uncomment UER_RCU_Saddlepoint_fixedTauHighData and comment UER_RCU_Saddlepoint_fixedTau
         % [eps_UER, eps_FER] = UER_RCU_Saddlepoint_fixedTauHighData(snr,n,R,N,T_vec) ;
         [eps_UER, eps_FER]  = UER_RCU_Saddlepoint_fixedTau(snr,n,R,N,T_vec,s_vec) ; 
         xAxis = T_vec ; 
      elseif(METHOD_SEL == 4)
         %% UER Thr. with RCU
         disp('CRC Gauss Approx. Start')
         delta = 0:1:min(k_vec(rr)-1,20); % Number of CRC bits
         k = n * R; 
         Rvec_delta = log(2) * (k + delta)/n;
         eps_FER = RCU_NormalApprox(snr,Rvec_delta,n); 
         eps_UER = 2.^(-delta).*eps_FER; 
         xAxis = delta ; 
      end
      
      if(eps_UER(end) > UER_target ) % Search for the trade-off parameter that gives UER_target
         xAxis_sel = xAxis(end) ; 
      elseif(eps_UER(1) < UER_target)
         xAxis_sel = xAxis(1); 
      else
         xAxis_sel = interp1(eps_UER,xAxis,UER_target,'linear') ; 
      end

      if(METHOD_SEL == 1 || METHOD_SEL == 4) % Number of CRC bits has to be an integer. 
      % (This condition can be overcome by using a randomized strategy, if needed this condition may be removed for a tighter bound!) 
         xAxis_sel = ceil(xAxis_sel) ; 
      end
      eps_FER_sel = interp1(xAxis,eps_FER,xAxis_sel,'linear','extrap') ; 

      % Decide if SNR should be increased or decreased
      if(eps_FER_sel < FER_target)
         snrdB_vec = [snrdB_vec(1),snrdB] ; 
      else
         snrdB_vec = [snrdB,snrdB_vec(2)] ; 
      end
      if(abs(diff(snrdB_vec)) < 0.02)
         SNR_out(rr) = snrdB_vec(2) ;
         break ; 
      end
   end
   toc
end

%SNR normalization according to the Eb/N0 definition in the paper.
SNR_out_lin = 10.^(SNR_out/10);  
SNR_out = log10(SNR_out_lin.*(1/R)./2)*10  ;

plot(n_vec,SNR_out) 
xlabel('n')
ylabel('Eb/N0')
grid on ; 
title_string = ['Method number:', num2str(METHOD_SEL)] ; 
title(title_string)

% Save the results
if(METHOD_SEL == 1)
   delta_SNR = SNR_out ; 
   save_script = ['biAWGN_DeltaBit_lowError_R',num2str(R),'.mat'] ; 
elseif(METHOD_SEL == 2)
   Forney_SNR = SNR_out ; 
   save_script = ['biAWGN_Forney_lowError_R',num2str(R),'.mat'] ; 
elseif(METHOD_SEL == 3)
   thr_RCU_SNR = SNR_out ; 
   save_script = ['biAWGN_Threshold_RCU_lowError_R',num2str(R),'.mat'] ; 
elseif(METHOD_SEL == 4)
   NormalApprox_SNR = SNR_out ; 
   save_script = ['biAWGN_NormalApprox_lowError',num2str(R),'.mat'] ; 
   TextMat = [n_vec(:), NormalApprox_SNR(:)]; 
   text_name = ['CRC_NormalApprox','.txt']; 
   writematrix(TextMat,text_name,'Delimiter',' ') ;
end
save(save_script) ; 
