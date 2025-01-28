clc; clear all; 
close all; 

addpath('biAWGN_functions')
rng(1)

k_vec = 64 ; 
n_vec = 128;
R = k_vec/n_vec; % rate (in BPCU)
snrdB_vec = [1.6 :0.4: 4] ; 
% snrdB_vec = 2.53
%SNR normalization according to the Eb/N0 definition in the paper.
snrdB_vec_lin = 10.^(snrdB_vec./10);  
snrdB_vec = log10(snrdB_vec_lin.*(1/R)./2).*10  ;

METHOD_SEL = 1; % 1: Theorem 2 (CRC bound), 2: Forney's bound, 3: Theorem 3 (Threshold RCU bound), 4: Normal Approx. 

UER_target = 1e-5 ; 
FER_target = 1e-3 ; 

N = 100 .* 1e6 ; 
for rr = 1:length(n_vec) 
   tic
   n = n_vec(rr); 
   k = k_vec(rr) ; 
   for SNR_counter = 1:length(snrdB_vec)
      snrdB = snrdB_vec(SNR_counter) ; 
      snr = 10^(snrdB/10);
      if(METHOD_SEL == 1)
         %% outer-code bound 
         disp('CRC Bit Sim Start')
         delta = 7
         Rvec_delta = log(2) * (k + delta)/n; % rate after CRC bits are added. 
         eps_FER(SNR_counter) = rcu_saddle_biawgn(snr,Rvec_delta,n);
         eps_UER(SNR_counter) = 2.^(-delta).*eps_FER(SNR_counter); 
      % elseif(METHOD_SEL == 2)
      %    disp('Forney Sim Start')
      %    %% Forney's bound 
      %    T_vec = 0:0.02:1; % Forney's bound parameter
      %    epsFER_Forney = nan(length(T_vec),1);
      %    epsUER_Forney = nan(length(T_vec),1);
      %    for ii = 1:length(T_vec)
      %       [eps_FER(ii),eps_UER(ii)] = ForneyAchievError_biAWGN(n,R,snrdB,T_vec(ii));  
      %    end 
      elseif(METHOD_SEL == 3)
         %% UER Thr. with RCU
         disp('UER RCU Sim Start')
         T_vec = -0.7450 ; 
         [eps_UER(SNR_counter), eps_FER(SNR_counter)] = UER_RCU_Saddlepoint_fixedTauHighData(snr,n,R,N,T_vec) ;
         % s_vec = [0.6 : 0.2 : 2.4] ; 
         % [eps_UER, eps_FER]  = UER_RCU_Saddlepoint_fixedTau(snr,n,R,N,T_vec,s_vec) ; 
      % elseif(METHOD_SEL == 4)
      %    %% UER Thr. with RCU
      %    disp('CRC Gauss Approx. Start')
      %    delta = 7 ; 
      %    k = n * R; 
      %    Rvec_delta = log(2) * (k + delta)/n;
      %    eps_FER = RCU_NormalApprox(snr,Rvec_delta,n); 
      %    eps_UER = 2.^(-delta).*eps_FER; 
      % end
   end
   toc
end



semilogy(snrdB_vec,eps_UER)
hold on ;
semilogy(snrdB_vec,eps_FER)
ylabel('UER and FER')
xlabel('Eb/N0')
ylim([1e-6, 1])
xlim([1.5, 4])
grid on ; 
title_string = ['Method number:', num2str(METHOD_SEL)] ; 
title(title_string)

% Save the results
if(METHOD_SEL == 1)
   save_script = ['biAWGN_DeltaBit_R',num2str(R),'_SNRvsErrProb.mat'] ; 
elseif(METHOD_SEL == 2)
   save_script = ['biAWGN_Forney_R',num2str(R),'_SNRvsErrProb.mat'] ; 
elseif(METHOD_SEL == 3)
   save_script = ['biAWGN_Threshold_RCU_R',num2str(R),'_SNRvsErrProb.mat'] ; 
elseif(METHOD_SEL == 4)
   NormalApprox_SNR = SNR_out ; 
   save_script = ['biAWGN_NormalApprox',num2str(R),'_SNRvsErrProb.mat'] ; 
   TextMat = [n_vec(:), NormalApprox_SNR(:)]; 
   text_name = ['CRC_NormalApprox','.txt']; 
   writematrix(TextMat,text_name,'Delimiter',' ') ;
end
save(save_script) ; 
