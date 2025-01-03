clc; clear all; 
close all; 

addpath('BlockPhase_functions')

%% Simulation inits.
R_main = 1 ; % Rate (BPCU)
s = 1 ; % s optimization may be skipeed since its impact is negligable for this case
k = 120 ; % Number of information bits in a block
n = k./R_main; % Blocklength

k_sel = k  ; 
np = 10 ; % Number of pilot symbols in a block
nc = n ; 

R = (n.*R_main)./nc ; % Effective rate (pilot symbols do not carry any information)

UER_target = 1e-5 ; 
FER_target = 1e-3 ; 

N = 1e5; 
theta_N = 40; 

METHOD_SEL = 2 ; % 1:Delta-Bit, 2: Thr. RCU

% CRC bit search for METHOD_SEL = 1
% This search over number of CRC bits is not needed for this type of a simulation. 
delta_vec = [0,2,4,6,7,8,10,12,14] ; 
% Threshold init. for METHOD_SEL = 2
T_vec = linspace(-1,2,32); 
T_vec = [linspace(-10,-2,5), T_vec, linspace(3,10,5)] ; 


if(METHOD_SEL == 1)
   eps_FixedTh_UER = nan(length(delta_vec),theta_N);
   eps_FixedTh_FER = nan(length(delta_vec),theta_N);
else
   eps_FixedTh_UER = nan(length(T_vec),theta_N);
   eps_FixedTh_FER = nan(length(T_vec),theta_N);
end
% SNR interval that the SNR search will be conducted
snrdB_vec = [-7,8]; 

%% SNR Search
while 1
   tic
   %% Channel estimation
   snrdB = mean(snrdB_vec) ; 
   snr = 10^(snrdB/10);
% To simplify the simulations, we first evaluate error probabilities for
% different rotations -> theta - \hat{theta}. Then sample theta and
% \hat{theta} and then find the average error probability averaged over rotation. 

   theta = 0/180*pi;   %Phase noise 
   [h_real, h_imag]= pol2cart(theta,1); % Polar coordinates for the phase noise
   h_vec = h_real + 1i.*h_imag ; 
   h_vec = repmat(h_vec,theta_N,1) ; 
   thetaHat_main_vec = linspace(0,60,theta_N)./180*pi ; % Phase rotation
   [h_hat_real, h_hat_imag]= pol2cart(thetaHat_main_vec,1);  
   h_hat_vec = h_hat_real + 1i.*h_hat_imag  ; 
   for th_cnt = 1 : theta_N % For each rotation in 0-60 degrees 
      rng(th_cnt)
      eps_tot_DeltaRCU_temp = nan(1,length(delta_vec)) ; 
      eps_UER_DeltaRCU_temp = nan(1,length(delta_vec)); 
      tic
      h = h_vec(th_cnt) ; 
      h_hat = h_hat_vec(th_cnt); 
      if(METHOD_SEL == 1)
         % CRC outer-code bound (Theorem 2)
         disp('Delta Bit Sim Start')
         k = nc * R; 
         R_delta = (k + delta_vec)/nc;
         for ii = 1 : length(R_delta) 
            eps_tot_DeltaRCU_temp(ii) = DeltaBitMetod_QPSK(snr,h,h_hat,nc,R_delta(ii),N) ;
            eps_UER_DeltaRCU_temp(ii) = 2.^(-delta_vec(ii)).*eps_tot_DeltaRCU_temp(ii); 
         end
         eps_FixedTh_FER(:,th_cnt) = eps_tot_DeltaRCU_temp ;
         eps_FixedTh_UER(:,th_cnt) = eps_UER_DeltaRCU_temp ; 
      elseif(METHOD_SEL == 2)
         % UER Thr. with RCU (Theorem 3)
         disp('UER RCU Sim Start')
         [eps_FixedTh_UER(:,th_cnt), eps_FixedTh_FER(:,th_cnt)] = ThrMetod_QPSK(snr,h,h_hat,nc,R,T_vec,N,s) ; 
      end
   end
   %% Average over phase rotation (theta- \hat{theta})
   theta = unifrnd(-pi,pi,1e5,1) ; % phase is uniformly selected
   [~,thetaHat_sim_vec]= BlockPhase_ChEst(snr,np,theta) ;  % Channel estimation
   thetaDiff = abs(theta-thetaHat_sim_vec) ; % Rotation
   thetaDiff_temp = abs(thetaDiff - 2.*pi) ; % Rotation normalized in 2*pi step 1
   thetaDiff = min([thetaDiff,thetaDiff_temp],[],2) ; % Rotation normalized in 2*pi step 2
   
   eps_FixedTh_UER = eps_FixedTh_UER + sort(1e-12.*rand(size(eps_FixedTh_UER)),'descend') ; 
   eps_FixedTh_FER = eps_FixedTh_FER + sort(1e-12.*rand(size(eps_FixedTh_FER)),'descend') ; 
   
   for kk = 1 : size(eps_FixedTh_UER,1)
     eps_UER(kk) = mean(interp1(thetaHat_main_vec,eps_FixedTh_UER(kk,:),thetaDiff,'linear','extrap')) ;
     eps_FER(kk) = mean(interp1(thetaHat_main_vec,eps_FixedTh_FER(kk,:),thetaDiff,'linear','extrap')) ;
   end   
%% Search for undetected error probability target.
   if (METHOD_SEL == 1)
      xAxis = delta_vec ; 
   else
      xAxis = T_vec ; 
   end

   if eps_UER(end) > 1e-5
      xAxis_sel = xAxis(end); 
   elseif eps_UER(1) < 1e-5
      xAxis_sel = xAxis(1); 
   else
      xAxis_sel = interp1(eps_UER,xAxis,UER_target,'linear') ; 
   end

   if(xAxis_sel <= xAxis(1))
      eps_FER_sel = eps_FER(1) ; 
   elseif(xAxis_sel >= xAxis(end))
      eps_FER_sel = eps_FER(end) ; 
   else
      eps_FER_sel = interp1(xAxis,eps_FER,xAxis_sel,'linear') ; 
   end
%% SNR adjustment
   if(eps_FER_sel < FER_target) % If SNR is too high decrease it, 
      snrdB_vec = [snrdB_vec(1),snrdB] ; 
   else % If SNR is too low increase it, 
      snrdB_vec = [snrdB,snrdB_vec(2)] ; 
   end
   if(abs(diff(snrdB_vec)) < 0.1)
      SNR_out = snrdB_vec(2) ; 
      toc
      break ; 
   end
   toc
end

%% Save the data
if(METHOD_SEL == 1)
   delta_SNR = SNR_out ; 
   save_script = ['SNRSearch_QPSK_RCU_DeltaBit_n',num2str(n),',_k',num2str(k_sel),'.mat'] ; 
else
   thr_RCU_SNR = SNR_out ; 
   save_script = ['SNRSearch_QPSK_RCU_Threshold_n',num2str(n),'_k',num2str(k_sel),'_s',num2str(s),'.mat'] ; 
end
save(save_script) ; 
