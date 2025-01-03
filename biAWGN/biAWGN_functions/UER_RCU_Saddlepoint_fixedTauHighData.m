function [eps_und, eps_tot]  = UER_RCU_Saddlepoint_fixedTauHighData(snr,n,rate,N,T)

% INPUTS:
% snr   = Value of the snr in linear scale (no dB!)
% rate = Values of Rate at which the bound will be computed
% n     = Blocklength
% T     = Parameter to trade-off between UER and FER (\lambda in paper).  T = -inf -> UER = FER
%
% OUTPUT:
% eps_und     = Undetected error probability 
% eps_tot     = Total error probability 

s = 2.5;  
debug =0 ; 
if debug == 1
   SNR_db = 3; 
   snr = DB_to_lin(SNR_db); 
   n = 64 ; 
   rate = 0.5 ; 
   T = 0.5 ; 
   N = 1e4 ; 
end
LH_func = @(x,y) 1/sqrt(2.*pi) .* exp(-(y-x).^2/2) ; 
LH_f1 = @(y) 1/sqrt(2.*pi) .* exp(-(y-sqrt(snr)).^2/2) ; 
LH_f2 = @(y) 1/sqrt(2.*pi) .* exp(-(y+sqrt(snr)).^2/2) ; 
LLR = @(x,y) sum(log(LH_func(x,y)),2) ; 
qMeanLog_xbar = @(Y,s) sum(log(0.5.* 1/sqrt(2.*pi).*(exp(-(s/2)*(Y+sqrt(snr)).^2) + exp(-(s/2)*(Y-sqrt(snr)).^2) ) ),2) ; 
% omega_CGF = @(yi,tau) sum(log(0.5.*(LH_func(sqrt(snr),yi).^(tau) + LH_func(-sqrt(snr),yi).^(tau) )),2)  ;
omega_Deriv1 = @(y,tau) sum( (LH_f1(y).^tau .* log(LH_f1(y)) + LH_f2(y).^tau .* log(LH_f2(y)) ) ./ (LH_f1(y).^tau + LH_f2(y).^tau) , 2 ) ;
omega_Deriv2 = @(y,tau) sum( LH_f1(y).^tau .* LH_f2(y).^tau .* (log(LH_f2(y)) - log(LH_f1(y)) ).^2 ./ (LH_f2(y).^tau + LH_f1(y).^tau).^2  ,2) ; 
InfDens = @(X,Y,s) sum(log(2) - (s/2)*(Y-X).^2 - log(exp(-(s/2)*(Y+sqrt(snr)).^2) + exp(-(s/2)*(Y-sqrt(snr)).^2)),2); %compute information density samples
InfDensTerm = @(X,Y,s) log(2) - (s/2)*(Y-X).^2 - log(exp(-(s/2)*(Y+sqrt(snr)).^2) + exp(-(s/2)*(Y-sqrt(snr)).^2));


% xiBar = (randi([0,1],N_CGF,n).*2 - 1).*sqrt(snr) ;
% LF_xy = LH_func(xiBar,yi) ; 

Zmax = 10 ; 
rho_current = NaN ;
rho_vec_in1 = -0.9 ; 
rho_vec_in2 = 5 ; 
while isnan(rho_current)
   rho_vec = linspace(rho_vec_in1,rho_vec_in2,200);
   s_vec=1./(1+rho_vec);
   E0prime_vec = nan(size(rho_vec));
   for ii=1:length(rho_vec)
       E0prime_vec(ii) = E0prime(snr,rho_vec(ii),s_vec(ii),Zmax);
   end
   rem_idx = find(isnan(E0prime_vec));
   E0prime_vec(rem_idx) = [];
   rho_vec(rem_idx) = [];
   E0prime_vec = E0prime_vec + sort(1e-10.*randn(size(E0prime_vec)),'descend') ;
   rate_nats = rate.*log(2);
   rho_current = interp1(E0prime_vec,rho_vec,rate_nats);
   tauSel = 1./(1+rho_current) ; 
   rho_vec_in1 = rho_vec_in1 + 5.9  ; 
   rho_vec_in2 = rho_vec_in2 + 5.9 ; 
end

batchNumLim = ceil(N/1e5) ; 
N = round(N/1e5) ; 
for batchNum = 1 : batchNumLim 
   xi = (randi([0,1],N,n).*2 - 1).*sqrt(snr) ;
   yi = xi + randn(N,n) ; 
   W_xi = LLR(xi,yi) ; 
   id_tau = InfDens(xi,yi,tauSel);

   omega_tau_d1 = omega_Deriv1(yi,tauSel) ; 
   omega_tau_d2 = omega_Deriv2(yi,tauSel) ;  
   omega_tau_d2(omega_tau_d2<0) = 0 ; 
   preExp_term1 = exp(tauSel.*(W_xi - omega_tau_d1 ) +0.5.*tauSel.^2.*omega_tau_d2 ) ;
   preExp_term2 = 0.5.*erfc((W_xi-omega_tau_d1 + tauSel.*omega_tau_d2)./(sqrt(2.*omega_tau_d2))) ; 
   pep = exp(-id_tau).*preExp_term1.*preExp_term2 ; 
  
   preExp2_term1 = exp(tauSel.*(max(W_xi,(n.*T + qMeanLog_xbar(yi,s))/s) - omega_tau_d1 ) +0.5.*tauSel.^2.*omega_tau_d2 ) ;
   preExp2_term2 = 0.5.*erfc((max(W_xi,(n.*T + qMeanLog_xbar(yi,s))/s)-omega_tau_d1 + tauSel.*omega_tau_d2)./(sqrt(2.*omega_tau_d2))) ; 
   pep2 = exp(-id_tau).*preExp2_term1.*preExp2_term2 ; 
   id =  InfDens(xi,yi,s) ; 
   R = rate ; 
   M = 2.^(n.*R) ; 
   % [id]= info_dens_biawgn(snr,n,1,N) ;
   eps_tot1(batchNum,:) =  mean( exp( - max(0, -log(pep) - log(M-1) - log(double(id>n*T))  ) ),1);   
   eps_tot2(batchNum,:) = mean(id <= n*T,1) ; 
   eps_und(batchNum,:) = mean( exp( - max(0, -log(pep2) - log(M-1)) ),1);  
end
eps_tot = mean(eps_tot1,1) + mean(eps_tot2,1) ; 
eps_und = mean(eps_und,1) ; 
eps_und = eps_und + sort(1e-10.*rand(size(eps_und)),'descend') ; 
  
end

% eps_erasure = mean(idMaxVec <= n.*T,1) ; 
% eps_und = eps_tot - eps_erasure ; 
function res=E0prime(snr,rho,s,Zmax) 
  
 f0 = @(z) exp(-z.^2/2)/sqrt(2*pi) .* exp(-rho*(log(2) - log(1+exp(-2*s*(z*sqrt(snr)+snr)))));
 
 f1= @(z) exp(-z.^2/2)/sqrt(2*pi) .*(log(2) - log(1+exp(-2*s*(z*sqrt(snr)+snr)))).* exp(-rho*(log(2) - log(1+exp(-2*s*(z*sqrt(snr)+snr)))));
 
 res= integral(f1, -Zmax, Zmax)/integral(f0,-Zmax,Zmax);
 
end
