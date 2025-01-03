function [epsFER_Forney,epsUER_Forney] = ForneyAchievError_biAWGN(n,rate,SNR_list,T)

% n is blocklength
% 'SNR_list' is the list of SNR and should be in dB
% 'rate' is rate in bits per channel use
% 'T' is the parameter to trade-off between UER and FER. For T=0 FER=UER.  
 

quant_bitNum = 12;

quan_upLim = 11 ; 
quan_downLim = -11 ; 
quan_lims = fliplr(linspace(quan_upLim,quan_downLim,2.^quant_bitNum-1)) ; 
bit_seq = 1:2.^quant_bitNum; 

data_num =1e6; 
q_k = 0.5; 



for LOOPER = 1:length(SNR_list)
   SNR = 10.^(SNR_list(LOOPER)./10) ; 
   for kk = 1 : length(quan_lims)+1
      if(kk == 1)
         transitionProb(kk) = 1-qfunc(quan_lims(kk) - sqrt(SNR)) ; 
      elseif(kk == length(quan_lims)+1)
         transitionProb(kk) = qfunc(quan_lims(kk-1) - sqrt(SNR)) ; 
      else
         transitionProb(kk) = qfunc(quan_lims(kk-1) - sqrt(SNR)) - qfunc(quan_lims(kk) - sqrt(SNR))  ; 
      end
   end

   p_list = linspace(0.01,1,15) ; 
   E1_RT_max = -inf ; 
   for ii = 1 : length(p_list)
      p_sel = p_list(ii) ; 
      s_list = linspace(0,p_sel,15) ; 
      for jj = 1 : length(s_list)
         s_sel = s_list(jj) ; 
         E1_RT = Eval_E0(q_k,transitionProb,s_sel,p_sel) - p_sel.*rate - s_sel.*T ; 
         if(E1_RT_max <= E1_RT)
            E1_RT_max = E1_RT ; 
            s = s_sel ; 
            p = p_sel ; 
         end
      end
   end
   E2_RT = E1_RT_max + T ; 
   epsFER_Forney(LOOPER) = 2.^(-n.*E1_RT_max) ; % FER: Frame error rate
   epsUER_Forney(LOOPER) = 2.^(-n.*E2_RT) ; % UER: Undetected error probability
end

end

function E0_sp = Eval_E0(q_k,transitionProb,s,p)
   transitionProb_1 = transitionProb; 
   transitionProb_0 = fliplr(transitionProb) ; 
   for jj = 1 : length(transitionProb)
      E0_term1 = (transitionProb_1(jj).^(1-s).*q_k + transitionProb_0(jj).^(1-s).*q_k) ;
      E0_term2 = (transitionProb_1(jj).^(s/p).*q_k + transitionProb_0(jj).^(s/p).*q_k).^(p) ; %(q_k.*flipProb.^(s/p) + q_k.*(1-flipProb).^(s/p)).^(p) ; 
      E0_temp(jj) = (E0_term1 .* E0_term2) ; 
   end
   E0_sp = -log2(sum(E0_temp)) ; 
end


% [epsFER_Gallager,epsUER_Gallager] = GallagerAchievError(64,0.5,0.1,2)
% transitionProb= [flipProb.^2, flipProb.*(1-flipProb), flipProb.*(1-flipProb), (1-flipProb).^2 ]; 





