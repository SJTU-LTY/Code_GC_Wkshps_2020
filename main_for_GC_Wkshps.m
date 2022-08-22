%% Prepare the workspace && Set parameters 
% Simulation for GC-Wkshps 2020, DOI: 10.1109/GCWkshps50303.2020.9367450
% Part of the simulation for IEEE JSAC, DOI: 10.1109/JSAC.2022.3145914
% Contact: tianya@sjtu.edu.cn
% Last change: 22/08/2022

clc; clear all; close all
path(path, '.\LDPC')

%% Basic setting
Ka_list     = 60:2:70;     % number of active UDs
M           = 16;          % number of BS antennas
L           = 200;         % Length of pilots
SNR         = 6; % Es/N0
num_iter    = 50;

modu        = 2; % QPSK 
CS_length   = 12;
code_length = 1e3;
tot_length  = code_length-L;
K_tot       = 2^CS_length; % number of potential UDs in the cell

data_length   = 100;

for idx_Ka = 1:length(Ka_list)
    
    p_fa        = 0;
    p_md        = 0;
    NMSE1       = 0;
    pe_ad       = 0;
    
    K_a         = Ka_list(idx_Ka);
    cols        = 200;   %  No. cloumns in the LDPC check matrix
    rows        = cols - data_length;   %  No. rows in the LDPC check matrix
    P_dB        = SNR + 10*log10(code_length/(L+ cols/modu));
        
    for count = 1:num_iter

%% Data generation
        B_1 = randi([0 1],K_a,CS_length); 
        B_2 = randi([0 1],K_a,data_length);
        B   = [B_1,B_2];  % size [K_a,tot_length]
        idx_active  = randperm(K_tot,K_a)';  % collision is not considered
     
%% LDPC encoding
        H_check = genH(rows,cols);
        [u,P,rearranged_cols] = ldpc_encode(B_2,H_check);
%% Mudulation \ Zero padding \ Interleaver \ AMP estimation
        H = (randn(K_a,M)+1i*randn(K_a,M))/sqrt(2); 
        A = (1 + 1i)/sqrt(2)*exp(2*pi*1i*rand(L,K_tot))/sqrt(L);  %　Gaussion pilots
        X_modu = 2 * u - 1; % BPSK  size [K_a,cols]
        if modu==2           
            X_modu = 1/sqrt(2)*(X_modu(:,1:2:cols-1) + sqrt(-1)*X_modu(:,2:2:cols)); % QPSK  size [K_a,cols/2]
        end           
        
        X_modu(:,end+1:tot_length) = 0; %　zero padding 
        for idx=1:size(X_modu,1)
            X_modu(idx,:) = randintrlv(X_modu(idx,:),idx_active(idx));  % rand interleave
        end        

%% Channel Estimation
        Y_1 = awgn(A(:,idx_active)*H, P_dB,'measured');
        [NMSE,x,~,~,~,idx_active_hat] = AMP_decoder(Y_1,A,H,idx_active);
        H_est = x(idx_active_hat,:).'; % size [M,K]
        Pe    = length(setdiff(idx_active,idx_active_hat)) / K_a;
              
%% LDPC decoding
        Y_2 = awgn(H.'*X_modu, P_dB, 'measured'); 
        P_n = 10.^(-P_dB/10); 

        idx_next = 1:length(idx_active_hat);
        [vhat,vhat_no_sic,idx_recv,~,~] = LDPC_SIC (Y_2,P_n,H_check,H_est,idx_next,idx_active_hat,modu);
  
%% error count     
        if isempty(vhat)
            fa_LDPC      = 0;
            md_LDPC      = 1;
        else         
            uhat         = extract_mesg(vhat,rearranged_cols); %　size [K,cols-rows]
            [md_LDPC,fa_LDPC] = error_count(B_2,uhat);
        end

        p_fa         = p_fa + fa_LDPC/num_iter;
        p_md         = p_md + md_LDPC/num_iter;
        NMSE1        = NMSE1 + NMSE/num_iter;
        pe_ad        = pe_ad + Pe/num_iter;
        
    fprintf ('%d %.4f %.4f %.4f %.4f\n',count,md_LDPC,fa_LDPC,NMSE,Pe); 
    end
    fprintf ('%d %.4f %.4f %.4f %.4f\n',Ka_list(idx_Ka),p_md,p_fa,NMSE1,pe_ad);   
end
  

 




