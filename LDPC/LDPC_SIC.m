function [Vhat,vhat_no_sic,idx_recv,count,idx_next] = LDPC_SIC (Y_2,noise_var,H_check,H_channel,idx_next,idx_active_hat,modu)
% SIC部分始终存在问题，SIC与原版本相比不成功  5/19

[~,tot_length] = size (Y_2);
[~,cols] = size (H_check);

if modu==1
    [vhat,iter,idx_true] = ...
                ldpc_decode_BPSK(Y_2,noise_var,H_check,H_channel(:,idx_next),idx_active_hat);   % single decoding
            % 由于数据填0并按填0前的长度归一化，实际的等效噪声能量与原噪声存在尺度变换 2021/5/18
elseif modu==2
    [vhat,iter,idx_true] = ...
                ldpc_decode_QPSK(Y_2,noise_var,H_check,H_channel(:,idx_next),idx_active_hat);   % single decoding
elseif modu==4
    [vhat,iter,idx_true] = ...
                ldpc_decode_16QAM(Y_2,noise_var,H_check,H_channel(:,idx_next),idx_active_hat);  % single decoding

end
    vhat_no_sic = vhat(idx_true,:);
    tmp_idx     = idx_active_hat(idx_next);
    tmp_H       = H_channel(:,idx_next);
    idx_recv    = [];
    idx_copy    = idx_active_hat;
    Vhat        = [];
    count       = 1;
    while (1)
        if isempty(idx_true)
            break;
        end
        pos = [];
        for idx=1:length(idx_true)
            [~,pos_tmp] = ismember(tmp_idx(idx_true(idx)),idx_copy);
            idx_copy(pos_tmp) = 0;
            pos = [pos;pos_tmp];
        end
        idx_next      = setdiff(idx_next,pos);   % original position 
        idx_recv      = [idx_recv,pos'];   % original position in idx_active_hat
        Vhat          = [Vhat;vhat(idx_true,:)];
        if  isempty(idx_next)
            break;
        end
         
        if modu == 1
            tmp = 2 * vhat(idx_true,:)-1;
            tmp(:,end+1:tot_length) = 0;
            for idx=1:size(tmp,1)
                tmp(idx,:) = randintrlv(tmp(idx,:),tmp_idx(idx_true(idx)));
            end
            Y_2          = Y_2 - tmp_H(:,idx_true) * tmp;
%             Y_2          = Y_2 - tmp_H(:,idx_true) * tmp;
            tmp_H        = H_channel(:,idx_next);
            tmp_idx      = idx_active_hat(idx_next);   
            [vhat,iter,idx_true] =  ldpc_decode_BPSK(Y_2,noise_var,H_check,H_channel(:,idx_next),tmp_idx);

        elseif modu == 2
            tmp =  1/sqrt(2)*(2*vhat(idx_true,1:2:cols-1)-1 + sqrt(-1)*(2*vhat(idx_true,2:2:cols)-1));
            tmp(:,end+1:tot_length) = 0;
            for idx=1:size(tmp,1)
                tmp(idx,:) = randintrlv(tmp(idx,:),tmp_idx(idx_true(idx)));
            end     
            Y_2          = Y_2 - tmp_H(:,idx_true) * tmp;
            tmp_H        = H_channel(:,idx_next);
            tmp_idx      = idx_active_hat(idx_next);   
            [vhat,iter,idx_true] =  ldpc_decode_QPSK(Y_2,noise_var,H_check,H_channel(:,idx_next),tmp_idx);    
        elseif modu == 4
            tmp = QAM16_modu (vhat(idx_true,:));
            tmp(:,end+1:tot_length) = 0;
            for idx=1:size(tmp,1)
                tmp(idx,:) = randintrlv(tmp(idx,:),tmp_idx(idx_true(idx)));
            end 
            Y_2          = Y_2 - tmp_H(:,idx_true) * tmp;
            tmp_H        = H_channel(:,idx_next);
            tmp_idx      = idx_active_hat(idx_next);  
        end    
      count   = count+1;
    end          
end
