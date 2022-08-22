function [vhat,iter,idx_true]=ldpc_decode_BPSK_nointrlv(rx_waveform,noise_var,H_check,H_channel)
% Optimized the cycle  some sum terms can be calculated in advance
% The calculation of message Q has benn optimized  11/16

dim    = size(H_check);
rows   = dim(1);
cols   = dim(2);   % L
[M,N]  = size(H_channel);
[h1i,h1j]=find(H_check==1);

% add_num = 0;
% H_add = randn(M,add_num)+1i*randn(M,add_num);
% H_channel = [H_channel,H_add];
% H_channel = H_channel(:,randperm(N+add_num));
% N=N+add_num;

y  = rx_waveform;  % normalize 7/6
P  = 0.5*ones(cols,N,M);
% mu = zeros(cols,N,M);
% var= zeros(cols,N,M);
A  = zeros(cols,N,M);
R  = zeros(N,rows,cols);
Q  = zeros(N,rows,cols);
L  = zeros(N,cols);
zero = zeros(rows,N);
sum_R = zeros(N,cols);
noise_var = noise_var/p_signal;  % normalize  7/7

for iter=1:30
    iter;
    
    %% 计算R
    for i=1:rows
        colind=find(h1i==i);
        colnum=length(colind);
        for j=1:colnum
            multi = 1;
            for k=1:colnum
                if k~=j
                    %除观测h1j(colind(j))的其它观测传递给校验i的消息
                    multi = multi.*tanh(Q(:,i,h1j(colind(k)))/2);   %size N 
                end
            end
            
            % 约束1：避免因atanh()函数造成R成为±inf，inf+(-inf)=NaN
%             multi_strict = (1-10^(-16)).*(multi>1-10^(-16))+multi.*(multi<=1-10^(-16)...
%                             & multi>=10^(-16)-1)+(10^(-16)-1).*(multi<10^(-16)-1);
                
            % 校验i传递给观测h1j(colind(j))的消息
            R(:,i,h1j(colind(j))) = 2*atanh(multi);   %size N
            R(R==inf)             = 19.07;
            R(R==-inf)            = -19.07;     % MATLAB Belief Propagation Decoding help   7/6
        end
    end
    
    %% 计算Q
    
     sum_A = sum(A,3); % size [cols,N] 7/6
    for j=1:cols
        rowind=find(h1j==j);
        rownum=length(rowind);
        sum_R(:,j) = sum(R(:,h1i(rowind),j),2);  % 11/16 
        for i=1:rownum
            Q_tmp = sum_R(:,j) - R(:,h1i(rowind(i)),j);   % 除校验h1i(rowind(i))的其它校验传递给变量j的消息  11/16

            % 变量j传递给校验h1i(rowind(i))的消息
            Q(:,h1i(rowind(i)),j) = Q_tmp + sum_A(j,:)';   %size N  7/6
            L(:,j) = Q(:,h1i(rowind(i)),j)+ R(:,h1i(rowind(i)),j);  % 变量j的本地LLR  size N
        end
    end

    
    %% 判决
    vhat = L>0;
    vhat = vhat*1;  %size  [N,cols]
    if (mod(H_check*vhat',2)==zero)
        if sum(sum(vhat))~=0
            break;
        end
    end 
    
    
    %% 计算A
    for i=1:M
        mu_tot  =  H_channel(i,:)*(2*P(:,:,i)-1).';  % size L 注意维度匹配
        var_tot = abs(H_channel(i,:)).^2 * 4 * ((1-P(:,:,i)).*P(:,:,i)).' + noise_var;  % size L 注意维度匹配
        for j=1:N
            %除变量j的其他变量传递给观测i的均值和方差
            mu_tmp = mu_tot -  (H_channel(i,j)*(2*P(:,j,i)-1)).';  % 后者 size L 注意维度匹配
            var_tmp = var_tot - (abs(H_channel(i,j))^2 * 4 * (1-P(:,j,i)).*P(:,j,i)).';  % 后者 size L 注意维度匹配
%             mu(:,j,i) = mu_tmp;
%             var(:,j,i) = var_tmp;   %debug阶段查看  引入不必要内存开销
            %观测i传递给变量j的消息
            A(:,j,i) = 4./var_tmp.'.*real(H_channel(i,j)'*(y(i,:).'-mu_tmp.'));  % 注意维度匹配  size L
        end
    end

    
    %% 计算P
    A_tot = sum(A,3);  % size [cols,N]  7/6
    for i=1:N
        for j=1:M
            sum_A = A_tot(:,i)-A(:,i,j);  % 7/6 
            tmp   = sum_A + sum_R(i,:).';  % 注意维度匹配
            
            % 约束2：避免P出现NaN  tmp过大将导致 exp(tmp)/(1+exp(tmp)) = inf/(inf) = NaN
            tmp_strict = 400*(tmp>400) + tmp.*(tmp<=400 & tmp>=-400) - 400*(tmp<-400);
            
            P(:,i,j) = exp(tmp_strict)./(1+exp(tmp_strict));   %size L
        end
    end
    
end

tmp      = mod(H_check*vhat',2); 
idx_true = find(sum(tmp==0,1)==size(tmp,1));  % return messgaes' index passed by LDPC check 7/6

