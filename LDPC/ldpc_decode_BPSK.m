function [vhat,iter,idx_true]=ldpc_decode_BPSK(rx_waveform,noise_var,H_check,H_channel,idx_active_hat)
% Optimized the cycle  some sum terms can be calculated in advance
% The calculation of message Q has benn optimized  11/16
% 更新为填0后交织版本 5/18/2021

dim    = size(H_check);
rows   = dim(1);
cols   = dim(2);   % L
[M,N]  = size(H_channel);
[h1i,h1j]=find(H_check==1);

%　creat alphabet
tabel = repmat(1:cols,N,1);
tabel(:,end+1:size(rx_waveform,2)) = 0;
for idx=1:N
    tabel(idx,:) = randintrlv(tabel(idx,:),idx_active_hat(idx));
end
[tli,tlj]=find(tabel~=0);  % 坐标为交织后位置，坐标内容为交织前位置

y  = rx_waveform;
P  = 0.5*ones(cols,N,M);
% mu = zeros(cols,N,M);
% var= zeros(cols,N,M);
A  = zeros(cols,N,M);
R  = zeros(N,rows,cols);
Q  = zeros(N,rows,cols);
L  = zeros(N,cols);
zero = zeros(rows,N);
sum_R = zeros(N,cols);

for iter=1:30
    
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
        
    %% 计算A  A和P均为解交织后的结果，即代入计算时为交织前的数组
    for i=1:M
        for j=1:size(y,2)  % tot_length
           rowind = find(tlj==j);  % positions of users
           rownum = length(rowind);
           mu_tot  = H_channel(i,tli(rowind))*diag(2*P(tabel(tli(rowind),j),tli(rowind),i)-1);  % scalar
           %　diag 可以将提取矩阵行和列交集转化为提取矩阵对应位置的元素
           % tli(rowind) 非0数据的用户
           % tabel(tli(rowind),j) 用户交织前的数据位置
           var_tot = abs(H_channel(i,tli(rowind))).^2 * ...
                    4 * (diag(1-P(tabel(tli(rowind),j),tli(rowind),i)).*diag(P(tabel(tli(rowind),j),tli(rowind),i))) + noise_var;  % scalar
            for k=1:rownum
                mu_tmp = mu_tot - H_channel(i,tli(rowind(k)))*...
                        (2*P(tabel(tli(rowind(k)),j),tli(rowind(k)),i)-1);
                var_tmp = var_tot - abs(H_channel(i,tli(rowind(k))))^2 * 4 * ...
                            (1-P(tabel(tli(rowind(k)),j),tli(rowind(k)),i)) * P(tabel(tli(rowind(k)),j),tli(rowind(k)),i);
                A(tabel(tli(rowind(k)),j),tli(rowind(k)),i) = ...
                        4./var_tmp.'.*real(H_channel(i,tli(rowind(k)))'*(y(i,j).'-mu_tmp.'));
            end
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
