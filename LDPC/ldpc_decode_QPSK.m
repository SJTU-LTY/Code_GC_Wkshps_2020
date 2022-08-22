function [vhat,iter,idx_true] = ...
    ldpc_decode_QPSK (Y_2,noise_var,H_check,H_channel,idx_active_hat)
% 更新为填0后交织版本 6/27/2021

[rows,cols] = size(H_check);    %　cols = L
[M,N]       = size(H_channel);  % size [M,K]
[h1i,h1j]   = find(H_check==1);
modu        = 2;

%　creat alphabet
tabel = repmat(1:cols/modu,N,1);
tabel(:,end+1:size(Y_2,2)) = 0;
for idx=1:N
    tabel(idx,:) = randintrlv(tabel(idx,:),idx_active_hat(idx));
end
[tli,tlj]=find(tabel~=0);  % 坐标为交织后位置，坐标内容为交织前位置

P_r         = 0.5*ones(cols/2,N,M);  %　QPSK modulation, at symbol
P_i         = 0.5*ones(cols/2,N,M);  %　QPSK modulation, at symbol
A           = zeros(cols/2,N,M);     %　QPSK modulation, at symbol  real and imag parts
R           = zeros(N,rows,cols); 
Q           = zeros(N,rows,cols);
L           = zeros(N,cols);
zero        = zeros(rows,N);
sum_R       = zeros(N,cols);
norm        = 1/sqrt(2); 

for iter = 1:30
    
    % 计算R
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
    
    % 计算Q
    sum_A = sum(A,3); % size [cols,N] 7/6
    for j=1:cols
        rowind=find(h1j==j);
        rownum=length(rowind);
        for i=1:rownum
            Q_tmp = 0;
            for k=1:rownum
                if k~=i
                    % 除校验h1i(rowind(i))的其它校验传递给变量j的消息
                    Q_tmp = Q_tmp + R(:,h1i(rowind(k)),j);   %size N;
                end
            end
            % 变量j传递给校验h1i(rowind(i))的消息
            if mod (j,2)~=0
                Q(:,h1i(rowind(i)),j) = Q_tmp + real(sum_A((j+1)/2,:))';   % 奇  size N  7/6  QPSK  \gamma, P 分实虚  Q, R 分奇偶。 
            else
                Q(:,h1i(rowind(i)),j) = Q_tmp + imag(sum_A(j/2,:))';   % 偶
            end 
            sum_R(:,j) = Q_tmp + R(:,h1i(rowind(i)),j);
            L(:,j) = Q(:,h1i(rowind(i)),j)+ R(:,h1i(rowind(i)),j);  % 变量j的本地LLR  size N
        end
    end
    
    % 判决
    vhat = L>0;  % 调制前
    vhat = vhat*1;  %size  [N,cols] 
    if (mod(H_check*vhat',2)==zero)
        if sum(sum(vhat))~=0
            break;
        end
    end
    
    % 计算A  仅此处用到了 H_channel, Y_2
    for i=1:M
        for j=1:size(Y_2,2)
           rowind = find(tlj==j);  % positions of users
           rownum = length(rowind);
           mu_tot = norm * H_channel(i,tli(rowind)) * ...
                diag(2*P_r(tabel(tli(rowind),j),tli(rowind),i)-1 + sqrt(-1)*(2*P_i(tabel(tli(rowind),j),tli(rowind),i)-1)); % scalar
           sigma_tot = 2 * diag(P_r(tabel(tli(rowind),j),tli(rowind),i)-P_r(tabel(tli(rowind),j),tli(rowind),i).^2 + ...
                P_i(tabel(tli(rowind),j),tli(rowind),i)-P_i(tabel(tli(rowind),j),tli(rowind),i).^2);  % size len(rowind)
           var_tot   = abs(H_channel(i,tli(rowind))).^2 * sigma_tot + noise_var;  % scalar
           
           for k=1:rownum
               mu_tmp  = mu_tot -  norm * H_channel(i,tli(rowind(k)))* ...
                   (2*P_r(tabel(tli(rowind(k)),j),tli(rowind(k)),i)-1 + ...
                   sqrt(-1)*(2*P_i(tabel(tli(rowind(k)),j),tli(rowind(k)),i)-1)); 
               var_tmp = var_tot - abs(H_channel(i,tli(rowind(k))))^2 * sigma_tot(k);  % 11/19/2021 sigma(rowind(k)) to sigma(k)
               A(tabel(tli(rowind(k)),j),tli(rowind(k)),i) = ...
                    2 * sqrt(2)./var_tmp.* (H_channel(i,tli(rowind(k)))'*(Y_2(i,j)-mu_tmp)); 
           end
        end
    end
    
    % 计算P
    A_tot = sum(A,3);  % size [cols,N]  7/6
    for i=1:N
        for j=1:M
            sum_A = A_tot(:,i)-A(:,i,j);  % 7/6 
            tmp_r = real(sum_A) + sum_R(i,1:2:cols-1).';  % 注意维度匹配 
            tmp_i = imag(sum_A) + sum_R(i,2:2:cols).';  % 注意维度匹配 
            
            % 约束2：避免P出现NaN  tmp过大将导致 exp(tmp)/(1+exp(tmp)) = inf/(inf) = NaN
            tmp_r = 100*(tmp_r>100) + tmp_r.*(tmp_r<=100 & tmp_r>=-100) - 100*(tmp_r<-100);
            tmp_i = 100*(tmp_i>100) + tmp_i.*(tmp_i<=100 & tmp_i>=-100) - 100*(tmp_i<-100);
            
            P_r(:,i,j) = exp(tmp_r)./(1+exp(tmp_r));   %size L
            P_i(:,i,j) = exp(tmp_i)./(1+exp(tmp_i));   %size L
        end
    end   
end

tmp      = mod(H_check*vhat',2); 
idx_true = find(sum(tmp==0,1)==size(tmp,1));  % return messgaes' index passed by LDPC check 7/6
% uhat     = extract_mesg(vhat,rearranged_cols);  %注意维度匹配   size  [N,cols-rows] 调制前 编码前

% L_scalar   = 100*(L>100) + L.*(L<=100 & L>=-100) - 100*(L<-100);
% P_final    = exp(L_scalar)./(1+exp(L_scalar));  % size [N,cols]
% mu_tot     = norm*(2*P_final(:,1:2:cols-1)-1 + sqrt(-1)*(2*P_final(:,2:2:cols)-1)); % size [N,cols/2]
% sigma_tot  = 2*(P_final(:,1:2:cols-1)-P_final(:,1:2:cols-1).^2 + P_final(:,2:2:cols)-P_final(:,2:2:cols).^2); % size [N,cols/2]


end