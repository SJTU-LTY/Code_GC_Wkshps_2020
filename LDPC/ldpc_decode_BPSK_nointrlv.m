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
    
    %% ����R
    for i=1:rows
        colind=find(h1i==i);
        colnum=length(colind);
        for j=1:colnum
            multi = 1;
            for k=1:colnum
                if k~=j
                    %���۲�h1j(colind(j))�������۲⴫�ݸ�У��i����Ϣ
                    multi = multi.*tanh(Q(:,i,h1j(colind(k)))/2);   %size N 
                end
            end
            
            % Լ��1��������atanh()�������R��Ϊ��inf��inf+(-inf)=NaN
%             multi_strict = (1-10^(-16)).*(multi>1-10^(-16))+multi.*(multi<=1-10^(-16)...
%                             & multi>=10^(-16)-1)+(10^(-16)-1).*(multi<10^(-16)-1);
                
            % У��i���ݸ��۲�h1j(colind(j))����Ϣ
            R(:,i,h1j(colind(j))) = 2*atanh(multi);   %size N
            R(R==inf)             = 19.07;
            R(R==-inf)            = -19.07;     % MATLAB Belief Propagation Decoding help   7/6
        end
    end
    
    %% ����Q
    
     sum_A = sum(A,3); % size [cols,N] 7/6
    for j=1:cols
        rowind=find(h1j==j);
        rownum=length(rowind);
        sum_R(:,j) = sum(R(:,h1i(rowind),j),2);  % 11/16 
        for i=1:rownum
            Q_tmp = sum_R(:,j) - R(:,h1i(rowind(i)),j);   % ��У��h1i(rowind(i))������У�鴫�ݸ�����j����Ϣ  11/16

            % ����j���ݸ�У��h1i(rowind(i))����Ϣ
            Q(:,h1i(rowind(i)),j) = Q_tmp + sum_A(j,:)';   %size N  7/6
            L(:,j) = Q(:,h1i(rowind(i)),j)+ R(:,h1i(rowind(i)),j);  % ����j�ı���LLR  size N
        end
    end

    
    %% �о�
    vhat = L>0;
    vhat = vhat*1;  %size  [N,cols]
    if (mod(H_check*vhat',2)==zero)
        if sum(sum(vhat))~=0
            break;
        end
    end 
    
    
    %% ����A
    for i=1:M
        mu_tot  =  H_channel(i,:)*(2*P(:,:,i)-1).';  % size L ע��ά��ƥ��
        var_tot = abs(H_channel(i,:)).^2 * 4 * ((1-P(:,:,i)).*P(:,:,i)).' + noise_var;  % size L ע��ά��ƥ��
        for j=1:N
            %������j�������������ݸ��۲�i�ľ�ֵ�ͷ���
            mu_tmp = mu_tot -  (H_channel(i,j)*(2*P(:,j,i)-1)).';  % ���� size L ע��ά��ƥ��
            var_tmp = var_tot - (abs(H_channel(i,j))^2 * 4 * (1-P(:,j,i)).*P(:,j,i)).';  % ���� size L ע��ά��ƥ��
%             mu(:,j,i) = mu_tmp;
%             var(:,j,i) = var_tmp;   %debug�׶β鿴  ���벻��Ҫ�ڴ濪��
            %�۲�i���ݸ�����j����Ϣ
            A(:,j,i) = 4./var_tmp.'.*real(H_channel(i,j)'*(y(i,:).'-mu_tmp.'));  % ע��ά��ƥ��  size L
        end
    end

    
    %% ����P
    A_tot = sum(A,3);  % size [cols,N]  7/6
    for i=1:N
        for j=1:M
            sum_A = A_tot(:,i)-A(:,i,j);  % 7/6 
            tmp   = sum_A + sum_R(i,:).';  % ע��ά��ƥ��
            
            % Լ��2������P����NaN  tmp���󽫵��� exp(tmp)/(1+exp(tmp)) = inf/(inf) = NaN
            tmp_strict = 400*(tmp>400) + tmp.*(tmp<=400 & tmp>=-400) - 400*(tmp<-400);
            
            P(:,i,j) = exp(tmp_strict)./(1+exp(tmp_strict));   %size L
        end
    end
    
end

tmp      = mod(H_check*vhat',2); 
idx_true = find(sum(tmp==0,1)==size(tmp,1));  % return messgaes' index passed by LDPC check 7/6

