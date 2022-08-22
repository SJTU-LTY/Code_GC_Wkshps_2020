function [NMSE,x,gamma,p_false,p_mis,in_detect] = AMP_decoder(y,A,H_origin,idx_active)
% Y = AH + Z

p_total     = 1; %已归一化处理
iter        = 20;
[~,K]       = size(A);
[L,M_r]     = size (y);
[K_a,~]     = size(H_origin);
x           = zeros(K,M_r);
r           = y;
tao_2       = sum( abs(r).^2 , 1 ) / L / p_total ;    % MMV  对每根天线处分别计算，后续求和  8/14
sigma_h2    = ones(K,1);
lambda      = K_a/K;

for  j = 1 : iter
    
    % MMSE denoiser
    
    x_hat   =  A' * r / sqrt(p_total) + x ;
    phi     = 1 ./ (  1 + (1-lambda)/lambda *  exp( sum( -bsxfun(@minus, 1./tao_2 ,1./bsxfun(@plus,sigma_h2,tao_2) ).* x_hat.* conj(x_hat),2 )+ ... 
     sum(log( 1+bsxfun(@rdivide,sigma_h2,tao_2)),2))); 
    phi     = real(phi) ;
    thres   = bsxfun(@rdivide, real(phi) .* sigma_h2 , bsxfun(@plus,sigma_h2 ,tao_2));  
    
    % denoiser
    x       = thres .* x_hat ;
                           
   % linear operation
    cc  = -(1-lambda)/lambda * bsxfun(@times,phi.^2 .* exp( sum(  -bsxfun(@minus, 1./tao_2 ,1./bsxfun(@plus,sigma_h2,tao_2)).* x_hat.* conj(x_hat),2 )+ ... 
                sum(log( 1+bsxfun(@rdivide,sigma_h2,tao_2) ),2)  ) , conj(x_hat)) .* -bsxfun(@minus, 1./tao_2 ,1./bsxfun(@plus,sigma_h2,tao_2))...
                .* bsxfun(@rdivide,sigma_h2 , bsxfun(@plus,sigma_h2 , tao_2)) .* x_hat ;      % 最先出现 NaN   7/7
                
    cc(isnan(cc)) = 0+0*1i; 
    r         =  y - sqrt(p_total)* A * x + 1/L * bsxfun(@times,sum( thres + abs(cc) ,1 ), r) ;
    
    tao_2_old   = tao_2 ;
    tao_2       = sum( abs(r).^2 , 1 ) / L / p_total ;
    
    % termination criterion NMSE
    if abs( norm(tao_2_old,'fro')^2 - norm(tao_2,'fro')^2 )/ norm(tao_2,'fro')^2 < 1e-3
        break;
    end
           
    if norm(tao_2,'fro') >norm(tao_2_old,'fro')
            break;
    end
end

x_power         = sort(sum(abs(x),2),'descend');
[gamma,~]       = mapminmax(x_power',0,1);     
gamma           = gamma';

soft_thres = log( prod(1+bsxfun(@rdivide,sigma_h2,tao_2),2) );
user_det   = diag(bsxfun(@minus, 1./tao_2 , 1./bsxfun(@plus,sigma_h2 , tao_2) ).* x_hat * x_hat');  

if ~sum(soft_thres==inf)
    in_detect = find (user_det>soft_thres);   % size [K,1]
else
    in_detect = find(gamma>0.5);
end
in_true     = idx_active;  % size: [Ka,1]

% false alarm
false1      = ismember(in_detect , in_true);
if isempty(in_detect) 
    p_false = sum(false1==0) / length(in_true) ;
else
    p_false = sum(false1==0) / length(in_detect) ;
end

% miss detection
mis1 = ismember(in_true, in_detect);
p_mis = sum(mis1==0) / length(in_true) ; 

% NMSE
NMSE = norm(x(idx_active,:)-H_origin,'fro')^2/norm(H_origin,'fro')^2;
end