function [ delta_hat, C_hat, significant ] = func_effect_twogroup( y0, y1, alpha_level )
%% TODO
%Dave Zachariah

%% Intialize
%Global variables
global n0
global n1
global vec0
global vec1
global y0_sum
global y1_sum
global y0_sq
global y1_sq


%Fixed
n0     = length(y0);
n1     = length(y1);
vec0 = ones(n0,1);
vec1 = ones(n1,1);
y0_sum = sum(y0);
y1_sum = sum(y1);
y0_sq  = sum(abs(y0).^2);
y1_sq  = sum(abs(y1).^2);

%Initial
delta_hat = 0;
mu_hat    = mean(y0);
lambda    = 0;
v0        = var(y0,1);
v1        = var(y1,1);
C_hat     = zeros(1,2);

%% Compute estimates

%Search minimizing mu
[mu_hat,rho,alpha,~,gamma] = func_1dsearch(mu_hat,sqrt(v0));


%Theta
v0        = (y0_sq + (mu_hat^2*n0) - 2*mu_hat*y0_sum) / n0;
v1        = (alpha + rho * gamma)/(n1 + rho*n1^2);
lambda    = rho*v1;
%theta_hat = [mu_hat, lambda, v0, v1]';

%Delta
delta_hat = lambda*(y1_sum - mu_hat*n1) / (lambda*n1 + v1);

%CRB
crb_hat = (lambda*v1)/(lambda*n1 + v1) + ...
           (  (lambda*n1)/(lambda*n1 + v1)  )^2 ...
           /( n0/v0  +  n1/(lambda*n1 + v1) );

%Confidence interval
C_hat(1) = delta_hat_ebm(k,m) - sqrt(1/alpha_level)*sqrt( crb_hat );
C_hat(2) = delta_hat_ebm(k,m) + sqrt(1/alpha_level)*sqrt( crb_hat );

%Significant at alpha-level
significant = ((C_hat(2) > 0) & (C_hat(1) > 0)) | ((C_hat(2) < 0) & (C_hat(1) < 0));
       
end


%%%%%%%%%%%%%%%%%%%%%%%
% Search minimizing mu
%%%%%%%%%%%%%%%%%%%%%%%

function [mu,rho,alpha,beta,gamma] = func_1dsearch(mu, s0)

global n0
global n1
%global vec0
%global vec1
global y0_sum
global y1_sum
global y0_sq
global y1_sq

%Initialize
bnd  = 6*s0;
M    = 200;
lvl  = 20;

%Refined line search
f_set = zeros(M,1);
for l = 1:lvl
    
    %Create set
    m_set = linspace( mu-bnd/l, mu+bnd/l, M);
    
    %Evaluate
    k = 1;
    for m = m_set
        
        alpha0 = y0_sq    + (m^2*n0) - 2*m*y0_sum;    %||y0 - m * vec0||^2;
        alpha  = y1_sq    + (m^2*n1) - 2*m*y1_sum;    %norm(y1 - m * vec1 )^2;
        beta   = y1_sum^2 + (m*n1)^2 - 2*m*n1*y1_sum; %sum( y1 - m * vec1 )^2;
        gamma  = alpha*n1 - beta;
        
        %TODO! beta >= alpha not always true
        if beta > alpha
            rho      = (beta - alpha)/gamma;           
        else
            rho      = 0;
        end
        
        f_set(k) =      n0 * log( alpha0 ) ...
                   +    n1 * log( alpha + rho*gamma ) ...
                   -(n1-1) * log( 1 + rho * n1 );
               
        k = k + 1;
    end
    
    %Find minimum
    [~, idx_min] = min(f_set);
    mu           = m_set(idx_min);
    
    
    %TEMP: plot
    %disp(mu)
    %figure
    %plot(m_set, f_set, 'LineWidth', 1.2), grid on, hold on
    %pause
    %close


end


%% Exit
alpha  = y1_sq    + (mu^2*n1) - 2*mu*y1_sum;    %norm(y1 - m * vec1 )^2;
beta   = y1_sum^2 + (mu*n1)^2 - 2*mu*n1*y1_sum; %sum( y1 - m * vec1 )^2;
gamma  = alpha*n1 - beta;

if beta > alpha
    rho = (beta - alpha)/gamma;
else
    rho = 0;
end


%TEMP: plot
%plot(mu,f_set(idx_min), '*')
%pause


end

