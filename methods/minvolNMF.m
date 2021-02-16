% Solving approximately
% 
% min_{W,H} ||M-WH||_F^2 + lambda' * logdet( W^TW + delta I) 
% 
% where W >= 0, H >= 0, ||H(:,j)||_2 <= 1 for all j. 
% 
% This is sovled by optimizing alternatively over W and H, 
% using a projected fast gradient method; see 
% Minimum-Volume Rank-Deficient Nonnegative Matrix Factorizations, 
% Valentin Leplat, Andersen M.S. Ang, Nicolas Gillis, 2018.
% 
% ****** Input ****** 
% X : m-by-n matrix to factorize 
% r : factorization rank r 
% options: 
%   - lambda will be used to set up lambda' (default: 0.1) 
%   - delta (default: 0.1) 
%   - maxtime (default: 100) 
%   - (W,H): initialization (default: SNPA) 
%   - display: =1 displays the iteration count, = 0 no display
%               (default: 0).  
% 
% ****** Outut ****** 
% (W,H) : low-rank approximation W*H of X, with W >= 0 of small volume 
%         (that is, logdet(W^TW + delta*I) 
%         and H(:,j) >= 0 and ||H(:,j)||_1 <= 1. 
% e     : evolution of the error 
%         ||M-WH||_F^2 + lambda * logdet( W^TW + delta I) 
% err1  = evolution of ||M-WH||_F^2
% err2  = evolution of logdet( W^TW + delta I) 
% lambda: 

function [W,H,e,err1,err2,etx,lambda] = minvolNMF(X,r,options)

if nargin <= 2
    options = [];
end
if ~isfield(options,'lambda')
    options.lambda=0.1;
end
if ~isfield(options,'delta')
    delta=0.1;
else
    delta=options.delta;
end
if ~isfield(options,'maxtime')
    options.maxtime=100;
end
if isfield(options,'W') && isfield(options,'H') 
    W = options.W;
    H = options.H;
else
    [K,H] = SNPA(X,r);
    W = X(:,K);
    if length(K) < r
        warning('SNPA recovered less than r basis vectors.');
        warning('The data poins have less than r vertices.');
        r = length(K);
        fprintf('The new value of r is %2.0d.\n',r);
    end
end
if isfield(options,'delta_iter')
    delta_iter = options.delta_iter;
else
    delta_iter = 1e-6;
end
if isfield(options,'inneriter') % number of updates of W and H, before the other is updated
    inneriter = options.inneriter;
else
    inneriter = 100;
end

% Initializations
normX2 = norm(X,'fro')^2;
WtW = W'*W;
WtX = W'*X;
err1(1) = max(0,normX2-2*sum(sum(WtX.*H))+sum(sum( WtW.*(H*H'))));
err2(1) = log( det (WtW  + delta*eye(r) ) );
lambda = options.lambda * max(1e-6,err1) / abs( err2 );
e(1) =  0.5*err1 + 0.5*lambda * err2 ;

% Main loop 
i=1;
etx = 0;
ety0 = 0;
tic
while etx(i) < options.maxtime
    i=i+1;
    % *** Update W ***
    XHt = X*H';
    HHt = H*H';
    Y = inv( ( W'*W + delta*eye(r) ) );
    A = lambda*Y + HHt;
    W = FGMqpnonneg(A,XHt,W,inneriter,delta_iter); 
    % *** Update H ***
    [H,WtW,WtX] = FGMfcnls(X,W,H,inneriter,delta_iter);
    
    % *** Saving the error without taking into account the error
    % computation time *** %
    etz = toc;
    err1(i) = max(0, normX2 - 2*sum(sum( WtX.*H ) )  + sum(sum( WtW.*(H*H') ) ) );
    err2(i) = log( det ( WtW + delta*eye(r) ) );
    e(i) = 0.5*err1(i) + 0.5*lambda * err2(i);
    ety1 = toc;
    etx(i) = etx(i-1) + etz - ety0; 
    ety0 = ety1;
end