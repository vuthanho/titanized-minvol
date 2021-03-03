
function [W,H,e,err1,err2,etx,lambda] = titanminvol(X,r,options)

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
if ~isfield(options,'inertial')
    inertial = false;
else
    inertial = options.inertial;
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
HHt=H*H';
XHt=X*H';
delta_eyer = delta*eye(r);
err1(1) = max(0,normX2-2*sum(sum(WtX.*H))+sum(sum( WtW.*HHt)));
err2(1) = log( det (WtW  + delta_eyer ) );
lambda = options.lambda * max(1e-6,err1) / abs( err2 );
P = inv( WtW + delta*eye(r) );
LpW = norm(HHt+lambda*P);
LpH = norm(WtW);
e(1) =  0.5*err1 + 0.5*lambda * err2 ;

% Main loop 
k=1;
etx = 0;
ety0 = 0;
Wold = W;
Hold = H;

% extrapolation sequence
alpha1=1;
alpha2=1;
tic
while etx(k) < options.maxtime
    k = k+1;
    % *** Update W ***   
    iter = 1; % inner iter counter
    % Stop if ||W^{k}-W^{k+1}||_F <= delta * ||W^{0}-W^{1}||_F
    eps0 = 0; eps = 1;
	Woldold = W;
    while iter <= inneriter && eps >= delta_iter*eps0
        alpha0 = alpha1;
        alpha1 = 0.5*(1+sqrt(1+4*alpha0^2));
        P = inv( WtW + delta_eyer );
        LW = norm(HHt+lambda*P); % New Lipschitz constant for W
        if inertial
			beta = min((alpha0-1)/alpha1 , 0.9999*sqrt(LpW/LW));
            Wextra = W + beta*(W-Wold);
            Wold=W;
            W = simplexProj(Wextra + 1/LW*(XHt - Wextra*(HHt+lambda*P)), 1e-16);
        else
            W = simplexProj(W + 1/LW*(XHt - W*(HHt+lambda*P)), 1e-16);
        end
        if iter == 1
            eps0 = norm(W-Woldold,'fro'); 
        end
        WtW = W'*W; % Pre-computing to save computation time
        LpW=LW;
        eps = norm(W-Woldold,'fro');
        iter = iter + 1;
    end
    LH = norm(WtW); % New Lipschitz constant for H
    WtX = W'*X; % Pre-computing to save computation time
   
    
    % *** Update H ***
    iter = 1; % inner iter counter
    % Stop if ||H^{k}-H^{k+1}||_F <= delta * ||H^{0}-H^{1}||_F
    eps0 = 0; eps = 1;
	Holdold = H;
    while iter <= inneriter && eps >= delta_iter*eps0
        alpha0 = alpha2;
        alpha2 = 0.5*(1+sqrt(1+4*alpha0^2));
					   
        if inertial
			beta = min((alpha0-1)/alpha2 , 0.9999*sqrt(LpH/LH));
            Hextra = H + beta*(H-Hold);
            Hold=H;
            H = max(1e-16,Hextra + 1/LH*(WtX - WtW*Hextra));
        else
            H = max(1e-16,H + 1/LH*(WtX - WtW*H));
        end
        if iter == 1
            eps0 = norm(H-Holdold,'fro'); 
        end
        eps = norm(H-Holdold,'fro');
        LpH=LH; 
        iter = iter + 1;
    end
    HHt = H*H'; % Pre-computing to save computation time
    XHt = X*H'; % Pre-computing to save computation time
    
    % *** Saving the error without taking into account the error
    % computation time *** %
    etz = toc;
    err1(k) = max(0, normX2 - 2*sum(sum( WtX.*H ) )  + sum(sum( WtW.*HHt ) ) );
    err2(k) = log( det ( WtW + delta_eyer) );
    e(k) = 0.5*err1(k) + 0.5*lambda * err2(k);
    ety1 = toc;
    etx(k) = etx(k-1) + etz - ety0; 
    ety0 = ety1;
end