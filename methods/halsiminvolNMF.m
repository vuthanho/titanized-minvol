
function [W,H,e,err1,err2,etx,lambda] = halsiminvolNMF(X,r,options)

IndWi = [];
for i =1:r
    IndWi = cat(1,IndWi,cat(2,1:(i-1),(i+1):r));
end

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
err1(1) = max(0,normX2-2*sum(sum(WtX.*H))+sum(sum( WtW.*HHt)));
err2(1) = log( det (WtW  + delta*eye(r) ) );
lambda = options.lambda * max(1e-6,err1) / abs( err2 );
P = inv( WtW + delta*eye(r) );
LpW = diag(HHt+lambda*P);
LW = LpW;
LH = diag(WtW);
e(1) =  0.5*err1 + 0.5*lambda * err2 ;

alpha1=1;
% Main loop
k=1;
etx = 0;
ety0 = 0;
tic
while etx(k) < options.maxtime
    k=k+1;
    alpha0 = alpha1;
    alpha1 = 0.5*(1+sqrt(1+4*alpha0^2));
    
    % *** Update W ***
    LpH = LH; % Last Lipschitz constant for H
    iter = 1; % inner iter counter
    % Stop if ||W^{k}-W^{k+1}||_F <= delta * ||W^{0}-W^{1}||_F
    eps0 = 0; eps = 1;
    Wold = W;
    while iter <= inneriter && eps >= delta_iter*eps0
        Wp = W; 
        for i = 1:r % For each column of W
            C = chol(inv( W'*W + delta*eye(r) ),'lower');
            LWi = (HHt(i,i)+lambda*C(i,:)*C(i,:)');
            if inertial
                Woldi = Wold(:,i);
                Wi = W(:,i);
                beta = min((alpha0-1)/alpha1*0.5 , 0.9999*sqrt(LpW(i)/LWi));
                Wextra = beta*(Wi-Woldi);
                W(:,i) = simplexProj(Wextra + W(:,i) + (XHt(:,i)-W*HHt(:,i)...
                -lambda*W(:,IndWi(i,:))*C(IndWi(i,:),:)*C(i,:)')/LWi,1e-16);
            else
                W(:,i) = simplexProj( W(:,i) + (XHt(:,i)-W*HHt(:,i)...
                -lambda*W(:,IndWi(i,:))*C(IndWi(i,:),:)*C(i,:)')/LWi,1e-16);
            end
        end
        if iter == 1
            eps0 = norm(W-Wp,'fro'); 
        end
        eps = norm(W-Wp,'fro');
        Wold = Wp;
        iter = iter + 1;
    end
    WtW = W'*W;
    LH = diag(WtW); % New Lipschitz constant for H
    WtX = W'*X;
    P = inv( WtW + delta*eye(r) );
    
    
    % *** Update H ***
    LpW = diag(HHt) + lambda*diag(P); % Last Lipschitz constant for W
    iter = 1; % inner iter counter
    % Stop if ||H^{k}-H^{k+1}||_F <= delta * ||H^{0}-H^{1}||_F
    eps0 = 0; eps = 1;
    Hold = H;
    while iter <= inneriter && eps >= delta_iter*eps0
        Hp = H;
        for i = 1:r % For each row of H
            if inertial
                Holdi = Hold(i,:);
                Hi = H(i,:);
                beta = min((alpha0-1)/alpha1*0.5 , 0.9999*sqrt(LpH(i)/LH(i)));
                Hextra = beta*(Hi-Holdi);
                H(i,:) = max(Hextra + H(i,:) + (WtX(i,:)-WtW(i,:)*H)/LH(i),1e-16);
            else
                H(i,:) = max(0,H(i,:) + (WtX(i,:)-WtW(i,:)*H)/LH(i));
            end
        end
        if iter == 1
            eps0 = norm(H-Hp,'fro'); 
        end
        eps = norm(H-Hp,'fro');
        Hold = Hp;
        iter = iter + 1;
    end
    HHt = H*H';
    XHt = X*H';
    
    % *** Saving the error without taking into account the error
    % computation time *** %
    etz = toc;
    err1(k) = max(0, normX2 - 2*sum(sum( WtX.*H ) )  + sum(sum( WtW.*HHt ) ) );
    err2(k) = log( det ( WtW + delta*eye(r) ) );
    e(k) = 0.5*err1(k) + 0.5*lambda * err2(k);
    ety1 = toc;
    etx(k) = etx(k-1) + etz - ety0; 
    ety0 = ety1;
end