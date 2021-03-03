%% Loading Data and initialization
clear all
close all
clc
%% Adding every files in the path

addpath(genpath(pwd))

dataset = "urban";
options.display = true;
options.inneriter = 10;
options.delta_iter = 1e-2;
options.save = false;

methods = {"titan\_minvol","minvol"};
symbols = {'-o','-s'};

[X,r,maxtime] = dataset_loader(dataset);
options.maxtime = maxtime;

%%
nX = norm(X,'fro')^2;
[m,n] = size(X);
options.lambda=0.001;

% Initialization

% [K,H] = SNPA(X,r);
% W = full(X(:,K));
% if length(K) < r
%     warning('SNPA recovered less than r basis vectors.');
%     warning('The data poins have less than r vertices.');
%     r = length(K);
%     fprintf('The new value of r is %2.0d.\n',r);
% end

W = rand(m,r);
H = FGMfcnls(X,W,[],100);   

scale = sum(W);
W = W./scale;
H = H.*scale';

options.W = W;
options.H = H;

k = 0;
methods_computed = {};


%% titanminvol
% Without extrapolation
if ~isempty(find([methods{:}] == "noextratitan\_minvol"))
    k = k+1;
    disp("computing TITANized-minvol without extrapolation...")
    options.inertial = false;
    [W,H,e,en,ep,etx,lambda] = titanminvol(X,r,options);
    es{k} = e;
    ens{k} = en;
    eps{k} = ep;
    etxs{k} = etx;
    methods_computed{k} = "noextratitan\_minvol";
end

% With extrapolation
if ~isempty(find([methods{:}] == "titan\_minvol"))
    k = k+1;
    disp("computing TITANized-minvol...")
    options.inertial = true;
    [W1,H1,e,en,ep,etx,lambda] = titanminvol(X,r,options);
    es{k} = e;
    ens{k} = en;
    eps{k} = ep;
    etxs{k} = etx;
    methods_computed{k} = "titan\_minvol";
end

%% halsiminvolNMF
% Without extrapolation
if ~isempty(find([methods{:}] == "hals"))
    k = k+1;
    disp("computing hals...")
    options.inertial = false;
    [W,H,e,en,ep,etx,lambda] = halsiminvolNMF(X,r,options);
    es{k} = e;
    ens{k} = en;
    eps{k} = ep;
    etxs{k} = etx;
    methods_computed{k} = "hals";
end

% With extrapolation
if ~isempty(find([methods{:}] == "hals\_extra"))
    k = k+1;
    disp("computing hals with TITAN extra...")
    options.inertial = true;
    [W,H,e,en,ep,etx,lambda] = halsiminvolNMF(X,r,options); 
    es{k} = e;
    ens{k} = en;
    eps{k} = ep;
    etxs{k} = etx;
    methods_computed{k} = "hals\_extra";
end

%% minvolNMF
if ~isempty(find([methods{:}] == "minvol"))
    k = k+1;
    disp("computing minvolNMF...")
    [W2,H2,e,en,ep,etx,lambda] = minvolNMF(X,r,options);
    es{k} = e;
    ens{k} = en;
    eps{k} = ep;
    etxs{k} = etx;
    methods_computed{k} = "minvol";
end

%% Displaying results
es_min = min([es{:}]);
ens_min = min([ens{:}]);
eps_min = min([eps{:}]);

set(0, 'DefaultAxesFontSize', 18);
set(0, 'DefaultLineLineWidth', 2);

figure; 
for kk = 1:length(methods)
    semilogy(etxs{kk},(es{kk}-es_min)./nX,symbols{kk}); hold on
end
ylim([-inf inf])
xlim([-inf inf])
legend(methods_computed{:})
xlabel('time (s)'); 
title('(||X-WH|| + lambda*logdet - min)/||X|| over time')
grid on
set(gca,'FontUnits','points','FontSize',18,'FontName','Times')

figure; 
for kk = 1:length(methods)
    semilogy(etxs{kk},(ens{kk}-ens_min)./nX,symbols{kk}); hold on
end
ylim([-inf inf])
xlim([-inf inf])
legend(methods_computed{:})
xlabel('time (s)'); 
title('(||X-WH|| - min)/||X|| over time')
grid on
set(gca,'FontUnits','points','FontSize',18,'FontName','Times')

figure; 
for kk = 1:length(methods)
    semilogy(etxs{kk},(eps{kk}-eps_min)./nX,symbols{kk}); hold on
end
ylim([-inf inf])
xlim([-inf inf])
legend(methods_computed{:})
xlabel('time (s)'); 
title('(lambda*logdet - min)/||X|| over time')
grid on
set(gca,'FontUnits','points','FontSize',18,'FontName','Times')

%% Homogenizing time sampling
if options.save 
    minmaxtime = options.maxtime;
    for kk = 1:k
        minmaxtime = min(minmaxtime,max(etxs{kk}));
    end
    new_time_sampling = 0:0.1:minmaxtime;

    for kk = 1:k
        sampled_es{kk} = interp1(etxs{kk},(es{kk}-es_min)./nX,new_time_sampling)';
        sampled_ens{kk} = interp1(etxs{kk},(ens{kk}-ens_min)./nX,new_time_sampling)';
        sampled_eps{kk} = interp1(etxs{kk},(eps{kk}-eps_min)./nX,new_time_sampling)';
        todat = cat(2,new_time_sampling',sampled_es{kk});
        save(pwd+"\saved_data\"+dataset+num2str(kk)+".dat", 'todat', '-ascii')
    end
end
