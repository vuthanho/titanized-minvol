function [X,r,maxtime] = dataset_loader(dataset)

% reviews

if dataset == "reviews"
    load('reviews.mat')
    X = dtm';
    clear dtm
    r = 5;
    maxtime = 30;
end

% sports

if dataset == "sports"
    load('sports.mat')
    X = dtm';
    clear dtm
    r = 7;
    maxtime = 30;
end

% 20news

if dataset == "20news"
    load('20news\test.data')
    X = sparse(test(:,2),test(:,1),test(:,3));
    clear('test')
    r=20;
    maxtime = 300;
end

% Terrain 

if dataset == "terrain"
    load('TerrainClean.mat','A')
    X = A';
    clear A
    r = 5;
    maxtime = 60;
end

% Pavia University

if dataset == "paviau"
    load('PaviaU.mat')
    X = reshape(paviaU,[610*340,103])';
    r = 9;
    maxtime = 90;
end

% San Diego

if dataset == "sandiego"
    load('SanDiego.mat')
    X = A';
    r = 7;
    maxtime = 120;
end

% Indian Pines

if dataset == "indian"
    load('Indian_pines_corrected.mat')
    X = reshape(indian_pines_corrected,[145*145,200])';
    r = 16;
    maxtime = 30;
end

% Urban

if dataset == "urban"
    load('Urban.mat')
    X = A';
    clear A
    r = 6;
    maxtime = 60;
end

% Mnist

% if dataset == "mnist"
    % X = processImagesMNIST('train-images.idx3-ubyte');
    % [m,n] = size(X);
    % X = X(:,randperm(n));
    % N = 30000;
    % X = X(:,1:N);
    % fprintf('Number of images kept in the dataset: %6d ...\n',N);
    % r = 30;
    % maxtime = 60;
% end

% Classic text mining

% if dataset == "classic"
    % load('classic.mat')
    % X = dtm';
    % clear dtm
    % r = 4;
    % maxtime = 300;
% end

% Samson

% if dataset == "samson"
    % load('samson_1.mat')
    % X = V;
    % r = 3;
% end

% Jasper 

% if dataset == "jasper"
    % load('jasperRidge2_R198.mat')
    % X = Y;
    % r=4;
% end

% Yale

% if dataset == "yale"
    % load('YaleB_32x32.mat')
    % X = fea'./max(fea(:));
    % [~,n] = size(X);
    % X = X(:,randperm(n));
    % clear fea gnd
    % r = 49;
% end

% la1

% if dataset == "la1"
    % load('la1.mat')
    % X = dtm';
    % clear dtm
    % r = 6;
    % maxtime = 300;
% end

X = X/max(X(:));
end

