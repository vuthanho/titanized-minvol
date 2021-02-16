function [p,permutedClusterInds] = clustPurity(clusterInds,labelInds)

% preconditions: 1) #clusters == #labels; 2) sum_ij(Xi,Cj) = sum_ij(Xi,Lj)
% cluster indices/labels = [1:numClusters]

numEls = length(clusterInds);
numClusters = length(unique(clusterInds));
numLabels = length(unique(labelInds));

% if ~(numClusters == numLabels)
%     error('Different number of clusters and labels');
% end

if ~(length(clusterInds)==length(labelInds))
    error('There must be the same number of elements in both vectors');
end

profitMatrix = zeros(numClusters,numLabels);
for i = 1:numClusters
    for j = 1:numLabels
        profitMatrix(i,j) = nnz(clusterInds==i & labelInds==j);        
    end
end

[bestMatchings,cost] = munkres(-profitMatrix);
permutedClusterInds = zeros(size(clusterInds));
for i=1:numClusters
    permutedClusterInds(clusterInds==i) = bestMatchings(i);    
end
p = abs(cost)/numEls;


