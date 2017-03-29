%%
 clear all; clc;
% AdjMatrix = [0 1 1 1 1;
%     1 0 1 1 1;
%     1 1 0 1 1;
%     1 1 1 0 1;
%     1 1 1 1 0];

AdjMatrix = [ 
     0     1     0     0     0;
     1     0     1     1     1;
     0     1     0     1     0;
     0     1     1     0     0;
     0     1     0     0     0];



N = 5;
M = 1000;

CacheSizeSamples = 100:100:500;

% CacheSizeSamples = 100;

PhitVec = zeros(N,length(CacheSizeSamples));

ZipfExponent = 0.8;
lambdaArray = 1:M;
lambdaArray = lambdaArray.^-ZipfExponent;
lambdaArray = lambdaArray/sum(lambdaArray);

lambdaArray = repmat(lambdaArray,N,1);

%%

for index = 1:length(CacheSizeSamples)
    fprintf('Iteration started for C = %d\n', CacheSizeSamples(index));

    TcVec = 10*ones(1,N);
    PnmArray = ones(N,M) * 0.5;
    TcVec_ = zeros(1,N);
    PnmArray_ = zeros(N,M);

    
    while norm(PnmArray-PnmArray_) > 0.1 || norm(TcVec_ - TcVec) > 0.001
        PnmArray_ = PnmArray;
        TcVec_ = TcVec;
        PnmArray = computePnmArray( AdjMatrix, TcVec, lambdaArray);
        TcVec = computeTc(AdjMatrix, lambdaArray, PnmArray, CacheSizeSamples(index));
    end
    for node = 1:N
        PhitVec(node,index) = sum(PnmArray(node,:).*lambdaArray(node,:))/sum(lambdaArray(node,:));
    end
end
%%
 hold on;
for index = 1:N
    plot(CacheSizeSamples,PhitVec(index,:));
    
end

%%

ZipfExponent = 0.8;
CatalogSize = 100;

Ncache = 5;
AdjMatrix = [0 1 1 1 1;
    1 0 1 1 1;
    1 1 0 1 1;
    1 1 1 0 1;
    1 1 1 1 0];


node(1).neighbors = [2 3 4 5];
node(2).neighbors = [1 2 3 4];
node(3).neighbors = [1 2 4 5];
node(4).neighbors = [1 2 3 5];
node(5).neighbors = [1 2 3 4];

NumRequests = Ncache*CatalogSize*5;
%%
ZipfExponent = 0.8;
CatalogSize = 100;

Ncache = 1;
AdjMatrix = 0;

node(1).neighbors = [];
NumRequests = Ncache*CatalogSize*5;
