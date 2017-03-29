clear all;
cache = LRU_Cache(10);

temp = cache.head;
while ~isempty(temp.Next)
    temp = temp.Next;
end

last = temp;
%%

clear all;clc
cache = LRU_Cache(10);
% printCache(cache);
emulate(cache,2);
% printCache(cache);

%%
clear all; clc;

ZipfExponent = 0.8;
CatalogSize = 1000;

NumRequests = 100000;

weights = 1:CatalogSize;
weights = weights.^-ZipfExponent;

CacheSizeSamples = [100 200:200:CatalogSize];


PhitVecArray = zeros(11,1000);
cacheSizeIndex = 0;
for run = 1:1
    Phit = [];
    ContentRequests = randsample(CatalogSize,NumRequests,'true',weights);
    for CacheSize = CacheSizeSamples
        cacheSizeIndex = cacheSizeIndex +1;
        CacheLRU = LRU_Cache(CacheSize,CatalogSize); % Cache is initially populated with most popular contents
        
        HitCount = 0;
        for index = 1:NumRequests
            HitCount = HitCount + emulate(CacheLRU,ContentRequests(index));
            
        end
        PhitVec = collectStats(CacheLRU); % 1xCatalogSize
        PhitVecArray(cacheSizeIndex,:) = PhitVec;
        Phit = [Phit HitCount/NumRequests]
        [run CacheSize sum(PhitVec.*weights)/sum(weights)]
    end
end

save('Simulation.mat')

%%
close;
figure; grid on;hold on;
plot(100:100:1000,PhitArray(1,:));
plot(100:100:1000,PhitArray(2,:),'r');
plot(100:100:1000,PhitArray(3,:),'k');
plot(100:100:1000,PhitArray(4,:),'g');
plot(100:100:1000,PhitArray(5,:),'c');
ylabel('P_{hit}');
xlabel('Cache Size; Catalog Size = 1000');

%%

clear all; clc;

ZipfExponent = 0.8;
CatalogSize = 1000000;

NumRequests = 100000;

weights = 1:CatalogSize;
weights = weights.^-ZipfExponent;

weights = weights/sum(weights);

N = 10:10:100000000;
s = zeros(1,size(N,2));
for index = 1:size(N,2)
    s(index) = (sum((1-weights).^N(index)));
end

CacheSize = 1000:1000:900000;
Tc = zeros(1,size(CacheSize,2));
Ph = zeros(1,size(CacheSize,2));

for index = 1:size(CacheSize,2)
    for ind = 1:size(s,2)
        if s(ind) < CatalogSize - CacheSize(index)
            Tc(index) = (ind);
            Ph(index) = sum(weights.*(1-(1-weights).^Tc(index)));
            break;
        end
    end
end

%%
clear all; 

ZipfExponent = 0.7;
CatalogSize = 1000;

weights = 1:CatalogSize;
weights = weights.^-ZipfExponent;

weights = weights/sum(weights);

Tc = 10:10:1000;

C = zeros(1,size(Tc,2));
Ph = zeros(1,size(Tc,2));

for index = 1:size(Tc,2)
    temp = 1-exp(-weights.*Tc(index));
    C(index)= sum(temp);
    Ph(index) = sum(weights.*temp);
end
plot(C,Ph);
grid on;










