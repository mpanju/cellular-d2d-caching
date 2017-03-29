clear all;
clc;

ZipfExponent = 0.8;
CatalogSize = 1000;

NumRequests = 10000;

weights = 1:CatalogSize;
weights = weights.^-ZipfExponent;

CacheSizeSamples = 100:100:1000;

ContentRequested = zeros(length(CacheSizeSamples),size(NumRequests,2));
ContentInLocalCache = zeros(length(CacheSizeSamples),size(NumRequests,2));
ContentInOtherCache = zeros(length(CacheSizeSamples),size(NumRequests,2));


i =1;
for CacheSize = CacheSizeSamples/2;

    tic
    
    cache(1) = LRU_Cache(CacheSize,CatalogSize);
    cache(2) = LRU_Cache(2*CacheSize,CatalogSize);
    
    ContentRequests = randsample(CatalogSize,NumRequests,'true',weights);
    
    HitCountSelf = 0;
    HitCountOther = 0;
    
    for index = 1:NumRequests

        CID = ContentRequests(index);
        
        ContentRequested(i,index) = CID;
       
        
        hit =  emulate(cache(1),CID);
        ContentInLocalCache(1,index) = hit;

        if(hit == 1)
            HitCountSelf = HitCountSelf+1;
            ContentInOtherCache(i,index) = isContentPresent(cache(2),CID);
        else
            h = emulate(cache(2),CID);
            if (h == 1)
               HitCountOther = HitCountOther +1; 
            end
            ContentInOtherCache(i,index) = h;
        end
        
    end
    
    stats(i).HitCountSelf = HitCountSelf;
    stats(i).HitCountOther = HitCountOther;
    stats(i).CacheSize = CacheSize;
    stats(i).statsHitVec_1 = collectStats(cache(1));
    stats(i).statsHitVec_2 = collectStats(cache(2));
    stats(i).RequestCount = NumRequests;
    fprintf('Simulation for cache size %d finished.\n',CacheSize);
    i = i+1;
    toc
end

%% Plot Simulation
close all;
SelfHitCounter = zeros(1,length(stats));
OtherHitCounter = zeros(1,length(stats));
TotalRequests = zeros(1,length(stats));
CacheSizeSamples = zeros(1,size(stats,2));

% StatsHitVec_1 = zeros(size(stats,2),CatalogSize);
% StatsHitVec_2 = zeros(size(stats,2),CatalogSize);
for j = 1:size(stats,2)
    SelfHitCounter(j) = stats(j).HitCountSelf;
    OtherHitCounter(j) = stats(j).HitCountOther;
    CacheSizeSamples(j) = stats(j).CacheSize;
    TotalRequests(j) = stats(j).RequestCount;
end

plot(CacheSizeSamples,SelfHitCounter./TotalRequests,'x'); hold on;
plot(CacheSizeSamples,OtherHitCounter./TotalRequests,'+r');
plot(CacheSizeSamples,OtherHitCounter./(TotalRequests-SelfHitCounter),'+k');


%% Plot theoretical


ZipfExponent = 0.8;
CatalogSize = 1000;

weights = 1:CatalogSize;
weights = weights.^-ZipfExponent;
weights = weights/sum(weights);

tc = 10:10:10000;
C = zeros(1,size(tc,2));
PHm = zeros(size(tc,2),size(weights,2));
Phit = zeros(1,size(tc,2));
Phit2 = zeros(1,size(tc,2));
for t = 1:size(tc,2)
    for index = 1:size(weights,2)
        pin = 0; pin_ = 1;
        while abs(pin_ - pin) > 1e-6
            pin_ = pin;
            pin = (1-exp(- tc(t).*weights(index).*(1+(1-pin))))./(1+(1-pin).*exp(-tc(t).*weights(index).*(1+(1-pin))));
        end
        PHm(t,index) = pin;
    end
    C(t)  = sum(PHm(t,:));
    Phit(t) = sum( PHm(t,:).*weights );
    Phit2(t) = sum(  PHm(t,:).*(1-PHm(t,:)).*weights );
end
plot(C,Phit); hold on;
plot(C,Phit2,'r');

%%

plot(C,Phit + Phit.*(1-Phit));hold on;

plot(CacheSizeSamples,(OtherHitCounter(1,:) + SelfHitCounter(1,:))./TotalRequests(1,:),'x'); hold on;
plot(CacheSizeSamples,(OtherHitCounter(2,:) +SelfHitCounter(2,:))./TotalRequests(2,:),'+r');

%%

plot(OtherHitCounter(1,:)./(TotalRequests(1,:)-SelfHitCounter(1,:))); hold on;
plot((SelfHitCounter(1,:) +OtherHitCounter(1,:))./(TotalRequests(1,:)),'r'); hold on;
plot(SelfHitCounter(1,:)./(TotalRequests(1,:)),'k'); hold on;
%%
hold on;

plot( CacheSizeSamples, (SelfHitCounter(1,:)+OtherHitCounter(2,:)) ./ (TotalRequests(1,:)+TotalRequests(2,:)-SelfHitCounter(2,:)) ,'+r');
plot( CacheSizeSamples, (SelfHitCounter(2,:)+OtherHitCounter(1,:)) ./ (TotalRequests(2,:)+TotalRequests(1,:)-SelfHitCounter(1,:)) ,'+r');

plot(CacheSizeSamples, ( OtherHitCounter(1,:)./ (TotalRequests(1,:) - SelfHitCounter(1,:) )) ,'x');
plot(CacheSizeSamples, ( OtherHitCounter(2,:)./ (TotalRequests(2,:) - SelfHitCounter(2,:) )) ,'x');

%%
plot(CacheSizeSamples,  SelfHitCounter(1,:) ./ TotalRequests(1,:) ,'or'); hold on;
plot(CacheSizeSamples,  SelfHitCounter(2,:) ./ TotalRequests(2,:) ,'or');


plot(CacheSizeSamples,  (SelfHitCounter(1,:) + OtherHitCounter(1,:) )./ TotalRequests(1,:) ,'+k');
plot(CacheSizeSamples,  (SelfHitCounter(2,:)  + OtherHitCounter(2,:) )./ TotalRequests(2,:) ,'+k');



plot(CacheSizeSamples,  (OtherHitCounter(1,:) )./ TotalRequests(1,:) ,'.');
plot(CacheSizeSamples,  (OtherHitCounter(2,:) )./ TotalRequests(2,:) ,'.');





