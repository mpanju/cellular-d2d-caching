clear all;
clc;

ZipfExponent = 0.8;
CatalogSize = 1000;

NumRequests = 10000;

weights = 1:CatalogSize;
weights = weights.^-ZipfExponent;

CacheSizeSamples = [10 50 100:100:1000];

ContentRequested = zeros(length(CacheSizeSamples),size(NumRequests,2));
ContentInLocalCache = zeros(length(CacheSizeSamples),size(NumRequests,2));
ContentInOtherCache = zeros(length(CacheSizeSamples),size(NumRequests,2));
WhichCache = zeros(length(CacheSizeSamples),size(NumRequests,2));



i =1;
for CacheSize = CacheSizeSamples;

    tic
    
    cache(1) = LRU_Cache(CacheSize,CatalogSize);
    cache(2) = LRU_Cache(CacheSize,CatalogSize);
    
    ContentRequests = randsample(CatalogSize,NumRequests,'true',weights);
    
    HitCountSelf = [ 0 0];
    HitCountOther = [0 0];
    RequestCount = [0 0];
    
    for index = 1:NumRequests
        c = randsample(2,1);
        d = invertCacheIndex(c);
       
        
        
        RequestCount(c) = RequestCount(c) +1;
        CID = ContentRequests(index);
        
        ContentRequested(i,index) = CID;
        ContentInLocalCache(i,index) = isContentPresent(cache(c),CID);
        ContentInOtherCache(i,index) = isContentPresent(cache(d),CID);
        WhichCache(i,index) = c;
        
        
        hit =  emulate(cache(c),ContentRequests(index));
        
        if(hit == 1)
            HitCountSelf(c) = HitCountSelf(c)+1;
        else
%             if ( isContentPresent(cache(d),CID) == 1 )
 
                h = emulate(cache(d),CID);
                if( h == 1)
                    HitCountOther(c) = HitCountOther(c) +1;
%                     fprintf('Somthing Wrong!\n');
                end
%             end
        end
        
    end
    
    stats(i).HitCountSelf = HitCountSelf;
    stats(i).HitCountOther = HitCountOther;
    stats(i).RequestCount = RequestCount;
    stats(i).CacheSize = CacheSize;
    stats(i).statsHitVec_1 = collectStats(cache(1));
    stats(i).statsHitVec_2 = collectStats(cache(2));
    fprintf('Simulation for cache size %d finished.\n',CacheSize);
    i = i+1;
    toc
end

% Plot Simulation
close all;
SelfHitCounter = zeros(2,size(stats,2));
OtherHitCounter = zeros(2,size(stats,2));
TotalRequests = zeros(2,size(stats,2));
CacheSizeSamples = zeros(1,size(stats,2));
StatsHitVec_1 = zeros(size(stats,2),CatalogSize);
StatsHitVec_2 = zeros(size(stats,2),CatalogSize);
for i = 1:2
    for j = 1:size(stats,2)
        SelfHitCounter(i,j) = stats(j).HitCountSelf(i);
        OtherHitCounter(i,j) = stats(j).HitCountOther(i);
        TotalRequests(i,j) = stats(j).RequestCount(i);
        CacheSizeSamples(j) = stats(j).CacheSize;
    end
end
for i = 1:size(stats,2)
   StatsHitVec_1(i,:) = stats(1).statsHitVec_1; 
   StatsHitVec_2(i,:) = stats(2).statsHitVec_2;
end

%%

lch1Vec = zeros(1,length(CacheSizeSamples));
lch2Vec = zeros(1,length(CacheSizeSamples));

och1Vec = zeros(1,length(CacheSizeSamples));
och2Vec = zeros(1,length(CacheSizeSamples));

rlVec1 =  zeros(1,length(CacheSizeSamples));
rlVec2 =  zeros(1,length(CacheSizeSamples));

mshVec1 = zeros(1,length(CacheSizeSamples));
mshVec2 = zeros(1,length(CacheSizeSamples));
i = 1;
for ci = CacheSizeSamples
   cl1 = find(WhichCache(i,:) == 1);
   cl2 = find(WhichCache(i,:) == 2);
   
   rlVec1(i) = length(cl1);
   rlVec2(i) = length(cl2);

   
   lch1Vec(i) = sum(ContentInLocalCache(i,cl1))/length(cl1);
   lch2Vec(i) = sum(ContentInLocalCache(i,cl2))/length(cl2);
     
   och1Vec(i) = sum(ContentInOtherCache(i,cl1))/length(cl1);
   och2Vec(i) = sum(ContentInOtherCache(i,cl2))/length(cl2);
  
   j = 1;
   missCount = 0; missHitCount = 0;
   for ji = cl1
       if ContentInLocalCache(i,ji) == 0
           missCount = missCount+1;
           if ContentInOtherCache(i,ji) == 1;
               missHitCount = missHitCount +1;
           end
       end
       j = j+1;
   end
   mshVec1(i) = missHitCount/missCount;
   
   
   j = 1;
   missCount = 0; missHitCount = 0;
   for ji = cl2
       if ContentInLocalCache(i,ji) == 0
           missCount = missCount+1;
           if ContentInOtherCache(i,ji) == 1;
               missHitCount = missHitCount +1;
           end
       end
       j = j+1;
   end
   mshVec2(i) = missHitCount/missCount;
   
   
   i = i+1;
end

hold on;

plot(CacheSizeSamples,lch1Vec);
plot(CacheSizeSamples,lch2Vec);


plot(CacheSizeSamples,och1Vec,'r');
plot(CacheSizeSamples,och2Vec,'r');


plot(CacheSizeSamples,mshVec1,'k');
plot(CacheSizeSamples,mshVec2,'k');

%%
plot(CacheSizeSamples,SelfHitCounter(1,:)./TotalRequests(1,:),'x'); hold on;
plot(CacheSizeSamples,SelfHitCounter(2,:)./TotalRequests(2,:),'+r');


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
for t = 1:size(tc,2)
    for index = 1:size(weights,2)
        pin = 0; pin_ = 1;
        while abs(pin_ - pin) > 1e-6
            pin_ = pin;
            pin = (1-exp(- tc(t).*weights(index).*(1+(1-pin))));
        end
        PHm(t,index) = pin;
    end
    C(t)  = sum(PHm(t,:));
    Phit(t) = sum( PHm(t,:).*weights );
end
plot(C,Phit);
%%
hold on;
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





