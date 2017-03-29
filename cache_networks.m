%% Setup
clear all;
clc;

ZipfExponent = 0.8;
CatalogSize = 1000;



Ncache = 5;
AdjMatrix = [0 1 0 0 0;
             1 0 1 1 1;
             0 1 0 1 0;
             0 1 1 0 0;
             0 1 0 0 0];



node(1).neighbors = [2];
node(2).neighbors = [1 3 4 5];
node(3).neighbors = [2 4];
node(4).neighbors = [2 3];
node(5).neighbors = 2;

NumRequests = Ncache*CatalogSize*5;

weights = 1:CatalogSize;
weights = weights.^-ZipfExponent;

% weights = ones(1,CatalogSize);

CacheSizeSamples = [10:10:100 200:100:1000];

%% Run simulation

i =1;
for CacheSize = CacheSizeSamples;

    tic
    
    for k = 1:Ncache
        cache(k) = LRU_Cache(CacheSize,CatalogSize);
    end
    
    for index = 1:Ncache
        node(index).nbrHitCount = zeros(1,length( node(index).neighbors));
    end
    
    ContentRequests = randsample(CatalogSize,NumRequests,'true',weights);
    
    HitCountSelf = zeros(1,Ncache);
    HitCountOther = zeros(1,Ncache);
    RequestCount = zeros(1,Ncache);
    
    for index = 1:NumRequests
        c = randsample(Ncache,1);


        RequestCount(c) = RequestCount(c) +1;
        CID = ContentRequests(index);
        
        
        hit =  emulate(cache(c),ContentRequests(index));
        
        if(hit == 1)
            HitCountSelf(c) = HitCountSelf(c)+1;
        else
            nbrs = node(c).neighbors;
            probeResult = zeros(1,length(nbrs));
            
            for j = 1:length(nbrs)
                probeResult(j) = isContentPresent(cache(nbrs(j)),CID);
            end
                
            if (sum(probeResult) ~= 0)
                d = datasample ( nbrs(probeResult == 1), 1);
                node(c).nbrHitCount ( node(c).neighbors == d) = node(c).nbrHitCount ( node(c).neighbors == d) + 1;
                HitCountOther(c) = HitCountOther(c) +1;
                h = emulate(cache(d),CID);
                if( h ==0)
                    fprintf('Somthing Wrong!\n');
                    return;
                end
            end
        end
        
    end
    
    stats(i).HitCountSelf = HitCountSelf;
    stats(i).HitCountOther = HitCountOther;
    stats(i).RequestCount = RequestCount;
    stats(i).CacheSize = CacheSize;
    stats(i).statsHitVec_1 = collectStats(cache(1));
    stats(i).statsHitVec_2 = collectStats(cache(2));
%     
%     stats(i).hitCountStruct.node_1 = node(1).nbrHitCount;
%     stats(i).hitCountStruct.node_2 = node(2).nbrHitCount;
%     stats(i).hitCountStruct.node_3 = node(3).nbrHitCount;
    fprintf('Simulation for cache size %d finished.\n',CacheSize);
    i = i+1;
    toc
end

%% Plot Simulation
%  close all; close;
SelfHitCounter = zeros(Ncache,size(stats,2));
OtherHitCounter = zeros(Ncache,size(stats,2));
TotalRequests = zeros(Ncache,size(stats,2));
CacheSizeSamples = zeros(1,size(stats,2));

for i = 1:Ncache
    for j = 1:size(stats,2)
        SelfHitCounter(i,j) = stats(j).HitCountSelf(i);
        OtherHitCounter(i,j) = stats(j).HitCountOther(i);
        TotalRequests(i,j) = stats(j).RequestCount(i);
         CacheSizeSamples(j) = stats(j).CacheSize;
    end
end

hold on; grid on;

for k = 1:Ncache
    plot(CacheSizeSamples, SelfHitCounter(k,:)./TotalRequests(k,:), '+b');
    plot(CacheSizeSamples, (SelfHitCounter(k,:)+OtherHitCounter(k,:))./TotalRequests(k,:), '+k');
    
end

%% Plot theoretical

tc = 10.^(1:0.1:4);

C = zeros(Ncache,size(tc,2));
Pnm = zeros(Ncache,size(weights,2));
Phit = zeros(Ncache,size(tc,2));
PhitT = zeros(Ncache,size(tc,2));

% weights = ones(1,CatalogSize);
weights = weights/sum(weights);

for t = 1:size(tc,2)
    for index = 1:size(weights,2)
        pin_ = zeros(1,Ncache);
        pin = ones(1,Ncache);
        
        while norm(pin_ - pin) > 1e-6
            pin_ = pin;
            for k = 1:Ncache
                LAMBDA = sum(AdjMatrix(k,:).*pin_) ;
                
                EXP_TERM = exp(- tc(t).*weights(index).*(1+LAMBDA));
                
                pin(k) = (1-EXP_TERM)./(1+(LAMBDA).*EXP_TERM);
            end
        end
        
        Pnm(:,index) = pin;
    end

    Nd = repmat(sum(AdjMatrix,2),1,length(weights));

    
    for k =  1:Ncache
        C(k,t)  = sum(Pnm(k,:));
        Phit(k,t) = sum( Pnm(k,:).* weights);
        
        PROD_TERM = ones(1,length(weights));
        
        for j = 1:Ncache
           if(AdjMatrix(k,j) == 1)
              PROD_TERM = PROD_TERM .* (1 - Pnm(j,:));  
           end
        end
        
        PhitT(k,t) = sum( (1- (1-Pnm(k,:)).*PROD_TERM ).* weights );
    end
end

hold on; grid on;
for i = 1:Ncache
    plot(C(i,:),Phit(i,:),'b');
    plot(C(i,:),PhitT(i,:),'k');
end

