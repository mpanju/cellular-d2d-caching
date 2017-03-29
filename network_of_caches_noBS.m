function stats = network_of_caches(CacheSizeSamples, weights)


ZipfExponent = 0.8;
CatalogSize = 1000;



Ncache = 5;
AdjMatrix = [0 1 1 1 1;
             1 0 1 1 0;
             1 1 0 1 0;
             1 1 1 0 0;
             1 0 0 0 0];



node(1).neighbors = [2 3 4 5];
node(2).neighbors = [1 3 4 5];
node(3).neighbors = [1 2 4 5];
node(4).neighbors = [1 2 3 5];
node(5).neighbors = [1 ];

NumRequests = Ncache*CatalogSize*5;

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
    TotalRequestCount = zeros(1,Ncache);
    TotalHitRate = zeros(1,Ncache);
    
    
    for index = 1:NumRequests
        c = randsample(Ncache,1);
        
        
        RequestCount(c) = RequestCount(c) +1;
        CID = ContentRequests(index);
        
        
        hit =  emulate(cache(c),ContentRequests(index));
        
        if(hit == 1)
            HitCountSelf(c) = HitCountSelf(c)+1;
        else
            nbrs = node(c).neighbors;
            d = datasample ( nbrs, 1);
            h = emulate(cache(d),CID);
            end
        
    end
    
    stats(i).HitCountSelf = HitCountSelf;
    stats(i).HitCountOther = HitCountOther;
    stats(i).RequestCount = RequestCount;
    stats(i).CacheSize = CacheSize;
    stats(i).CacheSize = CacheSize;
    for index = 1:Ncache
        TotalRequestCount(index) = getRequestCount( cache(index) );
        TotalHitRate(index) = getHitRate(cache(index));
    end
    stats(i).TotalRequestCount = TotalRequestCount;
    stats(i).TotalHitRate = TotalHitRate;
    
    
    %     fprintf('Simulation for cache size %d finished. ',CacheSize);
    i = i+1;
    toc
end



end