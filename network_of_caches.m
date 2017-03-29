function  stats = network_of_caches(CacheSize, weights, NumRequests, AdjMatrix)


CatalogSize = length(weights);

if length(CacheSize) ~= 1
    fprintf('CacheSize must be a positive integer\n');
end


if CacheSize > CatalogSize
    fprintf('CacheSize must be less than or equal to catalog size\n');
end

if size(AdjMatrix,1) == size(AdjMatrix,2)
    Ncache = size(AdjMatrix,1);
else
    fprintf('Adjacency Matrix is not a square matrix\n');
    return;
end

for index = 1:Ncache
    nodes(index).neighbors = find(AdjMatrix(index,:) == 1);
end



tic

for k = 1:Ncache
    cache(k) = LRU_Cache(CacheSize,CatalogSize);
end

for index = 1:Ncache
    nodes(index).nbrHitCount = zeros(CatalogSize,length( nodes(index).neighbors));
end

ContentRequests = randsample(CatalogSize,NumRequests,'true',weights);
TimeOfRequests = cumsum(exprnd( ones(1,NumRequests)/Ncache));
NodeOfRequests = randsample(Ncache,NumRequests,'true');


HitCountSelfVec = zeros(CatalogSize,Ncache);
HitCountOtherVec = zeros(CatalogSize,Ncache);
RequestCountVec = zeros(CatalogSize,Ncache);
TotalRequestCount = zeros(1,Ncache);
TotalHitRate = zeros(1,Ncache);
TotalRequestCountVec = zeros(Ncache,CatalogSize);
TotalHitCountVec = zeros(Ncache,CatalogSize);


for index = 1:NumRequests
    c = NodeOfRequests(index);
    CID = ContentRequests(index);
    ToR = TimeOfRequests(index);
    
    RequestCountVec(CID,c) = RequestCountVec(CID,c) +1;
    
    
    if( rem(index, round(NumRequests/5) ) == 0)
       fprintf('...%d ',20*round(index/round(NumRequests/5)));
    end
%     fprintf('Request No %d: CID = %d, ', index, CID);
    hit =  emulate(cache(c),ContentRequests(index),ToR);
    
    if(hit == 1)
        HitCountSelfVec(CID,c) = HitCountSelfVec(CID,c)+1;
%         fprintf('Found locally\n');
    else
        nbrs = nodes(c).neighbors;
        probeResult = zeros(1,length(nbrs));
        
        for j = 1:length(nbrs)
            probeResult(j) = isContentPresent(cache(nbrs(j)),CID);
        end
        
        if (sum(probeResult) ~= 0)
            d = datasample ( nbrs(probeResult == 1), 1);
            nodes(c).nbrHitCount ( CID,nodes(c).neighbors == d) = nodes(c).nbrHitCount ( CID, nodes(c).neighbors == d) + 1;
            HitCountOtherVec(CID,c) = HitCountOtherVec(CID,c) +1;
%             fprintf('Found in neighbor\n');
            h = emulate(cache(d),CID,ToR);
            if( h ==0)
                fprintf('Somthing Wrong!\n');
                return;
            end
        else
%            fprintf('Not Found\n'); 
        end
    end
    
end
HitCountOtherVec;
stats.HitCountSelfVec = HitCountSelfVec;
stats.HitCountOtherVec = HitCountOtherVec;
stats.RequestCountVec = RequestCountVec;
stats.CacheSize = CacheSize;

stats.TimerStructArray = [];
for index = 1:Ncache
    TotalRequestCount(index) = getRequestCount( cache(index) );
    TotalHitRate(index) = getHitRate(cache(index));
    TotalHitCountVec(index,:) = getHitCountVec( cache(index));
    TotalRequestCountVec(index,:) = getRequestCountVec( cache(index));
    stats.TimerStructArray = [stats.TimerStructArray ;getTimerStructArray(cache(index))];
end
stats.TotalRequestCount = TotalRequestCount;
stats.TotalHitRate = TotalHitRate;
stats.TotalRequestCountVec = TotalRequestCountVec;
stats.TotalHitCountVec = TotalHitCountVec;

stats.nodes = nodes;

stats.TimeOfRequests = TimeOfRequests;
stats.ContentRequests = ContentRequests;
stats.NodeOfRequests = NodeOfRequests;

toc

end



