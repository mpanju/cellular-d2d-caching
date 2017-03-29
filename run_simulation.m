clear all; clc;

ZipfExponent = 0.8;

 
CatalogSize = 100;
CacheSizeSamples = [ 10] ;

SIM_RUNS = 1;

AdjMatrix = [0];
% AdjMatrix = [0 1 0 0 0 0;
%              1 0 1 1 1 0;
%              0 1 0 0 0 0;
%              0 1 0 0 1 1;
%              0 1 0 1 0 0;
%              0 0 0 1 0 0];

 
%          AdjMatrix = [0 1 1 1 1;
%              1 0 1 1 1;
%              1 1 0 1 1;
%              1 1 1 0 1;
%              1 1 1 1 0];

if sum(sum(AdjMatrix - AdjMatrix')) > 0
    fprintf('Adjacency matrix is not symmetric\n');
    return
else
    Ncache = size(AdjMatrix,1);
end

NumRequests = Ncache*CatalogSize*10;

weights = 1:CatalogSize;
weights = weights.^-ZipfExponent;
% weights = ones(1,CatalogSize);
weights = weights./sum(weights);

%

for index = 1:length(CacheSizeSamples)
    
    LocalHitRate = zeros(length(CacheSizeSamples), Ncache);
    NeighborHitRate = zeros(length(CacheSizeSamples), Ncache);
    TotalHitRate = zeros(length(CacheSizeSamples), Ncache);
    
    LocalHitCountVec = zeros(CatalogSize, Ncache);
    LocalRequestCountVec = zeros(CatalogSize, Ncache);
    
    TotalHitCountVec = zeros( CatalogSize, Ncache);
    TotalRequestCountVec = zeros( CatalogSize, Ncache);
    
    A_LocalHitRate = zeros(1,Ncache);
    A_NeighborHitRate = zeros(1,Ncache);
    A_TotalHitRate = zeros(1,Ncache);
    
    A_TotalHitCount = zeros(1,Ncache);
    A_TotalRequestCount = zeros(1,Ncache);
    
    A_LocalHitCountVec = zeros(CatalogSize,Ncache);
    A_LocalRequestCountVec = zeros(CatalogSize, Ncache);
    
    A_TotalRequestCountVec = zeros(CatalogSize, Ncache);
    A_TotalHitCountVec = zeros(CatalogSize, Ncache);
    
     A_LocalHitCountOtherVec =  zeros(CatalogSize,Ncache);
    
    for i = 1:Ncache
        for j = 1:CatalogSize
            TimerStructArray(i,j).mean_ton = 0;
            TimerStructArray(i,j).mean_toff = 0;
            TimerStructArray(i,j).mean_tc = 0;
            TimerStructArray(i,j).std_tc = 0;
            t.time = [0];
            TimerStructArray(i,j).timeSeries = repmat(t,SIM_RUNS,1);
        end
    end
    
    for sim_run = 1:SIM_RUNS
        
        fprintf('\nStarted simulation for CacheSize %d, SIM_RUN = ( %d/%d): ', CacheSizeSamples(index), sim_run, SIM_RUNS);
        
        stats = network_of_caches( CacheSizeSamples(index) , weights , NumRequests, AdjMatrix);
        
        A_LocalHitRate = A_LocalHitRate + (sum(stats.HitCountSelfVec,1)./sum(stats.RequestCountVec,1));
        A_NeighborHitRate = A_NeighborHitRate + (sum(stats.HitCountOtherVec,1) ./ sum(stats.RequestCountVec,1) );
        A_TotalHitRate = A_TotalHitRate + stats.TotalHitRate;
        
        A_LocalHitCountVec = A_LocalHitCountVec + stats.HitCountSelfVec;
        A_LocalRequestCountVec = A_LocalRequestCountVec + stats.RequestCountVec;
        
        A_TotalHitCountVec = A_TotalHitCountVec + stats.TotalHitCountVec';
        A_TotalRequestCountVec = A_TotalRequestCountVec + stats.TotalRequestCountVec';

        A_LocalHitCountOtherVec = A_LocalHitCountOtherVec + stats.HitCountOtherVec;        
        
        for i = 1:Ncache
            for j = 1:CatalogSize
                TimerStructArray(i,j).timeSeries(sim_run).time = stats.TimerStructArray(i,j).time;
                [TON, TOFF, TC] = derive_ON_OFF_Tc( stats.TimerStructArray(i,j).time);
                TimerStructArray(i,j).mean_ton = TimerStructArray(i,j).mean_ton + mean(TON);
                TimerStructArray(i,j).mean_toff = TimerStructArray(i,j).mean_toff + mean(TOFF);
                TimerStructArray(i,j).mean_tc = TimerStructArray(i,j).mean_tc  + mean(TC);
                TimerStructArray(i,j).std_tc = TimerStructArray(i,j).std_tc + std(TC);
            end
        end
    end
    
    LocalHitRate(index,:) = A_LocalHitRate / SIM_RUNS;
    NeighborHitRate(index,:) = A_NeighborHitRate / SIM_RUNS;
    TotalHitRate(index,:) = A_TotalHitRate / SIM_RUNS;
    
    LocalHitCountVec = A_LocalHitCountVec / SIM_RUNS;
    LocalRequestCountVec = A_LocalRequestCountVec / SIM_RUNS;
    
    TotalHitCountVec = A_TotalHitCountVec / SIM_RUNS;
    TotalRequestCountVec = A_TotalRequestCountVec /SIM_RUNS;
 
    LocalHitCountOtherVec = A_LocalHitCountOtherVec /SIM_RUNS ;    
    
    nodes = stats.nodes;
    
    for i = 1:Ncache
        for j = 1:CatalogSize
            TimerStructArray(i,j).mean_ton = TimerStructArray(i,j).mean_ton / SIM_RUNS;
            TimerStructArray(i,j).mean_toff = TimerStructArray(i,j).mean_toff / SIM_RUNS;
            TimerStructArray(i,j).mean_tc = TimerStructArray(i,j).mean_tc / SIM_RUNS;
            TimerStructArray(i,j).std_tc = TimerStructArray(i,j).std_tc / SIM_RUNS;
        end
    end
    
    clear stats;
    
    stats.LocalHitCountVec = LocalHitCountVec;
    stats.LocalHitCountOtherVec =LocalHitCountOtherVec;
    stats.LocalRequestCountVec = LocalRequestCountVec;
    
    stats.TotalHitCountVec = TotalHitCountVec;
    stats.TotalRequestCountVec = TotalRequestCountVec;
    
    stats.TimerStructArray = TimerStructArray;
    stats.nodes = nodes;

    statsVec(index) = stats;

%     save(sprintf('SIM_%s_c%d_m%d.mat',datestr(now,'yyyy-mm-dd-HH-MM-SS'),CacheSizeSamples(index),CatalogSize))
end

% Calculate Hit Rate

HitLocalRateArray = zeros(Ncache,length(CacheSizeSamples));   
HitTotalRateArray = zeros(Ncache,length(CacheSizeSamples));   
for ondex = 1: Ncache     
    for index = 1:length(CacheSizeSamples) 
        HitLocalRateArray(ondex,index) = sum( statsVec(index).LocalHitCountVec(:,ondex) )./ sum( statsVec(index).LocalRequestCountVec(:,ondex) );   
        HitTotalRateArray(ondex,index) = sum( statsVec(index).LocalHitCountVec(:,ondex) + statsVec(index).LocalHitCountOtherVec(:,ondex) )./ sum( statsVec(index).LocalRequestCountVec(:,ondex) );   
    end   
end     
  
  
%% Analyize ON-OFF duration


PinArray = zeros(Ncache,CatalogSize);
P1Array = zeros(Ncache,CatalogSize);
P2Array = zeros(Ncache,CatalogSize);
meanTcArray = zeros(Ncache,CatalogSize);
stdTcArray = zeros(Ncache,CatalogSize);

for nodeID = 1:Ncache
    for contentID = 1:CatalogSize
        ts = TimerStructArray(nodeID,contentID).timeSeries(1);
        [ Pin, P1, P2, mean_tc, std_tc] = analyze_on_off(nodeID,contentID,stats);
        PinArray(nodeID,contentID) = Pin;
        P1Array(nodeID,contentID) =  P1;
        P2Array(nodeID,contentID) = P2;
        meanTcArray(nodeID,contentID) = mean_tc;
        stdTcArray(nodeID,contentID) = std_tc;
    end
end


%%

for ind = 1:CatalogSize
    [ton, toff, tc] = derive_ON_OFF_Tc( TimerStructArray(2,ind).timeSeries(1).time);
    [f,x] = hist(tc,200);
    dx = diff(x(1:2));
    bar(x,f/sum(f*dx));
    axis([mean(x)-5*std(x) mean(x)+5*std(x) 0 max(f/sum(f*dx))+0.1]);
    legend( sprintf('Content ID = %d',ind));
    pause(1);
end



%% Plot theoretical

CSamples = CacheSizeSamples;

M = CatalogSize;

PhitVec = zeros(Ncache,length(CSamples));

% 
% lambdaArray = 1:M;
% lambdaArray = lambdaArray.^-ZipfExponent;
% lambdaArray = lambdaArray/sum(lambdaArray);

lambdaArray = repmat(weights,Ncache,1);

for index = 1:length(CSamples)
    fprintf('Iteration started for C = %d\n', CSamples(index));

    TcVec = 10*ones(1,Ncache);
    PnmArray = ones(Ncache,M) * 0.5;
    TcVec_ = zeros(1,Ncache);
    PnmArray_ = zeros(Ncache,M);

    iter = 1;
    while norm(PnmArray-PnmArray_) > 1e-6 || norm(TcVec_ - TcVec) > 1e-3
        PnmArray_ = PnmArray;
        TcVec_ = TcVec;
        PnmArray = computePnmArray( AdjMatrix, TcVec_, lambdaArray);
        TcVec = computeTc(AdjMatrix, lambdaArray, PnmArray_, CSamples(index));
        if( rem(iter,5) == 0)
           fprintf('ITERATION %d; Pnm,Pnm_ = %f,%f; TcVec,TcVec_ =  %f,%f \n',iter,norm(PnmArray),norm(PnmArray_),norm(TcVec), norm(TcVec_));
        end
        iter = iter+1;
    end
    for node = 1:Ncache
        PhitVec(node,index) = sum(PnmArray(node,:).*lambdaArray(node,:))/sum(lambdaArray(node,:));
    end
end


%%

ColorProfile = [1 0 0;
                0 1 0;
                0 0 1;
                1 1 0;
                1 0 1];


% PhitArray = computePhitArray( AdjMatrix, TcVec, lambdaArray);
% PhitArray = computePhitArray( AdjMatrix, mean(meanTcArray,2), lambdaArray);

PnmArray = computePnmArray( AdjMatrix, TcVec, lambdaArray);
% PnmArray = computePnmArray( AdjMatrix, mean(meanTcArray,2), lambdaArray);

% errMatrix = abs((PhitArray - TotalHitCountVec'./TotalRequestCountVec')./(TotalHitCountVec'./TotalRequestCountVec'))*100;
errMatrix = abs((PnmArray - LocalHitCountVec'./LocalRequestCountVec')./(LocalHitCountVec'./LocalRequestCountVec'))*100;

%for index = 1:Ncache
    plot( errMatrix(2,:),'Color',ColorProfile(2,:)); hold on;
%     plot( errMatrix(index,:)); hold on;
%end
%% Plot Simulation

ColorProfile = [1 0 0;
                0 1 0;
                0 0 1;
                1 1 0;
                0 1 1];
figure; hold on;
for index = 1:Ncache
    

end

 hold on;
for index = 1:Ncache
    plot(CSamples,PhitVec(index,:));
    
end
%% Miss stream analysis


contentID = 1;
fprintf('\n');
for index = 1:Ncache
   
    n = length(stats.nodes(index).neighbors);
    
    nbrs = '';
    nbrhitcount = '';
    for i = 1:n
        nbrs = strcat(nbrs, sprintf('-%d',stats.nodes(index).neighbors(i)));
        nbrhitcount = strcat(nbrhitcount, sprintf(' %d',(stats.nodes(index).nbrHitCount(contentID,i))));
    end
    
   fprintf(sprintf('Node %d. Neighbors: %s\n',index,nbrs));
   fprintf(sprintf('Node %d. Neighbor Hit Count: %s\n',index,nbrhitcount));
   

   % Excess rate calculation
   PROD_TERM = weights(contentID);
   for i = 1:Ncache
       if(AdjMatrix(index,i) == 1)
           PROD_TERM = PROD_TERM * (1-PnmArray(index,contentID));
       end
   end   
   ts = stats.TimerStructArray(index,contentID).timeSeries(1).time;
   ts = ts(ts>0);
   
%    fprintf(sprintf('Node %d. ContentID=%d: \n',index,contentID));
   fprintf(sprintf('Node %d. ContentID=%d: Arrival Rate = %.4f, Request Rate = %.4f, Excess Rate = %f \n',index,contentID,1/mean(diff(ts)),weights(contentID) ,  PROD_TERM ));
   
   fprintf('\n\n');
end


%%


nodeId = 5;

Lprime = computeLambdaArray( PnmArray, lambdaArray, AdjMatrix);

plot(1:100, TotalRequestCountVec(:,nodeId)/sum(TotalRequestCountVec(:,nodeId)));
hold on; 
plot(1:100,LocalRequestCountVec(:,nodeId)/sum(LocalRequestCountVec(:,nodeId)),'r');



plot(1:CatalogSize, Lprime(nodeId,:)/sum(Lprime(nodeId,:)),'k' ); hold on;
plot(1:CatalogSize, lambdaArray(nodeId,:)/sum(lambdaArray(nodeId,:)),'--k' );


%%

plot(1:CatalogSize, abs ( ( TotalRequestCountVec(:,nodeId)-LocalRequestCountVec(:,nodeId)) ./LocalRequestCountVec(:,nodeId) ) );



















