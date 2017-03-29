function [Pin, P1 ,P2 ,mean_tc, std_tc] = analyze_on_off(nodeID,contentID,stats)

TimerStructArray = stats.TimerStructArray;
LocalHitCountVec = stats.LocalHitCountVec;
LocalRequestCountVec = stats.LocalRequestCountVec;
TotalHitCountVec = stats.TotalHitCountVec;
TotalRequestCountVec = stats.TotalRequestCountVec;

sim_runs = length(TimerStructArray(nodeID,contentID).timeSeries);

mean_ton = 0;
mean_toff = 0;
mean_tc = 0;
std_tc = 0;
sum_ton = 0;
sum_toff = 0;

for index = 1 : sim_runs
    ts = TimerStructArray(nodeID,contentID).timeSeries(index).time;
    [TON, TOFF, Tc] = derive_ON_OFF_Tc(ts);
    
    sum_ton = sum_ton + sum(TON);
    sum_toff = sum_toff + sum(TOFF);
    mean_ton = mean_ton + mean(TON);
    mean_toff = mean_toff + mean(TOFF);
    mean_tc = mean_tc + mean(Tc);
    std_tc = std_tc + std(Tc);
    
end


mean_ton = mean_ton / sim_runs;
mean_toff = mean_toff / sim_runs;
mean_tc = mean_tc / sim_runs;
std_tc = std_tc / sim_runs;
sum_ton = sum_ton /sim_runs;
sum_toff = sum_toff / sim_runs;

Pin  = sum_ton/(sum_ton + sum_toff);

P1 = LocalHitCountVec(contentID,nodeID)/LocalRequestCountVec(contentID, nodeID);
P2 = TotalHitCountVec(contentID,nodeID)/TotalRequestCountVec(contentID,nodeID);


% [Pin, P1 ,P2 ,mean_tc, std_tc]

end