function [TON, TOFF, Tc] = derive_ON_OFF_Tc( ts)

    
TON = [];
TOFF = [];
Tc = [];

if length(ts) <= 3
   TON = 0;
   TOFF = 0;
   Tc = 0;
   
   return;
end

prevTON = 0;

if ts(1) < 0
    prevTON = ts(2);
end


for index = 2: length(ts)-1
    
    if ts(index) < 0
        TON = [TON abs(ts(index))-prevTON ];
        TOFF = [TOFF ts(index+1)-abs(ts(index)) ];
        Tc = [Tc abs(ts(index))-ts(index-1)];
        prevTON = abs(ts(index+1));
    end
    
end
if isempty(Tc)
   Tc = abs(ts(length(ts)-1)); 
end

end