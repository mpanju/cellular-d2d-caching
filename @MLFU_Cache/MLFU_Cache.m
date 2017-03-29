classdef MLFU_Cache < handle
    properties

     data
     counter
     timeArrival
     timeRequest
     size
     ttlThreshold
    end

   properties %(SetAccess = private)
        last
        popularityProfile
        catalogSize

        state
        
        statsRequestCountVec
        statsHitCountVec
        
        Timer
   end    
    
    methods
        
        function cache = MLFU_Cache(cacheSize,catalogSize)
            if (nargin > 0)
                cache.size = cacheSize;
                cache.catalogSize = catalogSize;
                cache.state = zeros(1,catalogSize);

                
                cache.data = -1*(1:cacheSize);
                cache.counter = zeros(1,cacheSize);
                cache.timeArrival = -2*ones(1,cacheSize);
                cache.timeRequest = -1*(1:cacheSize);
                
                cache.statsHitCountVec = zeros(1,catalogSize);
                cache.statsRequestCountVec = zeros(1,catalogSize);
                
            end
            
            t.time = [];
            cache.Timer = repmat(t,1,catalogSize);
%             printCache(cache);
        end

        function setPopularityProfile(this,weights)
            this.popularityProfile = weights;
        end
        
        function setCatalogSize(this,catalogSize)
            this.catalogSize = catalogSize;
            this.statsHitCountVec = zeros(1,catalogSize);
            this.statsRequestCountVec = zeros(1,catalogSize);
        end
        
        function hit = emulate(this,CID,time)
            
            this.statsRequestCountVec(CID) = this.statsRequestCountVec(CID) +1;
            this.Timer(CID).time = [this.Timer(CID).time time];
            
            probe = find(this.data == CID, 1 );
            if length(probe) <= 1
                if isempty(probe) == 1 % Data not found in the cache.
%                     expired = find( this.timeArrival < (time -
%                     this.ttlThreshold)); Commented out to make the cache
%                     LFU
                    expired = 1:this.size;    
                    if isempty(expired) ~= 1 % There are some expired contents
                        arrivalRate = this.counter(expired)./(time - this.timeArrival(expired));
%                         this.counter(expired)./this.timeArrival(expired)
%                         min(arrivalRate(expired))
                        replace_ = find( arrivalRate == min(arrivalRate),1);
                        replace = expired(replace_);
                        this.data(replace) = CID;
                        this.counter(replace) =  1;
                        this.timeArrival(replace) = time;
                        this.timeRequest(replace) = time;
                    else % None of the contents have expired; Follow LRU. After making the cache LFU
                        This part should never run
                        
                        replace = find( this.timeRequest == min(this.timeRequest));
                        if(length(replace) ~= 1)
                            fprintf('\n\nLRU part not right \n');
                        end
                        this.data(replace) = CID;
                        this.counter(replace) =  1;
                        this.timeArrival(replace) = time;
                        this.timeRequest(replace) = time;
                    end
                    hit = 0;
                else % Data found in the cache.
                    this.data(probe) = CID;
                    this.counter(probe) = this.counter(probe) + 1;
                    this.timeRequest(this.data==CID) = time;
                    this.statsHitCountVec(CID) = this.statsHitCountVec(CID) +1;
                    
                    hit = 1;
                end
            else % Duplicate contents found
                fprintf('Duplicate contents found.\n');
                printCache(this)
            end
            %                             printCache(this);
            %             fprintf('head = %d last = %d \n',this.head.Data, this.last.Data);
        end
        
        function printCache(this)
            fprintf('\n');
            for index = 1:length(this.data)
                fprintf('%d ',this.data(index));
            end
            fprintf('\n');
            for index = 1:length(this.counter)
                fprintf('%d ',this.counter(index));
            end
            
            fprintf('\n');
            for index = 1:length(this.timeArrival)
                fprintf('%.2f ',this.timeArrival(index));
            end
                    
            fprintf('\n');
            for index = 1:length(this.timeRequest)
                fprintf('%.2f ',this.timeRequest(index));
            end
            fprintf('\n');
        end
        
        function result = isContentPresent(this,CID)
            
            result = this.state(CID);
        end
        
        
        function hr =  getHitRate(this)
            hr = sum(this.statsHitCountVec)./sum(this.statsRequestCountVec);
        end
        
        
        function trc =  getRequestCount(this)
            trc = sum(this.statsRequestCountVec);
        end
        
        
        function trcVec =  getRequestCountVec(this)
            trcVec = this.statsRequestCountVec;
        end
                

        function thcVec =  getHitCountVec(this)
            thcVec = this.statsHitCountVec;
        end 


        function tsa =  getTimerStructArray(this)
            tsa = this.Timer;
        end 
                        
    end
end