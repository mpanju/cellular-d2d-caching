classdef LFRU_Cache < handle
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
        
        function cache = LFRU_Cache(cacheSize,catalogSize,ttlT)
            if (nargin > 0)
                cache.size = cacheSize;
                cache.catalogSize = catalogSize;
                cache.state = zeros(1,catalogSize);
                cache.ttlThreshold = ttlT;
                
                cache.data = -1*(1:cacheSize);
                cache.counter = zeros(1,cacheSize);
                cache.timeArrival = -2*ttlT*ones(1,cacheSize);
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
            
            this.statsRequestCountVec(CID) = this.statsRequestCountVec(CID) +1; % Update Request counter
            this.Timer(CID).time = [this.Timer(CID).time time]; % Record time of request/arrival 
%             this.Timer(CID).time
            
            probe = find(this.data == CID, 1 );
            if length(probe) <= 1
                if isempty(probe) == 1 % Data not found in the cache.
                    expired = find( this.timeArrival < (time - this.ttlThreshold)); % Obtain list of expired contents

                    if isempty(expired) ~= 1 % There are some expired contents
                        replace_ = find( this.counter(expired) == min(this.counter(expired))); % Find the expired content with minimum count
%                         replace_ = find( this.counter(expired) == max(this.data(expired)),1); % Find the expired content with least popularity
                        replace_ = datasample(replace_,1);
                        replace = expired(replace_);
                        if this.data(replace) ~= 0
                            this.state(this.data(replace)) = 0;
                        end
                        this.Timer(this.data(replace)) = [this.Timer(this.data(replace)) -1*time]; % Record time of departure
                        this.data(replace) = CID; % insert requested content in the place of removed content and update conter and timer.
                        this.state(CID) = 1;
                        this.counter(replace) =  1;
                        this.timeArrival(replace) = time;
                        this.timeRequest(replace) = time;
                    else % None of the contents have expired; Follow LRU
                        replace = find( this.timeRequest == min(this.timeRequest)); % Find the LRU content
                        if(length(replace) ~= 1)
                            fprintf('\n\nLRU part not right \n');
                        end
                        this.Timer(this.data(replace)) = [this.Timer(this.data(replace)) -1*time]; % Record time of departure
                        if this.data(replace) ~= 0
                            this.state(this.data(replace)) = 0;
                        end
                        this.data(replace) = CID;% insert requested content in the place of removed content and update conter and timer.
                        this.counter(replace) =  1;
                        this.state(CID) = 1;                        
                        this.timeArrival(replace) = time;
                        this.timeRequest(replace) = time;
                    end
                    hit = 0;
                else % Data found in the cache.
                    this.data(probe) = CID; % This line is not really necessary
                    this.counter(probe) = this.counter(probe) + 1; % Update counter and timers
                    this.timeRequest(this.data==CID) = time;
                    this.statsHitCountVec(CID) = this.statsHitCountVec(CID) +1; % Record a hit event
                    
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