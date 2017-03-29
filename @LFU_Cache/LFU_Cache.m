classdef LFU_Cache < handle
    properties

     data
     counter
     size
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
        
        function cache = LFU_Cache(cacheSize,catalogSize)
            if (nargin > 0)
                cache.size = cacheSize;
                cache.catalogSize = catalogSize;
                cache.state = zeros(1,catalogSize);
                
                cache.data = zeros(1,cacheSize);
                cache.counter = zeros(1,cacheSize);

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
              if isempty(probe) == 1 % Data not found in the cache.
                  probe = find( this.counter == min(this.counter), 1 );
                  if this.data(probe) ~=0
                      this.state( this.data(probe)) = 0;
                  end
                  this.state( CID) = 1;

                  this.data(probe) = CID;
                  this.counter(probe) =  1;
                  hit = 0;
              else % Data found in the cache.
                  this.data(probe) = CID;
                  this.counter(probe) = this.counter(probe) + 1;
                  
                  this.statsHitCountVec(CID) = this.statsHitCountVec(CID) +1;
                  hit = 1;
              end
              %               printCache(this);
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