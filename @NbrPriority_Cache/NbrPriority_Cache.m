classdef NbrPriority_Cache < handle
    properties
        head
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
        
        function cache = NbrPriority_Cache(cacheSize,catalogSize)
            if (nargin > 0)
                cache.size = cacheSize;
                cache.head = dlnode(1);
                cache.catalogSize = catalogSize;
                cache.state = zeros(1,catalogSize);

                cache.statsHitCountVec = zeros(1,catalogSize);
                cache.statsRequestCountVec = zeros(1,catalogSize);
            end
            
            initContents = randperm(catalogSize,cacheSize);
            cache.head = dlnode(initContents(1));
            cache.state(initContents(1)) = 1;
            
            for k = 2:cacheSize
                %                 nextNode = dlnode(randsample(catalogSize,1,'true'));
                nextNode = dlnode(initContents(k));
                %    insertAfter(nextNode,head)
                nextNode.insertAfter(cache.head); % insert nextNode after cache.head
                cache.state(initContents(k)) = 1;
            end
            
            cache.last = cache.head;
            while ~isempty(cache.last.Next)
                cache.last = cache.last.Next;
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

        function hit = emulate(this,CID,time,rank)
            
            this.statsRequestCountVec(CID) = this.statsRequestCountVec(CID) +1;
%             printCache(this);
%             fprintf('head = %d last = %d \n',this.head.Data, this.last.Data);
            this.Timer(CID).time = [this.Timer(CID).time time];

            if rank < 1 || rank > this.size
                fprintf('InsertPosition is out of bounds\n');
            end
            
            InsertPosition = rank;
            
            if this.state(CID) == 1
                this.statsHitCountVec(CID) = this.statsHitCountVec(CID) + 1;
                
                % Find the content in the cache
                temp = this.head;
                ContentPosition = 1;
                while ~isempty(temp.Next) && temp.Data ~= CID
                    temp = temp.Next;
                    ContentPosition = ContentPosition +1;
                end
                
                if(temp.Data == CID)
                    if ContentPosition ~= InsertPosition % Do not do anything if insertposition is same as content's position

                        if ContentPosition ==1 
                            this.head = this.head.Next;
                        elseif ContentPosition == this.size
                            this.last = this.last.Prev;
                        end
                        
                        removeNode(temp);

                        if InsertPosition == 1
                            insertBefore(temp, this.head);
                            this.head = this.head.Prev;
                        elseif InsertPosition == this.size
                            insertAfter(temp, this.last);
                            this.last = temp;
                        else
                            temp2 = this.head;
                            for  ix = 2: InsertPosition -1
                                temp2 = temp2.Next;
                            end
                            
                            insertAfter(temp,temp2);
                        end
                        
                    end
                    hit = 1;
                else
                    fprintf('Something Wrong! Content not found. But, content state = 1\n');
                    return;
                end
            else
                this.state(CID) = 1;
                
                NewNode = dlnode(CID);
                
                if InsertPosition == 1
                    insertBefore(NewNode, this.head);
                    this.head = this.head.Prev;
                elseif InsertPosition == this.size
                    insertBefore(NewNode, this.last);
                else
                   
                    temp = this.head;
                    for  ix = 2: InsertPosition -1
                       temp = temp.Next; 
                    end
                    
                    insertAfter(NewNode,temp);
                end
                
                this.state(this.last.Data) = 0;
                this.Timer(this.last.Data).time = [this.Timer(this.last.Data).time -1*time];
                this.last = this.last.Prev;
                removeNode(this.last.Next);
                
                hit = 0;
            end
            
%             printCache(this);
%             this.state
%             fprintf('head = %d last = %d \n',this.head.Data, this.last.Data);
        end
        
        
        
        
%         function hit = emulate(this,CID,time)
%             
%             this.statsRequestCountVec(CID) = this.statsRequestCountVec(CID) +1;
%             temp = this.head;
% %             printCache(this);
% %             fprintf('head = %d last = %d \n',this.head.Data, this.last.Data);
%             this.Timer(CID).time = [this.Timer(CID).time time];
% 
%             if this.head.Data ~= CID
%                 while ~isempty(temp.Next) && temp.Data ~= CID
%                     temp = temp.Next;
%                 end
%                 
%                 if temp.Data == CID
%                     
%                     if isempty(temp.Next)
% %                         fprintf('__________head = %d last = %d \n',this.head.Data, this.last.Data);
%                         this.last = this.last.Prev;
% %                         fprintf('__________head = %d last = %d \n',this.head.Data, this.last.Data);
%                     end
%                     
%                     temp.removeNode();
%                     temp.insertBefore(this.head);
%                     this.head = temp;
%                     this.statsHitCountVec(CID) = this.statsHitCountVec(CID) +1;
%                     hit = 1;
%                 elseif isempty(temp.Next)
%                     this.state(CID) = 1;
%                     NewNode = dlnode(CID);
%                     insertBefore(NewNode,this.head);
%                     this.head = NewNode;
%                     this.state(this.last.Data) = 0;
%                     this.Timer(this.last.Data).time = [this.Timer(this.last.Data).time -1*time];
%                     this.last = this.last.Prev;
%                     removeNode(this.last.Next);
%                     hit = 0;
%                 end
%             else 
%                 this.statsHitCountVec(CID) = this.statsHitCountVec(CID) +1;
%                 hit = 1;
%             end
% %             printCache(this);
% %             fprintf('head = %d last = %d \n',this.head.Data, this.last.Data);
%         end
        
        function printCache(this)
            temp = this.head;
            
            fprintf('\n');
            while ~isempty(temp.Next) 
                fprintf('%d ',temp.Data);
                temp = temp.Next;
            end                
            fprintf('%d\n',temp.Data);
            
        end
        
        function result = isContentPresent(this,CID)
            
            result = this.state(CID);
%             temp = this.head;
%             while ~isempty(temp.Next) && temp.Data ~= CID
%                 temp = temp.Next;
%             end
%             
%             if temp.Data == CID
%                 result = 1;
%             else
%                 result = 0;
%             end
        end
        
        function hit = softEmulate(this,CID)
            
            this.statsRequestCountVec(CID) = this.statsRequestCountVec(CID) +1;
            temp = this.head;
            while ~isempty(temp.Next) && temp.Data ~= CID
                temp = temp.Next;
            end
            
            if temp.Data == CID
                this.statsHitCountVec(CID) = this.statsHitCountVec(CID) +1;
                hit = 1;
            else
                hit = 0;
            end            
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
