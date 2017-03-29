

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%












%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
cache = LRU_Cache(10);

temp = cache.head;
while ~isempty(temp.Next)
    temp = temp.Next;
end

last = temp;
%%

clear all;clc
cache = LRU_Cache(10);
% printCache(cache);
emulate(cache,2);
% printCache(cache);

%%
clear all; clc;

ZipfExponent = 0.8;
CatalogSize = 1000;

NumRequests = 100000;

weights = 1:CatalogSize;
weights = weights.^-ZipfExponent;

CacheSizeSamples = [50:50:CatalogSize];


PhitVecArray = zeros(11,1000);
cacheSizeIndex = 0;

%%
for run = 1:1
    PhitVec = [];
    ContentRequests = randsample(CatalogSize,NumRequests,'true',weights);
    for CacheSize = CacheSizeSamples
        cacheSizeIndex = cacheSizeIndex +1;
        CacheLRU = LRU_Cache(CacheSize,CatalogSize); % Cache is initially populated with most popular contents
        
        HitCount = 0;
        for index = 1:NumRequests
            HitCount = HitCount + emulate(CacheLRU,ContentRequests(index));
            
        end
        PhitVec = collectStats(CacheLRU); % 1xCatalogSize
        PhitVecArray(cacheSizeIndex,:) = PhitVec;
        PhitVec = [PhitVec HitCount/NumRequests]
        [run CacheSize sum(PhitVec.*weights)/sum(weights)]
    end
end

% save('Simulation.mat')

%% Plot theoretical for two interacting caches (Assuming Phit and Pin are same)
%clear all;clc;

ZipfExponent = 0.8;
CatalogSize = 1000;

weights = 1:CatalogSize;
weights = weights.^-ZipfExponent;
weights = weights/sum(weights);

tc = 10.^(1:0.1:4);

C = zeros(1,size(tc,2));
PHm = zeros(size(tc,2),size(weights,2));
missRatem = zeros(size(tc,2),size(weights,2));
PhitVec = zeros(1,size(tc,2));
missRate = zeros(1,size(tc,2));
for t = 1:size(tc,2)
    for index = 1:size(weights,2)
        pin = 0; pin_ = 1;
        while abs(pin_ - pin) > 1e-6
            pin_ = pin;
            pin = (1-exp(- tc(t).*weights(index).*(1+pin.*(1-pin))))./(1+pin.*(1-pin).*exp(-tc(t).*weights(index).*(1+pin.*(1-pin))));
        end
        PHm(t,index) = pin;
        
        missRatem(t,index) = weights(index)*(1-pin)*exp(-weights(index)*(1-pin)*tc(t)) / ( weights(index)*(1+pin*(1-pin)*exp(-weights(index)*(1-pin)*tc(t))) );
    end
    C(t)  = sum(PHm(t,:));
    PhitVec(t) = sum( PHm(t,:).*weights );
    missRate(t) = sum( missRatem(t,:).*weights );
end
plot(C,PhitVec);
plot(C,missRate,'k');
%% Plot theoretical for two interacting caches (Phit and Pin are not assumed to be the same)
%clear all;clc;

ZipfExponent = 0.8;
CatalogSize = 1000;

weights = 1:CatalogSize;
weights = weights.^-ZipfExponent;
weights = weights/sum(weights);

tc = 10.^(1:0.1:4);

C = zeros(1,size(tc,2));
PHm = zeros(size(tc,2),size(weights,2));
PINm = zeros(size(tc,2),size(weights,2));
PhitVec = zeros(1,size(tc,2));
Pin = zeros(1,size(tc,2));
for t = 1:size(tc,2)
    for index = 1:size(weights,2)
        pVec = [0 0];
        pVec_ = [1 1];
        pin_ = 1;
        pin = 0; ph = 0;
%         while abs(pin_ - pin) > 1e-3
        while norm(pVec_ - pVec)> 1e-3
            pVec_ = pVec;
            pin_ = pin;
            L = weights(index).*(1+(1-ph));
            pin = (1-exp(- tc(t).*L))./(1+(1-ph).*exp(-tc(t).*L));
            ph = (1 -(1+L*tc(t))*exp(- tc(t).*L))/(1+(1-ph)*exp(- tc(t).*L));
%             ph = (1 -exp(- tc(t).*L));
            pVec = [pin ph];
        end
        
        PHm(t,index) = ph;
        PINm(t,index) = pin;
    end
    C(t)  = sum(PINm(t,:));
    PhitVec(t) = sum( PHm(t,:).*weights );
    Pin(t) = sum( PINm(t,:).*weights );
end

%%
hold on;

plot(C,Pin,'r');
plot(C,Pin+PhitVec.*(1-Pin),'k');
plot(C,PhitVec.*(1-Pin));
grid on;

%%
clear all; clc


ZipfExponent = 0.8;
CatalogSize = 1000;

weights = 1:CatalogSize;
weights = weights.^-ZipfExponent;
weights = weights/sum(weights);

tc = 10.^(1:0.1:4);

C = zeros(1,size(tc,2));
PHm = zeros(size(tc,2),size(weights,2));
PINm = zeros(size(tc,2),size(weights,2));
PhitVec = zeros(1,size(tc,2));
Pin = zeros(1,size(tc,2));

for t = 1:size(tc,2)
    for index = 1:size(weights,2)
        
        eTerm = exp(-weights(index)*tc(t));
        ph = 1 - ( 0.5*eTerm + eTerm^2 *( 0.5*eTerm +  weights(index)*tc(t) )  );
        
        L = weights(index)*( 1 + ph*(1-ph));
        eTerm = exp(-L*tc(t));
        pin = (1 - eTerm) / (1 + (1-ph)*eTerm );
        
        
        PHm(t,index) = ph;
        PINm(t,index) = pin;
    end
    C(t)  = sum(PINm(t,:));
    PhitVec(t) = sum( PHm(t,:).*weights );
    Pin(t) = sum( PINm(t,:).*weights );
end

plot(C,Pin); hold on;
plot(C,PhitVec,'r');

%%

tc = 1:1000;
y1 = 1:length(tc);
y2 = 1:length(tc);

for index = 1:length(tc)

        eTerm = exp(-weights(index)*tc(t));
        ph = 1 - ( 0.5*eTerm + eTerm^2 *( 0.5*eTerm +  weights(index)*tc(t) )  );
        
        L = weights(index)*( 1 + ph*(1-ph));
        eTerm = exp(-L*tc(t));
        pin = (1 - eTerm) / (1 + ph*(1-ph)*eTerm );
        
        y1(index) = ph;
        y2(index) = pin;
end      

plot(y1); hold on;
plot(y2,'r');

%%


tc = 10.^(1:0.2:4);

C = zeros(Ncache,size(tc,2));
Pnm = zeros(Ncache,size(weights,2));
PhitVec = zeros(Ncache,size(tc,2));
PhitT = zeros(Ncache,size(tc,2));

weights = ones(1,CatalogSize);
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

    C(:,t)  = sum(Pnm,2);
    PhitVec(:,t) = sum( Pnm.* repmat(weights,Ncache,1),2 );
    PhitT(:,t) = sum( (1- (1-Pnm).^(Nd+1) ).*repmat(weights,Ncache,1) ,2 );
end

%%


p = [0:0.05:1];

figure; hold on;

LAMBDA = 1-p;                

for index = 0.1:0.5:10
EXP_TERM = exp(- index.*(1+LAMBDA));
                
f = (1-EXP_TERM)./(1+(LAMBDA).*EXP_TERM);

plot(p,f);

end

%%
 clear all; clc;
% AdjMatrix = [0 1 1 1 1;
%     1 0 1 1 1;
%     1 1 0 1 1;
%     1 1 1 0 1;
%     1 1 1 1 0];
AdjMatrix = 0;
N = 1;
M = 1000;

CacheSizeSamples = 100:100:900;

CacheSizeSamples = 100;

PhitVec = zeros(N,length(CacheSizeSamples));

ZipfExponent = 0.7;
lambdaArray = 1:M;
lambdaArray = lambdaArray.^-ZipfExponent;
lambdaArray = repmat(lambdaArray,N,1);



for index = 1:length(CacheSizeSamples)
    fprintf('Iteration started for C = %d\n', CacheSizeSamples(index));

    TcVec = 10*ones(1,N);
    PnmArray = ones(N,M) * 0.5;
    TcVec_ = zeros(1,N);
    PnmArray_ = zeros(N,M);

    
    while norm(PnmArray-PnmArray_) > 0.1 || norm(TcVec_ - TcVec) > 0.001
        PnmArray_ = PnmArray;
        TcVec_ = TcVec;
        PnmArray = computePnmArray( AdjMatrix, TcVec, lambdaArray);
        TcVec = computeTc(AdjMatrix, lambdaArray, PnmArray, CacheSizeSamples(index));
    end
    for node = 1:N
        PhitVec(node,index) = sum(PnmArray(node,:).*lambdaArray(node,:))/sum(lambdaArray(node,:));
    end
end

%%
figure; hold on;
for index = 1:N
    plot(CacheSizeSamples,PhitVec(index,:));
    
end



%%























