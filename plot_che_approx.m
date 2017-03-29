clear all; 

ZipfExponent = 0.8;
CatalogSize = 1000;

weights = 1:CatalogSize;
weights = weights.^-ZipfExponent;

weights = weights/sum(weights);

Nc = 10:10:10000;

C = zeros(1,size(Nc,2));
Ph = zeros(1,size(Nc,2));

for index = 1:size(Nc,2)
    temp = 1-(1-weights).^Nc(index);
    C(index)= sum(temp);
    Ph(index) = sum(weights.*temp);
end
plot(C,Ph,'r');
%%
hold on;
phitSimulation =  [0.2613    0.3773    0.4578    0.5216    0.5754    0.6224    0.6642    0.7024 0.7372    0.7692    0.7991    0.8274    0.8537    0.8786    0.9019    0.9238 0.9449    0.9648    0.9833    1.0000];

cacheSizeSamples = 50:50:1000;
plot(cacheSizeSamples,phitSimulation,'xr');