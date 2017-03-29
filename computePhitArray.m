function PhitArray = computePhitArray( AdjMatrix, TcVec, lambdaArray)

% AdjMatrix = [0 1 1 1 1;
%     1 0 1 1 1;
%     1 1 0 1 1;
%     1 1 1 0 1;
%     1 1 1 1 0];

N = size(lambdaArray,1);
M = size(lambdaArray,2);

PhitArray = 0.5* ones(N,M);
% TcVec = ones(1,N);


% ZipfExponent = 0.7;
% lambdaArray = 1:M;
% lambdaArray = lambdaArray.^-ZipfExponent;
% lambdaArray = repmat(lambdaArray,N,1);


for index = 1:M
    phit_ = zeros(1,N);
    phit = ones(1,N);
    
    while norm(phit_ - phit) > 1e-6
        phit_ = phit;
        
        for node = 1:N
            EXT_RATE = 0;
            for j = 1:N
                if(AdjMatrix(node,j) == 1)
                    EXT_RATE = EXT_RATE + (1 - phit_(j))/sum(AdjMatrix(j,:));
                end
            end
            
            LAMBDA = lambdaArray(node,index) *(1+EXT_RATE);
            EXP_TERM = exp( - LAMBDA*TcVec(node) );
            
            phit(node) = (1 - EXP_TERM);
            
        end
    end
    
    PhitArray(:,index) = phit;
end

end