function PnmArray = computePnmArray( AdjMatrix, TcVec, lambdaArray)

% AdjMatrix = [0 1 1 1 1;
%     1 0 1 1 1;
%     1 1 0 1 1;
%     1 1 1 0 1;
%     1 1 1 1 0];

N = size(lambdaArray,1);
M = size(lambdaArray,2);

PnmArray = 0.5* ones(N,M);
% TcVec = ones(1,N);


% ZipfExponent = 0.7;
% lambdaArray = 1:M;
% lambdaArray = lambdaArray.^-ZipfExponent;
% lambdaArray = repmat(lambdaArray,N,1);


for index = 1:M
    pin_ = zeros(1,N);
    pin = ones(1,N);
    
    while norm(pin_ - pin) > 1e-6
        pin_ = pin;
        
        for node = 1:N
            EXT_RATE = 0;
            for j = 1:N
                if(AdjMatrix(node,j) == 1)
                    EXT_RATE = EXT_RATE + (1 - pin_(j))/sum(AdjMatrix(j,:));
                end
            end
            
            LAMBDA = lambdaArray(node,index) *(1+EXT_RATE);
            EXP_TERM = exp( - LAMBDA*TcVec(node) );
            
            pin(node) = (1 - EXP_TERM)/ (1 + (LAMBDA/lambdaArray(node,index) - 1)*EXP_TERM);
            
        end
    end
    
    PnmArray(:,index) = pin;
end

end