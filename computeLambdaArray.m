function LambadaArray = computeLambdaArray( PnmArray,lambadaArray, AdjMatrix)
%
% PnmArray = diag([1 1.5 2 2.5 3 ]) * ones(5,10);
%
% % Ph should not change when units of M and C are changed
% %
% node = 1;
% AdjMatrix = [0 1 1 1 1;
%     1 0 1 1 1;
%     1 1 0 1 1;
%     1 1 1 0 1;
%     1 1 1 1 0];
%
%
%

%
% lambadaArray = ones(N,M);

N = size(PnmArray,1);
M = size(PnmArray,2);

LambadaArray = ones(N,M);

for node = 1:N
    for index =  1:M
        
        EXT_RATE = 0;
        for j = 1:N
            if(AdjMatrix(node,j) == 1)
                EXT_RATE = EXT_RATE + (1 - PnmArray(node,j))/sum(AdjMatrix(j,:));
            end
        end
        
        LambadaArray(node,index) =  lambadaArray(node,index) .* (1 + (EXT_RATE));
%         fprintf('Node %d, Content %d, EXT_RATE = %f\n',node,index,EXT_RATE);
    end
end

end