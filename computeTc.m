function TcVec = computeTc(AdjMatrix, lambdaArray, PnmArray, C)
%     PmVec =


% PnmArray = diag([1 1.5 2 2.5 3 ]) * ones(5,10);

% Ph should not change when units of M and C are changed
%

N = size(PnmArray,1);
M = size(PnmArray,2);
% C = 3;
% lambdaArray = ones(N,M);


%
TcVec = ones(1,N);

iteration_count = 0;

LambdaArray = computeLambdaArray(PnmArray, lambdaArray,AdjMatrix);

for node = 1:N
    Tc = 0;
    expSUM = 0;
%     fprintf('expSUM = %f, (M,N,C) = (%d,%d,%d), Tc = %d\n',expSUM, M, N, C, Tc);
    while (expSUM-C) < 0
        %         [M C M-C Tc iteration_count expSUM]
        
        Numerator = 1 - exp(- LambdaArray(node,:).*Tc);
        Denominator = 1+ (LambdaArray(node,:)./lambdaArray(node,:) - 1).*exp(- LambdaArray(node,:).*Tc);
        
        
        expSUM = sum( Numerator./Denominator );
        Tc = Tc +0.001;
        %         if( iteration_count <= 2)
%         fprintf('expSUM = %f, (M,N,C) = (%d,%d,%d), Tc = %d\n',expSUM, M, N, C, Tc);
        
        %         end
        iteration_count = iteration_count +1;
        %         pause(0.1);
    end
    TcVec(node) = Tc;
end

end