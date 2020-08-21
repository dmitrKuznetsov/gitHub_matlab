function [y,W] = my_Frost(s_in, J, K, mu)
%   frost1972.pdf
%   J - filter order
%   K - Number of microphones
%   s_in - input signal s_in(:,K)
    X = zeros(K*J,1);
    y = zeros(K*J,1);
    C = zeros(K*J,J);
    W = zeros(K*J,1);
    I = diag(ones(1,K*J));
    P = zeros(K*J,K*J);
    F = zeros(K*J,1);                                    
    for j = 1:J
        C(:,j) = [zeros(1,(j-1)*K) ones(1,K) zeros(1,J*K-j*K)].';
    end
    if(J == 1)
        FF = 1;
    else
    FF = fir1(J-1,0.99,'low',chebwin(J,30)).';
    end
    FF = [1; zeros(J-1,1)];
    % figure(111);
    % freqz(FF,1)
    F = C/(C.'*C)*FF;                           %!!!C* inv(C.'*C)*FF;
    P = I - C/(C.'*C)*C.';
    W = F;
    for iter = 1:length(s_in)
        X(K*J-K+1:K*J) = [];
        X = [s_in(iter, :).'; X];
        y(iter) = W.'*X;
        W = P*(W - mu*y(iter)*X) + F;
%         Constrains(iter) = sum(ceil(C.'*W*10000) ~= ceil(FF*10000));
%         если Constrains = 0, то ограничения выполняются; умножение на
%         10000 и округление чтобы ограничить точность матлаба до 5 знака
    end
end

