function [y,W] = my_nlms_Frost_2sposob(s_in, J, K)
%   frost1972.pdf
%   J - filter order
%   K - Number of microphones
%   s_in - input signal s_in(:,K)

mu = 1;
W = zeros(K,J);
X = zeros(K,J);
y = zeros(length(s_in),1);
if(J == 1)
    C = 1;
else
    C = fir1(J-1,0.99,'low',chebwin(J,30)).';
end

for ii = 1:length(W(:,1))
    W(ii,:) = 1/K*C.';
end

for iter = 1:length(s_in)
    X(:,J) = [];
    X = [s_in(iter, :).' X];
    y(iter) = sum(sum(W.*X));
    E = sum(sum(X.^2)) + 0.000000119209289;
    W = W + mu/E*(-y(iter))*X;
    
    W = W + 1/K*(C.'-ones(1,K)*W);

end

end

