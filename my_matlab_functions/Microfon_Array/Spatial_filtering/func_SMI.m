function [W] = func_SMI(X1,X2,L)
Nel = 2;
N = length(X1);

X_array = zeros(L, N-L);

for i = 0:L-1
    X_array(2*i+1,:) = X1(L-i+1:N-i).';
    X_array(2*i+2,:) = X2(L-i+1:N-i).';
end 
Rxx = 1/(N-L+1)*(X_array*X_array.');
Rxx_inv = inv(Rxx);
C = zeros(2*L,L);
for j = 1:L
    C(:,j) = [zeros(1,(j-1)*2) ones(1,2) zeros(1,L*2-j*2)].';
end

% F = fir1(L-1,0.99,'low',chebwin(L,30)).';
F = [1; zeros(L-1,1)];

W = Rxx_inv*C*inv(C.'*Rxx_inv*C)*F;

end

