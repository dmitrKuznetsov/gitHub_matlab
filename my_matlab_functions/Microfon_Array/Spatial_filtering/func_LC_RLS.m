function [y, h] = func_LC_RLS(s_in,K,J,lamda)

    x = zeros(K*J,1);
    y = zeros(length(s_in),1);
    alpha = zeros(length(s_in),1);
    C = zeros(K*J,J);
    I = diag(ones(1,K*J));
    for j = 1:J
        C(:,j) = [zeros(1,(j-1)*K) ones(1,K) zeros(1,J*K-j*K)].';
    end
%     if(J == 1)
%         FF = 1;
%     else
%         FF = fir1(J-1,0.99,'low',chebwin(J,30)).';
%     end
    FF = [1; zeros(J-1,1)];



    R_inv = 1000*I; 
    Gamma = R_inv*C; 
    Psi_inv = inv(C'*Gamma);
    h = Gamma*Psi_inv*FF;

    for iter = 1:length(s_in)
        x(K*J-K+1:K*J) = [];
        x = [s_in(iter, :).'; x];
        g = R_inv*x/(   lamda + x'*R_inv*x   );
        R_inv = 1/lamda*(   R_inv - g*x'*R_inv   );

        v = C'*g;
        nu = x'*Gamma;
        Ij = Psi_inv*v/(   1 - nu*Psi_inv*v  );
        Psi_inv = lamda*(Psi_inv+Ij*nu*Psi_inv);
        Gamma = 1/lamda*(Gamma-g*nu);
        y(iter) = h'*x;
        alpha(iter) = 0 - y(iter);
        h = h + (g-lamda*Gamma*Ij)*conj(alpha(iter));
    end

end
































