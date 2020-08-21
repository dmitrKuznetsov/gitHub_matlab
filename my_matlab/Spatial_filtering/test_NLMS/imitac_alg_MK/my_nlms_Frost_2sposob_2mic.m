function [y,W_out] = my_nlms_Frost_2sposob_2mic(s_in, J,mu)


    X1 = zeros(J,1);
    X2 = zeros(J,1);
    y = zeros(length(s_in),1);

    C = [1; zeros(J-1,1)];

    W1 = 1/2*C;
    W2 = 1/2*C;


    for iter = 1:length(s_in)
        X1 = [s_in(iter, 1); X1(1:end-1)];
        X2 = [s_in(iter, 2); X2(1:end-1)];
        
        y(iter) = W1.'*X1 + W2.'*X2;
        energy = X1.'*X1 + X2.'*X2 + 0.000000119209289;
        w = mu/energy*(-y(iter));
        W1 = W1 + w*X1;
        W2 = W2 + w*X2;
        
        V = 1/2*( C - (W1 + W2) );
        W1 = W1 + V;
        W2 = W2 + V;
    end
    
    W_out = [W1; W2];


end

