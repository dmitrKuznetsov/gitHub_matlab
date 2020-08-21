function [y,W_out] = my_nlms_Frost_2sposob_2mic_imitac(s_in, J, mu)


    X1 = zeros(J,1);
    X2 = zeros(J,1);
    y = zeros(length(s_in),1);
    
    C = [1; zeros(J-1,1)];

    W1 = C/2;
    W2 = C/2;


    for iter = 1:length(s_in)
        if(iter == 11)
            stop_bit = 1;
        end
        X1 = [s_in(iter, 1); X1(1:end-1)];
        X2 = [s_in(iter, 2); X2(1:end-1)];
        
        sum1 = 0;
        sum2 = 0;
        for tapCnt = 1 : J
            sum1 = sum1 + W1(tapCnt)*X1(tapCnt);
            sum2 = sum2 + W2(tapCnt)*X2(tapCnt);
        end
        y(iter) = sum1 + sum2;
        err = - (sum1 + sum2);
        energy = X1.'*X1 + X2.'*X2 + 0.000000119209289;
        w = (err*mu)/energy;
        
        for tapCnt = 1 : J
            W1(tapCnt) = W1(tapCnt) + w*X1(tapCnt);
            W2(tapCnt) = W2(tapCnt) + w*X2(tapCnt);
        end

        for tapCnt = 1 : J
            v = (C(tapCnt) - W1(tapCnt) - W2(tapCnt))/2;
            W1(tapCnt) = W1(tapCnt) + v;
            W2(tapCnt) = W2(tapCnt) + v;
        end

    end
    
    W_out = [W1; W2];


end

