function [y,W_out] = my_nlms_Frost_2sposob_2mic_imitac_single(s_in, J, mu)


    X1 = single( zeros(J,1));
    X2 = single( zeros(J,1));
    y =  single( zeros(length(s_in),1));
    
    C = single([1; zeros(J-1,1)]);

%     if(J == 1)
%         C = 1;
%     else
%         C = single( fir1(J-1,0.99,'low',chebwin(J,30)).');
%     end        
    W1 = single( C/2 );
    W2 = single( C/2 );


    for iter = 1:length(s_in)
        if(iter == 27)
            stop_bit = 1;
        end
        X1 = [s_in(iter, 1); X1(1:end-1)];
        X2 = [s_in(iter, 2); X2(1:end-1)];
        
        sum1 = single( 0 );
        sum2 = single( 0 );
        for tapCnt = 1 : J
            sum1 = sum1 + W1(tapCnt)*X1(tapCnt);
            sum2 = sum2 + W2(tapCnt)*X2(tapCnt);
        end
        y(iter) = sum1 + sum2;
        err = - (sum1 + sum2);
        energy = single( X1.'*X1 + X2.'*X2 + 0.000000119209289 );
        w = single(  (err*mu)/energy );
        
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

