function [out] = func_CIC_N4(y,R,nbit)
    int1 = 0;
    int2 = 0;
    int3 = 0;
    int4 = 0;  
    kk = 1;
    last_c1 = 0;
    last_c2 = 0;
    last_c3 = 0;
    last_int4 = 0;
    last2_c1 = 0;
    last2_c2 = 0;
    last2_c3 = 0;
    last2_int4 = 0;
    out = zeros(1,length(y)/R);

    for i=1:length(y)
      int1 = int1 + y(i) ;
      int1 = mod2n(int1,nbit);

      
      int2 = int2 + int1 ;
      int2 = mod2n(int2,nbit);


      int3 = int3 + int2 ;
      int3 = mod2n(int3,nbit);


      int4 = int4 + int3 ;
      int4 = mod2n(int4,nbit);


      if (mod(i,R)==0)
        comb1 = int4 - last2_int4 ;
        comb1 = mod2n(comb1,nbit);


        comb2 = comb1 - last2_c1 ;
        comb2 = mod2n(comb2,nbit);


        comb3 = comb2 - last2_c2 ;
        comb3 = mod2n(comb3,nbit);


        comb4 = comb3 - last2_c3 ;
%         mod2n(comb4,nbit)
        out(kk) = mod2n(comb4,nbit);
        last2_int4 = last_int4 ;
        last_int4 = int4 ;
        last2_c1 = last_c1 ;
        last_c1 = comb1;
        last2_c2 = last_c2 ;
        last_c2 = comb2;
        last2_c3 = last_c3 ;
        last_c3 = comb3;
        kk = kk + 1 ;
      end
    end
end

