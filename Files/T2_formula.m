function val = T2_formula(T2,A,TE)

fitFunc = @(T2,A,x) A*exp(-x./T2);
val = fitFunc(T2,A,TE);