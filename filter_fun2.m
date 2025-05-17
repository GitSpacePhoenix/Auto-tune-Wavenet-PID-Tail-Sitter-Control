function [yem, yen] = filter_fun2(c, d, s, N, Z_eval, YE_eval)

yem = zeros (N,1);
yen = zeros (N,1);

    for j=1 : N
        yem(j) = c(j,:)* Z_eval(N*(s-1)+1+(j-1):N*(s-1)+j,:)';
    end

    for j = 1 : N
        yen(j) =  d(j,:)* YE_eval(N*(s-1)+1+(j-1):N*(s-1)+j,:)';
    end
end