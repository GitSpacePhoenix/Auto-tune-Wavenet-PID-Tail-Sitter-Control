function P_para=order_decent(Part_values, K, M, N, E1, C_IIR, D_IIR, Np)
 [E, Index]=sort(E1);
  for j=1:Np
    a=Part_values(K*(Index(j)-1)+1:K*Index(j),1);
    b=Part_values(K*(Index(j)-1)+1:K*Index(j),2);
    w=(Part_values(K*(Index(j)-1)+1:K*Index(j),3:M+2)).';
    c=Part_values(K*(Index(j)-1)+1:K*Index(j),M+3:M+3+(M-1));
    d=Part_values(K*(Index(j)-1)+1:K*Index(j),M+4+(M-1):M+3+(M-1)+M);
    p=Part_values(K*(Index(j)-1)+1:K*Index(j),M+3+(M-1)+M+1:M+3+(M-1) + M+2);
    i=Part_values(K*(Index(j)-1)+1:K*Index(j),M+3+(M-1)+M+2+1:M+3+(M-1) + M+2*2);
    dc=Part_values(K*(Index(j)-1)+1:K*Index(j),M+3+(M-1)+M+2*2+1:M+3+(M-1) + M+3*2);
    
    P_para(K*(j-1)+1:K*j,:)=[a b w.' c d p i dc];
  end
end