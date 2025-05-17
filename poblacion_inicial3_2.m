function P_para=poblacion_inicial3_2(k_neuron, N, m, P_paraartNumber, C_IIR, D_IIR)
% Neuronas(k)  entradas(N) salidas(m) Particulas (P_paraartNumber) a(kxm)    b(kxm)   w(nxk)
%     k=2;        N=3;          m=1;       P_paraartNumber=10
%Numero de P_paraarticulas de P_paraSO                          [iter E1_best(1,partic) E_extended_best(1,partic) yr]
Coef_a= [3 3 3 3 3 3 3 3 3 3 3 3 3 3 8 8 8 8 8 8 8 8 8];  %Vectores de dilatacion
Coef_b= [2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 2 2 2 5 5 5 5];  % Estos acotan los límites del npumero aleatorio
Coef_w= [2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2];  %
Coef_c= [5 5 5 5 5 5 5 5 5 5 5 5 5 3 8 8 8 8 8 8 5 5 5];  %    
Coef_d= [4 4 4 4 4 4 4 4 4 4 4 8 8 8 8 8 8 8 8 8 5 5 5];  %    
Coef_Kp1=120*[14 14 14 14 14 14 14 16 16 11 4 9 9 9 9 8 8 8 8 8 5 5 5];
Coef_Kp2=11*[14 14 14 14 14 14 14 14 14 14 14 14 9 9 9 8 8 8 8 8 5 5 5];
Coef_Ki1=-70.8*[5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 8 8 8 8 5 5 5];
Coef_Ki2=-70.8*[5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 8 8 8 8 5 5 5];
Coef_Kd1=5*[1.3 1.3 1.1 1.1 0.9 0.9 0.9 0.9 0.9 0.9 0.9 8 8 8 8 8 8 8 8 8 5 5 5];
Coef_Kd2=3*[2 2 2 2 2 2 2 2 2 2 2 2 2 2 28 28 8 8 5 5 5];
k=max(k_neuron,max(C_IIR, D_IIR));
P_para=zeros(P_paraartNumber*k, 2*m+N); n=1; da=1; db=1; dw=1; dc=1; dd=1; dKp=1; dKi=1; dKd=1;


%%%%%%%%%%%%%%%%Diataciones%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while n<=P_paraartNumber*k %Salto de bloque a bloque en A
    c=da*k; %Incremento P_paraara el ciclo  
    coef1=Coef_a(c/k); %Coeficiente correspondiente al vector de dilatciones 
    for i=n:c %Llenado vertical de datos por cada neurona hasta completar todas en la población
        %for j=1:m %llenadp horizontal de datos segun el numero de entradas
            P_para(i,1)=-coef1+2*coef1*rand; % -5 + (5+5)*rand || Selección aleatoria de cada dato
        %end
    end
    da=da+1;
    
    %P_paraesos%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cb=db*k; %Incremento P_paraara el%%%%%%%%%%%%%%%%%%% ciclo  
    coef1=Coef_b(cb/k);
    for i=n:cb
        %for j=m+1:2*m
            P_para(i,2)=-coef1+2*coef1*rand;
        %end
    end
    db=db+1;
    %%%%%%%%%%%%%%%Escalamientos%%%%%%%%%%%%%%%%%%%%
    cw=dw*k; %Incremento P_paraara el ciclo  
    coef1=Coef_w(cw/k);
    for i=n:cw
        for j=3: m+2
            P_para(i,j)=-coef1+2*coef1*rand;
        end
    end 
    dw=dw+1;
    %%%%%%%%%%%%%%Adelanto y atraso%%%%%%%%%%%%%%%%%%%%%%%%
    cc=dc*k;
    coef1=Coef_c(cc/k);
    for i=n:cc
       for j=m+3:m+	3+(m-1)
            P_para(i,j)=-coef1+2*coef1*rand;
       end
    end
    dc=dc+1;
    %%%%%%%%%%%Coeficientes de atraso%%%%%%%%%%%%%%%%%%%%%
    cd=dd*k;
    coef1=Coef_d(cd/k);
    for i=n:cd
       for j=m + 4 + (m-1):m+3+(m-1)+m
          P_para(i,j)=-coef1+2*coef1*rand; 
       end
    end
    dd=dd+1;
     %%%%%%%%% Coeficientes proporcionales %%%%%%%%%%%%%%%
     cKp=dKp*k;
     coef1=Coef_Kp1(cKp/k);
     for i=n:cKp
         if rem(i,2) == 1
            coef1 = Coef_Kp1(cKp/k); 
         else
             coef1 = Coef_Kp2(cKp/k);
         end
         
         for j=m+3+(m-1) + m+1:m+3+(m-1) + m+2
            %P_para(i,j)=-coef1+2*coef1*rand; 
            P_para(i,j)=-coef1*rand; 
         end
     end
     
     dKp=dKp+1;
%     %%%%%%% Coeficientes integrativos %%%%%%%%%%%%%%%%%%
     cKi=dKi*k;
     
     for i=n:cKi
         for j=m+3+(m-1) + m+2+1:m+3+(m-1) + m+2*2
             if rem(i,2) == 1
                coef1 =  Coef_Ki1(cKi/k);
             else
                 coef1 =  Coef_Ki2(cKi/k);
             end
             %P_para(i,j)=-coef1+2*coef1*rand; 
             P_para(i,j)=coef1*rand;
         end
    end
     dKi=dKi+1;
     %     %%%%%%% Coeficiente derivativo %%%%%%%%%%%%%%%%%%%%
     cKd=dKd*k;
     for i=n:cKd
         if rem(i,2)==1
             coef1=Coef_Kd1(cKd/k);
         else
             coef1=Coef_Kd2(cKd/k);
         end
         for j=m+3+(m-1) + m+2*2+1:m+3+(m-1) + m+3*2
            %P_para(i,j)=-coef1+2*coef1*rand; 
            P_para(i,j)=coef1*rand; 
         end
     end
     dKd=dKd+1;
     
n=n+k; %Incremento de 3 en tres P_paraarra recorrrido vertical
end
for  i=1:P_paraartNumber
   P(i,1)=i;
end
P_para;