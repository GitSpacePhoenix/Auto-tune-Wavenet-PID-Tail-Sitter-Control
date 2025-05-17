function [Gbest, Pbest0, E1_best, total_error]=WPSO_PID2(P_para, K, M, N, MaxIter, yr, C_IIR, D_IIR, Z, YE, u, act_function, epsilon)
iter=1; Np=([size(P_para)]*[1 0]')/K; 
X=P_para; Kc=2;
persistent  partic Pbest V Vpre Vgbest Vpbest YE_eval YE_eval_best;
if isempty(Pbest)
    partic=1; Pbest=X; Vpre=zeros(size(X)); Vpbest=zeros(size(X)); Vgbest=zeros(size(X)); V=zeros(size(X));
    %YE_eval=repmat(YE,Np,1); YE_eval_best=repmat(YE,Np,1); %repmat(YE,Np,1);
    YE_eval=zeros(Np*M, D_IIR); YE_eval_best=zeros(Np*M, D_IIR);
end
%Añadidos provisionales
%Cmax=[1.9 1.8 2.0]; repmant(YE,Np,1)
Cmin=[0.4 0.0 0.003]; Cmax=[1.2 3.2 1.0];
maxp = max(K, max(C_IIR, D_IIR));
u_pre = u;
%Cmin_PID=[0.0 0.00 0.00]; Cmax_PID=[0.0 0.0 0.0];
% M3emoria de datos en entrada 
u_mem = zeros(Np, N);
u_mem_best = zeros(Np, N);
Gbest=Pbest(maxp*(partic-1)+1:maxp*partic,:);
%[E_extended,   E1]=wavenet_eval3(u, M, N, K, P_para, yr, C_IIR, D_IIR, Z, YE, epsilon, u_pre, y_ref0);
%%%While line%% 
while iter<MaxIter
W=[((Cmin(3)-Cmax(3))/(MaxIter-1))*(iter-1)]+Cmax(3);
C1=[((Cmin(1)-Cmax(1))/MaxIter)*iter]+Cmax(1); %Actualización de la primer constante de acelaracion
C2=[((Cmax(2)-Cmin(2))/MaxIter)*iter]+Cmin(2);
Np=([size(X)]*[1 0]')/maxp; 
for s=1:Np
    %-----------Datos de actualización de la red-------------%
    Vpre(maxp*(s-1)+1:maxp*s,:)=((2*W*rand-W)*V(maxp*(s-1)+1:maxp*s,:)); %2*W*rand-W
    Vpbest(maxp*(s-1)+1:maxp*s,:)=(2*C1*rand-C1)*(Pbest(maxp*(s-1)+1:maxp*s,:)-X(maxp*(s-1)+1:maxp*s,:)); %2*C1*rand-C1
    Vgbest(maxp*(s-1)+1:maxp*s,:)=(2*C2*rand-C2)*(Gbest-X(maxp*(s-1)+1:maxp*s,:)); %2*C2*rand-C2
    V(maxp*(s-1)+1:maxp*s,:)=Vpre(maxp*(s-1)+1:maxp*s,:)+Vpbest(maxp*(s-1)+1:maxp*s,:)+Vgbest(maxp*(s-1)+1:maxp*s,:);
    %V(K*(s-1)+1:K*s,1:2*M+N+C_IIR+D_IIR)=((2*W*rand-W)*V(K*(s-1)+1:K*s,1:2*M+N+C_IIR+D_IIR))+(2*C1*rand-C1)*(Pbest(K*(s-1)+1:K*s,1:2*M+N+C_IIR+D_IIR)-X(K*(s-1)+1:K*s,1:2*M+N+C_IIR+D_IIR))+(2*C2*rand-C2)*(Gbest(:,1:2*M+N+C_IIR+D_IIR)-X(K*(s-1)+1:K*s,1:2*M+N+C_IIR+D_IIR));
    X(maxp*(s-1)+1:maxp*s,:)=X(maxp*(s-1)+1:maxp*s,:)+V(maxp*(s-1)+1:maxp*s,:);
    %Actuales
    a=X(maxp*(s-1)+1:maxp*s-(maxp-K),1);  %Matriz de asociado la particula correspondiente
    b=X(maxp*(s-1)+1:maxp*s-(maxp-K),2); %Mattriz de b asociado a la particula correspondiente
    w=(X(maxp*(s-1)+1:maxp*s-(maxp-K),3:M+2)).'; % Quiar lo que posiblemente se añadió
    c=X(maxp*(s-1)+1:maxp*s-(maxp-C_IIR),M+3:M+	3+(M-1))';
    d=X(maxp*(s-1)+1:maxp*s-(maxp-D_IIR), M + 4 + (M-1):M+3+(M-1)+M)';
    p=X(maxp*(s-1)+1:maxp*s-(maxp-Kc),M+3+(M-1) + M+1:M+3+(M-1) + M+2);
    i=X(maxp*(s-1)+1:maxp*s-(maxp-Kc),M+3+(M-1) + M+2+1:M+3+(M-1) + M+2*2);
    dc=X(maxp*(s-1)+1:maxp*s-(maxp-Kc),M+3+(M-1) + M+2*2+1:M+3+(M-1) + M+3*2);
    %Mejores
        a_best=Pbest(maxp*(s-1)+1:maxp*s-(maxp-K),1);  %Matriz de asociado la particula correspondiente
        b_best=Pbest(maxp*(s-1)+1:maxp*s-(maxp-K),2); %Mattriz de b asociado a la particula correspondiente
        w_best=(Pbest(maxp*(s-1)+1:maxp*s-(maxp-K),3:M+2)).'; % Quiar lo que posiblemente se añadió
        c_best=Pbest(maxp*(s-1)+1:maxp*s-(maxp-C_IIR),M+3:M+	3+(M-1))';
        d_best=Pbest(maxp*(s-1)+1:maxp*s-(maxp-D_IIR), M + 4 + (M-1):M+3+(M-1)+M)';
        p_best=Pbest(maxp*(s-1)+1:maxp*s-(maxp-Kc),M+3+(M-1) + M+1:M+3+(M-1) + M+2);
        i_best=Pbest(maxp*(s-1)+1:maxp*s-(maxp-Kc),M+3+(M-1) + M+2+1:M+3+(M-1) + M+2*2);
        dc_best=Pbest(maxp*(s-1)+1:maxp*s-(maxp-Kc),M+3+(M-1) + M+2*2+1:M+3+(M-1) + M+3*2);

     for k=1:K
                 %for m=1:M
                     Tau(k)=sqrt((u(1)-b(k))^2+(u(2)-b(k))^2+(u(3)-b(k))^2+(u(4)-b(k))^2)/a(k);
                     Tau_best(k)=sqrt((u(1)-b_best(k))^2+(u(2)-b_best(k))^2+(u(3)-b_best(k))^2+(u(4)-b_best(k))^2)/a_best(k);
                     switch act_function
                        case 'morlet'
                            psi(k)=cos(0.5*Tau(k))*exp(-0.5*Tau(k)^2); %Morlet
                            psi_best(k)=cos(0.5*Tau_best(k))*exp(-0.5*Tau_best(k)^2); %Morlet_best
                        case 'rasp1'
                              psi(k)=(Tau(k))/(Tau(k)^2+1)^2; %Rasp1 
                              psi_best(k)=(Tau_best(k))/(Tau_best(k)^2+1)^2; %Rasp1 best
                         case 'rasp2'
                             psi(k)=(Tau(k))*cos(Tau(k))/(Tau(k)^2+1); %Rasp2 
                             psi_best(k)=(Tau_best(k))*cos(Tau_best(k))/(Tau_best(k)^2+1); %Rasp2 best
                         case 'slog1'
                                psi(k)=1/(1+exp(-Tau(k)+1)) - 1/(1+exp(-Tau(k)+3)) - 1/(1+exp(-Tau(k)-3)) + 1/(1+exp(-Tau(k)-1));
                                psi_best(k)=1/(1+exp(-Tau_best(k)+1)) - 1/(1+exp(-Tau_best(k)+3)) - 1/(1+exp(-Tau_best(k)-3)) + 1/(1+exp(-Tau_best(k)-1));
                        case 'slog2'
                                psi(k)=3/(1+exp(-Tau(k)-1)) - 3/(1+exp(-Tau(k)+1)) - 3/(1+exp(-Tau(k)-3)) + 3/(1+exp(-Tau(k)-3));
                                psi_best(k)=3/(1+exp(-Tau_best(k)-1)) - 3/(1+exp(-Tau_best(k)+1)) - 3/(1+exp(-Tau_best(k)-3)) + 3/(1+exp(-Tau_best(k)-3));
                        case 'poliwog1'
                                psi(k)=Tau(k)*exp(-(1/2)*Tau(k)^2); %Polywog1
                                psi_best(k)=Tau_best(k)*exp(-(1/2)*Tau_best(k)^2); %Polywog1best
                        case 'poliwog2'
                                psi(k)=(Tau(k)^3 - 3*Tau(k))*exp(-(1/2)*Tau(k)^2); %Polywog2
                                psi_best(k)=(Tau_best(k)^3 - 3*Tau_best(k))*exp(-(1/2)*Tau_best(k)^2); %Polywog2best
                        case 'poliwog3'
                                psi(k)=(Tau(k)^4 -6*Tau(k)^2 + 3)*exp(-(1/2)*Tau(k)^2); %Polywog3
                                psi_best(k)=(Tau_best(k)^4 -6*Tau_best(k)^2 + 3)*exp(-(1/2)*Tau_best(k)^2); %Polywog3best
       
                         otherwise
                                psi(k) = Tau(k);     
                                psi_best(k) = Tau_best(k);     
                    end
                      
                      %Psi(k) = (1/sqrt(a(k)))*psi(k);
                      %Psi_best(k) = (1/sqrt(a_best(k)))*psi_best(k);
                 %end
              end
             z_hat(s,:)=psi*w'; z_hat_best(s,:)=psi_best*w_best';
             
             %Z_eval(M*(s-1)+1:M*s,2)=Z(:,1); Z_eval_best(M*(s-1)+1:M*s,2)=Z(:,1); %Resultados de instante anmteior en segunda columna
            for j=0:C_IIR-2
                Z_eval(M*(s-1)+1:M*s,C_IIR-j)=Z(:,C_IIR-j-1);  %Guardado de instantes enteriores para el filtro   
                Z_eval_best(M*(s-1)+1:M*s,C_IIR-j)=Z(:,C_IIR-j-1); %Guardado de mejores instantes anteiores
            end
            Z_eval(M*(s-1)+1:M*s,1)=z_hat(s,1:M)'; Z_eval_best(M*(s-1)+1:M*s,1)=z_hat_best(s,1:M)';
            %Cambio de   ,2)=YE(:,1) por ,:)=YE
            YE_eval(M*(s-1)+1:M*s,:)=YE; YE_eval_best(M*(s-1)+1:M*s,:)=YE;
            %Cada que aparezca un Z o un YE actual habra un Z best y YE best
            %respectuvamente 7pr7caur
            %uc=PID_discrete(p, i, dc, u, epsilon);
            %uc_best=PID_discrete(p_best, i_best, dc_best, u, epsilon);
            uc=PID_discrete(p, i, dc, u, epsilon);
            uc_best=PID_discrete(p_best, i_best, dc_best, u, epsilon);
            u_mem(s,:) = uc;
            u_mem_best(s,:) = uc_best;
            %v=(uc./2); v_best=(uc_best./2);
            v = 3; v_best = 3;
            %for j=1:N 
            %    %Fue cambiado el (s,:) por (N*(s-1)+j,:)
            %    yem(s,j)=(Z_eval(N*(s-1)+j,:)*c(j,:)')*(uc*ones(M,1));  yem_best(s,j)=(Z_eval_best(N*(s-1)+j,:)*c_best(j,:)')*(uc_best*ones(M,1));
            %    yen(s,j)=(YE_eval(N*(s-1)+j,:)*d(j,:)')*(v*ones(M,1));  yen_best(s,j)=(YE_eval_best(N*(s-1)+j,:)*d_best(j,:)')*(v_best*ones(M,1));
            %end %Cada instante actual tiene su Pbest
            [yem_sum, yen_sum] = filter_fun2(c, d, s, N, Z_eval, YE_eval);
            [yem_sum_best, yen_sum_best] = filter_fun2(c_best, d_best, s, N, Z_eval_best, YE_eval_best);
            yem(s,:) = yem_sum*sum(uc); yem_best(s,:) = yem_sum_best*sum(uc_best);
           yen(s,:) = yen_sum*sum(v); yen_best(s,:) = yen_sum_best*sum(v_best);
            y_hat_eval(s,:)=yem(s,:)+yen(s,:);  y_hat_eval_best(s,:)=yem_best(s,:)+yen_best(s,:);
%             if iter==MaxIter-1 && s==Np
%             disp('YE');
%             end
            for j=0:D_IIR-2
                YE_eval(M*(s-1)+1:M*s,D_IIR-j)=YE(:,D_IIR-j-1);  %Guardado de instantes enteriores para el filtro   
                YE_eval_best(M* (s-1)+1:M*s,D_IIR-j)=YE(:,D_IIR-j-1); %Guardado de mejores instantes anteiores
            end
%             if iter==MaxIter-1 && s==Np
%             disp('YE');
%             end
            YE_eval(M*(s-1)+1:M*s,1)=y_hat_eval(s,:); YE_eval_best(M*(s-1)+1:M*s,1)=y_hat_eval_best(s,:); %Guardado de instantes actuales
        
    error=(y_hat_eval(s,:)-yr); error_best=(y_hat_eval_best(s,:)-    yr);
    E1(1,s)=(1/2)*(error(1).*error(1))+(1/2)*(error(2).*error(2))+(1/2)*(error(3).*error(3))+(1/2)*(error(4).*error(4));
    E1_best(1,s)=(1/2)*(error_best(1).*error_best(1))+(1/2)*(error_best(2).*error_best(2))+(1/2)*(error_best(3).*error_best(3))+(1/2)*(error_best(4).*error_best(4)); 
    if E1(1,s)<E1_best(1,s)
        Pbest(maxp*(s-1)+1:maxp*s,:)=X(maxp*(s-1)+1:maxp*s,:);
    end

      %mejores partículas F:\Diseño estructural\Dibijos an
      %SC\Costillas2\frame1 maxp*(partic-1)+1:maxp*partic
    end
[value, partic]=min(E1_best);
Gbest=Pbest(maxp*(partic-1)+1:maxp*partic,:);
 Pbest0=Pbest;
    if mod(iter,3)==0
        X(maxp*Np+1:(maxp*Np+(6*maxp)),:)=poblacion_inicial3_2(K, N, M, 6, C_IIR, D_IIR); %Segundo incremento de población
        %X(maxp*Np+1:(maxp*Np+(4*maxp)),:) = Gbest_dispertion(Gbest, W, 4, M, maxp);
        V(maxp*Np+1:(maxp*Np+(6*maxp)),:)=zeros(size(X(maxp*Np+1:(maxp*Np+(6*maxp)),:)));
        Pbest(maxp*Np+1:(maxp*Np+(6*maxp)),:)=X(maxp*Np+1:(maxp*Np+(6*maxp)),:);
        YE_eval(M*Np+1:M*Np+6*M,1:D_IIR)=zeros(6*M,D_IIR); YE_eval_best(M*Np+1:M*Np+6*M,1:D_IIR)=zeros(6*M,D_IIR);
    end
iter=iter+1;
end
total_error=(y_hat_eval_best(partic,:)-yr)*(y_hat_eval_best(partic,:)-yr)';
end