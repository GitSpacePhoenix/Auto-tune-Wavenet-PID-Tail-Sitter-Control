function [y_and_u]=wavenetMIMO_IIR_v5_2(x)
    K = 2; act_function = 'poliwog2'; N=4; M=4; MaxIter=9; Np=13; C_IIR=2; D_IIR=2; K2=2;
    %x(2) --> z_ref  x(4) --> theta_ref x(6) --> z   x(7) --> theta
    y_ref0(1)=x(3); y_ref0(2)=x(3); y_ref0(3)=x(2); y_ref0(4)=x(2); y(1)=x(5); y(2)=x(5); y(3)=x(6); y(4)=x(6); time=floor(100*x(9)+1);
      %[u, z, theta, phi]    K = 6; act_function = 'poliwog3'; N=4; M=4; MaxIter=13; Np=13; C_IIR=2; D_IIR=2; K2=2;
    persistent ti; persistent Z; persistent YE; persistent epsilon Global_epsilon; persistent u_pre E_func; 
    persistent P_para; persistent u Global_u; persistent Gbest Global_PID Gbest_PID Global_Gbest PID_ctrl error;
    yr=[x(5) x(6) x(7) x(8)];
    if isempty(P_para)
       ti=8; YE=zeros(N,D_IIR); YE(:,1)=yr; Z=zeros(N,C_IIR); epsilon=zeros(M,M); error=2; u_pre=0.2*ones(1,M); u=0.2*ones(1,N); Gbest=zeros(K,2*M+N+C_IIR+D_IIR+6); 
       P_para=poblacion_inicial3_2(K,N, M, Np, C_IIR, D_IIR); 
        end
    %--------Inicialización de matricesPr
    
    %u(1,1:2)=(0.32*sin(0.05*time)+0.6679 +0.1*(rand-rand))*ones(1,2);
    %u(1,3:4)=0*ones(1,2);
    for i=0:M-2 %Registro de error de relulación (control) en instantes anteriores
       epsilon(M-i,:)=epsilon(M-i-1,:);
    end
   maxp = max(K, max(C_IIR, D_IIR));
    %epsilon(1,1:2)=y_ref0(1,1:2)-[x(6) x(6)];     epsilon(1,3:4)=[x(7) x(7)]-y_ref0(1,3:4);
    epsilon(1,1:2)=[x(5) x(5)]-y_ref0(1,1:2);     epsilon(1,3:4)=y_ref0(1,3:4)-[x(6) x(6)];
    % if abs(error)>0.02
     [Gbest, Pbest, E1, error]=WPSO_PID2(P_para, K, M, N, MaxIter, yr, C_IIR, D_IIR, Z, YE, u_pre, act_function, epsilon); 
     P_para=order_decent(Pbest, maxp, M, N, E1, C_IIR, D_IIR, Np); %Función que re-ordena la selección elitista
    % end
    a=Gbest(1:maxp-(maxp-K),1); b=Gbest(1:maxp-(maxp-K),2); w=Gbest(1:maxp-(maxp-K),3:M+2).';         
    c=Gbest(1:maxp-(maxp-C_IIR),M+3:M+	3+(M-1))';  d=Gbest(1:maxp-(maxp-D_IIR),M + 4 + (M-1):M+3+(M-1)+M)';
    p=Gbest(1:maxp-(maxp-K2),M+3+(M-1) + M+1:M+3+(M-1) + M+2);
    i=Gbest(1:maxp-(maxp-K2),M+3+(M-1)+M+2+1:M+3+(M-1) + M+2*2);
    dc=Gbest(1:maxp-(maxp-K2),M+3+(M-1) + M+2*2+1:M+3+(M-1) + M+3*2);
    Global_Gbest(:,:,time) = Gbest;
    %Tau = zeros(K,1); psi = zeros (K,1);
    Global_epsilon(time,:) = epsilon(1,:);
    
    Global_u(time,:) = u;
    %u=PID_discrete(Kp, Ki, Kd, u_pre, epsilon);
    for k=1:K
            Tau(k)=sqrt((u(1)-b(k))^2+(u(2)-b(k))^2+(u(3)-b(k))^2+(u(4)-b(k))^2)/a(k);
            switch act_function
                case 'morlet'
                    psi(k)=cos(0.5*Tau(k))*exp(-0.5*Tau(k)^2); %Morlet
                case 'rasp1'
                    psi(k)=(Tau(k))/(Tau(k)^2+1)^2; %Rasp1
                case 'rasp2'
                    psi(k)=(Tau(k)*cos(Tau(k)))/(Tau(k)^2+1); %Rasp2
                case 'slog1'
                    psi(k)=1/(1+exp(-Tau(k)+1)) - 1/(1+exp(-Tau(k)+3)) - 1/(1+exp(-Tau(k)-3)) + 1/(1+exp(-Tau(k)-1));
                case 'slog2'
                    psi(k)=3/(1+exp(-Tau(k)-1)) - 3/(1+exp(-Tau(k)+1)) - 3/(1+exp(-Tau(k)-3)) + 3/(1+exp(-Tau(k)-3));
                case 'poliwog1'
                    psi(k)=Tau(k)*exp(-(1/2)*Tau(k)^2); %Polywog1
                case 'poliwog2'
                    psi(k)=(Tau(k)^3 - 3*Tau(k))*exp(-(1/2)*Tau(k)^2); %Polywog2
                case 'poliwog3'
                    psi(k)=(Tau(k)^4 -6*Tau(k)^2 + 3)*exp(-(1/2)*Tau(k)^2); %Polywog3
                otherwise
                    psi(k) = Tau(k);     
            end
            %Psi(k) = (1/sqrt(a(k)))*psi(k);
    end
    %assignin('base','u',u); assignin('base','epsilon',epsilon);
    z_hat=psi*w'; % señal re-construida
    u=PID_discrete(p, i, dc, u_pre, epsilon);
    %u = PD_discrete(p, dc, u, epsilon);
    for i=0:C_IIR-2
        Z(:,C_IIR-i)=Z(:,C_IIR-i-1);
    end
     Z(:,1)=z_hat'; v=3; 
    Global_u(time,:) = u;
    %for i=1:N
    %    yem(1,i)=(Z(i,:)*c(i,:)')*(u*ones(M,1));
    %    yen(1,i)=(YE(i,:)*d(i,:)')*(v*ones(M,1));
    %end
    [yem_sum, yen_sum] = filter_fun2(c, d, 1, N, Z, YE); % s = 1 en este caso
    yem = yem_sum*sum(u);
    yen = yen_sum.*sum(v);
    y_hat=yem+yen;
    errorf = y_hat'-yr;
    E_func(time) = (1/2)*(errorf*errorf');
    error = E_func(time);
    y_and_u(1,1:M)=y_hat; y_and_u(1,M+1:2*M)=u;
    u_pre=u;
    %if time>(D_IIR-1)
        for i=0:D_IIR-2
           YE(:,D_IIR-i)=YE(:,D_IIR-i-1);
        end  
    %end
    YE(:,1)=y_hat';
    
     if time == 1280
         assignin('base','params',Global_Gbest);
         assignin('base','epsilon_plot',Global_epsilon);
         assignin('base','u_plot',Global_u);
         assignin('base','func_error',E_func);
     end
end