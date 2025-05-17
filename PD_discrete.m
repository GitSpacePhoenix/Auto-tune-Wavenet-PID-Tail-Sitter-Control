function u=PD_discrete(Kp, Kd, u_pre, epsilon)
       Ts=0.01; 
       Td=[Kd(1,1)/Kp(1,1) Kd(1,1)/Kp(1,1); Kd(2,1)/Kp(2,1) Kd(2,1)/Kp(2,1)];
       K0=[Kp(1,1)*(1+Td(1,1)/Ts) Kp(1,1)*(1+Td(1,1)/Ts) Kp(2,1)*(1+Td(2,1)/Ts) Kp(2,1)*(1+Td(2,1)/Ts)];
       K1=[-Kp(1,1)*(1+2*(Td(1,1)/Ts)) -Kp(1,1)*(1+2*(Td(1,1)/Ts)) -Kp(2,1)*(1+2*(Td(2,1)/Ts)) -Kp(2,1)*(1+2*(Td(2,1)/Ts))];
       K2=[Kp(1,1)*Td(1,1)/Ts Kp(1,1)*Td(1,1)/Ts Kp(2,1)*Td(2,1)/Ts Kp(2,1)*Td(2,1)/Ts];
       out=u_pre+K0.*epsilon(1,:)+K1.*epsilon(2,:)+K2.*epsilon(3,:);
       u(1:2)=min(0.9, max(0.08, out(1:2))); u(3:4)=min(45, max(-45, out(3:4)));
end