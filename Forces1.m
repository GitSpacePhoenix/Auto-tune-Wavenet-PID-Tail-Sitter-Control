function [F, Tau]=Forces1(VT, alpha, beta, theta, phi, psi, Tr1, Tr2, Sr, Sl)
%alpha=0; beta=0;
%%%%%%%%%%%%%%%%%ASUMIENDO  HOVER%%%%%%%%%%%%%%%%%
f1=Motor(Tr1); f2=Motor(Tr2); l=0.2;%syms f1 f2 TauMi1 TauMi2 Sr Sl; %Parámetros de motores
u=f1(1)+f2(1);   %Vector de fuerzas de motores
F_barMotors=[u; 0; 0];  %Fuerzas que contribuyen a la tracción en el eje z
Tau_Motors=[(f2(1)-f1(1))*l;0;f1(2)-f2(2)];                                                                       
 %%%%%%%%%%%%%ASUMIENDO MARCO VUELO CRUCERO%%%%%%%%%%%%%%%%%%%
MAC=0.22; q=1/2*0.9691*(VT^2); S=0.00493782; b=0.2;
% q=120; S=0.460; MAC=26; b=1.11; %Constantes físicas
[CL, CD, CY, Cl, Cm, Cn]=Aero_Coefs(alpha, beta, Sr, Sl); %Función donde se calculan los coeficientes aerodinámicos
Lift=(CL)*q*S; %Fuerza de levantamiento
Drag=(CD)*q*S; %Fuerza de arrastre
Fay=(CY)*q*S;  %Fuerza lateral
F_barWing=[(-Drag*cos(deg2rad(alpha))+Lift*sin(deg2rad(alpha))); Fay; (Drag*sin(deg2rad(alpha))+Lift*cos(deg2rad(alpha)))];     %Aplicando nueva convención de ejes
assignin('base','F_barMotors',F_barMotors); assignin('base','F_barWing',F_barWing);
Lbar=(Cl)*q*S*b;
Mbar=(Cm)*q*S*MAC;
Nbar=(Cn)*q*S*b;
Tau_Wing=[Lbar; Mbar; Nbar];
assignin('base','Tau_Motors',Tau_Motors); assignin('base','Tau_Wing',Tau_Wing);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%VALORES TOTALES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Suma de fuerzas y momentos totales sin rotación
F_tot=F_barWing+F_barMotors; Tau=Tau_Motors+Tau_Wing;
 R_be=[cos(theta)*cos(psi) cos(theta)*sin(psi) -sin(theta);
     -cos(phi)*sin(psi)+sin(phi)*sin(theta)*cos(psi) cos(phi)*cos(psi)+sin(phi)*sin(theta)*sin(psi) sin(phi)*cos(theta);
     sin(phi)*sin(psi)+cos(phi)*sin(theta)*cos(psi) -sin(phi)*cos(psi)+cos(phi)*sin(theta)*sin(psi) cos(phi)*cos(theta)];
F=R_be*F_tot; %Fueerzas totales rotadas