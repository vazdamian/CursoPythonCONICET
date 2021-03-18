clc;                        % Clears the screen
clear all;
format long

M0005=dlmread('E:\TP_Python_CONICET\f3a_CP0005_Ox3-T26Hr86_MODIFICADO.DTA');

area=0.6*pi*3;
Rs=13000;           %ohm*cm2
Y0=1.1e-4;          %F/cm2
a=0.85;             %alfa del CPE
Icor=0.42e-9;       %A/cm2
ba=10;              %V/dec   
bc=0.09;            %V/dec
k=Y0/(gamma(1-a));
Rp=ba*bc/(2.3*Icor*(ba+bc));
f=1/(1-a);

Ec0005=-0.1020515;
Itexp0005=M0005(:,4)'/area;
Vtexp0005=M0005(:,3)'-Ec0005;
Vfexp0005=Vtexp0005-Rs*Itexp0005;

sr=0.000005;            % Polarization scanning rate V/sec
Vapex=-0.3;             % Potencial Total máximo  del vértice de la polarización triangular 
h=8;                    % Step size sec.
tapex=-Vapex/sr+2*h;    % Tiempo al Potencial Total máximo o del vértice de la polarización triangular

tida=0:h:tapex; Vida=-sr*(tida-2*h); Vida(1)=0; Vida(2)=0;
tvuelta=tapex+h:h:2*tapex-2*h; Vvuelta=2*Vapex+sr*(tvuelta-2*h);
t=[tida tvuelta];
Vt=[Vida Vvuelta];

                 %step size sec.
%t=[0 0 taux];
% Vt=sr*t-2*h*sr;
Vt(1)=0;
Vt(2)=0;
Vc=zeros(1,length(t));
Vc(1)=0; 
Vc(2)=0;
If=zeros(1,length(t));
If(1)=0;                       % initial If ****** chequear
If(2)=0;
suma=0;                     % inicialmente nula

for i=3:(length(t)-1)    % calculation loop
    suma=0;
    for j=2:i-1
    suma=suma+(Vc(j)-Vc(j-1))*((t(i)-t(j))^(-a));
    end
         
    If(i)=Icor*(10^(Vc(i)/ba)-10^(-Vc(i)/bc));
    
    Rpi=ba*bc/(2.3*Icor*(bc*10^(Vc(i)/ba)+ba*10^(-Vc(i)/bc)));
   
    dVc=((Vt(i+1)-Vc(i))/Rs-If(i)-k*suma)/(k*(h^(-a))+1/Rs+1/Rpi); 

    Vc(i+1) = Vc(i)+dVc;  % main equation
end

Vcpe=Vc(3:length(t));
Vtot=Vt(3:length(t));
It=(Vtot-Vcpe)/Rs;

m=length(It);

parametros0005=zeros(m,8); parametros0005(1,:)=[sr Ec0005 Rs Icor ba bc Y0 a];
CP0005_O2_f3a=[Vtot' It' abs(It)' parametros0005]; dlmwrite('CP0005_O2_f3a.txt',CP0005_O2_f3a,'\t');

figure, semilogx(abs(It),Vtot,'b',abs(Itexp0005),Vtexp0005,'r')

