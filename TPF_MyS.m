clc
close all

global a b X Y Z k1 k2 k3 k4 k5 T t e e2 h q f;
%==============================================
%PARAMETROS Generales 
%==============================================
a=0.06; %concentracion de BrO3^-
b=0.03; %concentracion de MA+BrMA
X=0.0003 ;%concentracion de HBrO2
Y= 4;%Concentracion de Br^-
Z= 0.0002;%Concentracion de Me^(n+1)+ 
k1=2;
k2=3*(10^6);
k3=40;
k4=3000;
k5=0.6;
h=0.5;
T=1000;
%==============================================
%PARAMETROS DE LAS ECUACIONES 
%==============================================
e=(k5*b)/k3*a;
e2=2*k4/k2;
x=X*((2*k4)/(k3*a));
y=Y*((k2)/(k3*a));
z=Z*(((k3*a)^2)/(k4*k5*b));
q=(2*k1*k4)/(k2*k3);
f=2*h;

t=linspace(0,1/(k5*b),T);

%=======================================================================
%Condiciones iniciales
%=======================================================================
y0=[X,Z];% Condiciones iniciales
%R=[x,z]; %variables en forma de vector

%=========================================================================
% PUNTOS DE EQUILIBRIO
%=========================================================================
pe = fsolve(funb,y0); % funa() = 0
display (pe); % Muestro los puntos de equilibrio


%========================================================================
%Analisis del sistema no lineal
%=========================================================================
%graficos temporales

[t_sc,y_sc]=ode45 (func(),t,y0);
figure (1)
subplot(2,1,1);
plot (t_sc,y_sc(:,1),'LineWidth',2);
title('x en funcion del tiempo');
hold on

subplot(2,1,2);
plot (t_sc,y_sc(:,2),'LineWidth',2);
title('z en funcion del tiempo');
hold on

% Dia gramas de fase
figure (2)

plot (y_sc(:,1),y_sc(:,2),'LineWidth',2);
hold on;
plot(1.5e-04,1.5e-04,'k*');
plot(X,Z,'r*');


title ('x vs z');

%=========================================================================
% LINEALIZACION DEL SISTEMA
%=========================================================================
syms x z
R = [x,z];
A = jacobian(funa(),[x,z]); % Jacobiano del sistema
A = (subs(A,[x,z],pe)); % Jacobiano evaluado en punto de equilibrio
dF_lin = formula(A*(R-pe')); % Sistema lineal
display (dF_lin);



%=========================================================================
% Funcion de transferencia
%=========================================================================
a11=348488687213804437736221687163170654224307095346403021641576166292278513/95093584416526732986134454391201571865134545632064563707473302839925669888 ;
a12=348488687213804437736221687163170654224307095346403021641576166292278513/95093584416526732986134454391201571865134545632064563707473302839925669888 ;
a21=0;
a22=0;

A=[a11 a12; a21 a22];
pe=[1e-03*0.15;1*1e-03*0.15];
funst=@(t,R) [a11*(R(1)-pe(1))+a12*(R(2)-pe(2));
             a21*(R(1)-pe(1))+a22*(R(2)-pe(2))];
         
[tm,Rs]=ode45 (funst,0:(1/0.0018):1000,[3;5]);

figure(3)
subplot(2,1,1);
plot (tm,Rs(:,1),'LineWidth',2);
title('x en funcion del tiempo');

hold on

subplot(2,1,2);
plot (tm,Rs(:,2),'LineWidth',2);
title('z en funcion del tiempo');


[num,den] = ss2tf (A,pe,[1 1],0);
H= tf(num,den) %Funcion de transferencia 

figure (4)
pzmap (H)
[polos,zeros] = pzmap (H);

%==========================================================================
%FUNCIONES
%==========================================================================

%Para linealizar
function dFa=funa()
global q f x z
syms x z
dFx=x*(1-x)-f*z*(x-q)/(x+q);
dFz=x-z;
dFa=[dFx;dFz];
end

%Para buscar los pe
function dF =funb()
global q f
dF = @(R) [R(1)*(1-R(1))-f*R(2)*(R(1)-q)/(R(1)-q); 
        R(1)-R(2)];
end
%Para graficar
function dFc =func()
global q f
dFc = @(t,R) [R(1)*(1-R(1))-f*R(2)*(R(1)-q)/(R(1)-q); 
        R(1)-R(2)];
end

