%Limpiar todo
clear all; clc; clf;

%cargar variables
load('PanelSolarSeguidor.mat')

%inicializacion
simout=sim("PanelSolarSeguidor.slx");
Simulation_Time=10
selector=0;


   %---parametros---%
    % panel %
    m = 24.5;       % Masa [Kg]
    w = 0.996;      % Ancho [m]
    l = 2.18;       % Largo [m]
    d = 0.035;      % Profundidad [m]
    A = w*l;        % Área [m^2]
    beta = pi/4;    % Ángulo de elevación [rad]
    Kd = 5;         % Constante de amortiguamiento [(N*m)/(rad/s)]
    beta2 = beta*180/pi;

    % motor %  
    Kf = 1.8350;    % Constante de FEM trasera [V/(rad/s)]
    Kt = 1.8350;    % Constante de Torque [N*m/A]
    L = 0.1042;     % Inductancia [H]
    R = 3.05;       % Resistencia [Ohm]
    Kg= 20;         % Razón entre engranajes

   %---entradas---%
    %rampa
    pendiente_rampa=10000;
    initial_pendiente=0;
    start_time_pendiente=0;
    
    %escalon
    steptime_step=10000;
    intial_step=0;
    final_step=40;
    sampletime_step=0;

    %onda
    amplitude_sinewave=1000;
    bias_sinewave=0;
    frec_sinewave=1;
    phase_sinewave=0;
    sampletime_sinewave=0;

    %constante
    const=10000;
Simulation_Time=10
simout=sim("PanelSolarSeguidor.slx");

%iniciar simulacion
sim('PanelSolarSeguidor.slx',Simulation_Time)


%Definir valores
Xmin = -2; %Limite negativo x para el plano
Xmax = 2; %Limite positivo x para el plano
Ymin = -2; %Limite negativo y para el plano
Ymax = 2; %Limite positivo y para el plano
Zmin = -2; %Limite negativo y para el plano
Zmax = 2; %Limite positivo y para el plano

%Definir valores animacion
frames = 10000;
t_pausa = 100;
t_pausaInicio = 0.1;

%Definir angulo
theta = simout.angulo.data;
theta_time = simout.angulo.time;
n = numel(theta)

%Vertices
a=[w/2,-l/2,0];
b=[w/2,l/2,0];
d=[w/2,-l/2,0.5];
c=[w/2,l/2,d];
f=[-w/2,l/2,0];
g=[-w/2,l/2,d];
h=[-w/2,-l/2,d];
e=[-w/2,-l/2,0];



Ry= [cosd(beta2) 0 sind(beta2);
        0      1      0    ; 
    -sind(beta2) 0 cosd(beta2)];

%Obtener puntos rotados en x
    a = Ry * [a(1);a(2);a(3)];
    b = Ry * [b(1);b(2);b(3)];
    c = Ry * [c(1);c(2);c(3)];
    d = Ry * [d(1);d(2);d(3)];
    e = Ry * [e(1);e(2);e(3)];
    f = Ry * [f(1);f(2);f(3)];
    g = Ry * [g(1);g(2);g(3)];
    h = Ry * [h(1);h(2);h(3)];

%Matriz de puntos
x = [a(1) b(1) c(1) d(1) a(1) e(1) f(1) g(1) h(1) e(1) h(1) d(1) h(1) g(1) c(1) b(1) f(1)];
y = [a(2) b(2) c(2) d(2) a(2) e(2) f(2) g(2) h(2) e(2) h(2) d(2) h(2) g(2) c(2) b(2) f(2)];
z = [a(3) b(3) c(3) d(3) a(3) e(3) f(3) g(3) h(3) e(3) h(3) d(3) h(3) g(3) c(3) b(3) f(3)];

%Dibujar Panel_original
Panel_original = line(x,y,z); %Dibujar panel
Panel_original.Color = 'blue'; %Color
%Panel_original.LineStyle = '--'; %Tipo de linea
Panel_original.LineWidth = 3; %Grosor

%Formato plano 
axis([Xmin,Xmax,Ymin,Ymax,Zmin,Zmax]); %Tamaño de plano carteciano
axis square; %Cuadricula cuadrada
grid on %Mostrar cuadricula
view(3); %Vista isometrica
xlabel('x'); %Eje x
ylabel('y'); %Eje y
zlabel('z'); %Eje z

%------------------------Animacion----------------------%

pause(t_pausaInicio) 
delete(Panel_original)

for i=1:n
    %Definir Matrice de rotacion
    Rz = [cosd(theta(i)) -sind(theta(i)) 0 ;
          sind(theta(i))  cosd(theta(i)) 0 ; 
                0                 0      1];
    
    %Obtener puntos rotados en z
    ar = Rz * [a(1);a(2);a(3)];
    br = Rz * [b(1);b(2);b(3)];
    cr = Rz * [c(1);c(2);c(3)];
    dr = Rz * [d(1);d(2);d(3)];
    er = Rz * [e(1);e(2);e(3)];
    fr = Rz * [f(1);f(2);f(3)];
    gr = Rz * [g(1);g(2);g(3)];
    hr = Rz * [h(1);h(2);h(3)];
    
    %Vertices panel rotado z
    xr = [ar(1) br(1) cr(1) dr(1) ar(1) er(1) fr(1) gr(1) hr(1) er(1) hr(1) dr(1) hr(1) gr(1) cr(1) br(1) fr(1)];
    yr = [ar(2) br(2) cr(2) dr(2) ar(2) er(2) fr(2) gr(2) hr(2) er(2) hr(2) dr(2) hr(2) gr(2) cr(2) br(2) fr(2)];
    zr = [ar(3) br(3) cr(3) dr(3) ar(3) er(3) fr(3) gr(3) hr(3) er(3) hr(3) dr(3) hr(3) gr(3) cr(3) br(3) fr(3)];
    
    %Dibujar panel rotado 
    Panel_rotado = line(xr,yr,zr); %Dibujar panel
    Panel_rotado.Color = 'blue'; %Color
    %Panel_rotado.LineStyle = '--'; %Tipo de linea
    Panel_rotado.LineWidth = 3; %Grosor

    
    if i~=n
        pause(theta_time(i+1)-theta_time(i))
    end
    delete(Panel_rotado);
end

%Definir Matrice de rotacion
Rz = [cosd(theta(n)) -sind(theta(n)) 0 ;
      sind(theta(n))  cosd(theta(n)) 0 ; 
            0                 0      1];

%Obtener puntos rotados en z
ar = Rz * [a(1);a(2);a(3)];
br = Rz * [b(1);b(2);b(3)];
cr = Rz * [c(1);c(2);c(3)];
dr = Rz * [d(1);d(2);d(3)];
er = Rz * [e(1);e(2);e(3)];
fr = Rz * [f(1);f(2);f(3)];
gr = Rz * [g(1);g(2);g(3)];
hr = Rz * [h(1);h(2);h(3)];

%Vertices panel rotado z
xr = [ar(1) br(1) cr(1) dr(1) ar(1) er(1) fr(1) gr(1) hr(1) er(1) hr(1) dr(1) hr(1) gr(1) cr(1) br(1) fr(1)];
yr = [ar(2) br(2) cr(2) dr(2) ar(2) er(2) fr(2) gr(2) hr(2) er(2) hr(2) dr(2) hr(2) gr(2) cr(2) br(2) fr(2)];
zr = [ar(3) br(3) cr(3) dr(3) ar(3) er(3) fr(3) gr(3) hr(3) er(3) hr(3) dr(3) hr(3) gr(3) cr(3) br(3) fr(3)];

%Dibujar panel rotado 
Panel_rotado = line(xr,yr,zr); %Dibujar panel
Panel_rotado.Color = 'blue'; %Color
%Panel_rotado.LineStyle = '--'; %Tipo de linea
Panel_rotado.LineWidth = 3; %Grosor
%{
ti=datetime
t=ti;
tf=ti+seconds(Simulation_Time)
i=1;

while t<tf
    %while between(ti,tf)
    while (t-ti)<seconds(theta_time(i))
        % (t-ti)>seconds(theta_time(i))
    
        %Definir Matrice de rotacion
        Rz = [cosd(theta(i)) -sind(theta(i)) 0 ;
              sind(theta(i))  cosd(theta(i)) 0 ; 
                        0                 0  1];
        patata="hola"

        %Obtener puntos rotados en z
        ar = Rz * [a(1);a(2);a(3)];
        br = Rz * [b(1);b(2);b(3)];
        cr = Rz * [c(1);c(2);c(3)];
        dr = Rz * [d(1);d(2);d(3)];
        er = Rz * [e(1);e(2);e(3)];
        fr = Rz * [f(1);f(2);f(3)];
        gr = Rz * [g(1);g(2);g(3)];
        hr = Rz * [h(1);h(2);h(3)];
    
        %Vertices panel rotado z
        xr = [ar(1) br(1) cr(1) dr(1) ar(1) er(1) fr(1) gr(1) hr(1) er(1) hr(1) dr(1) hr(1) gr(1) cr(1) br(1) fr(1)];
        yr = [ar(2) br(2) cr(2) dr(2) ar(2) er(2) fr(2) gr(2) hr(2) er(2) hr(2) dr(2) hr(2) gr(2) cr(2) br(2) fr(2)];
        zr = [ar(3) br(3) cr(3) dr(3) ar(3) er(3) fr(3) gr(3) hr(3) er(3) hr(3) dr(3) hr(3) gr(3) cr(3) br(3) fr(3)];
    
        %Dibujar panel rotado 
        Panel_rotado = line(xr,yr,zr); %Dibujar panel
        Panel_rotado.Color = 'blue'; %Color
        %Panel_rotado.LineStyle = '--'; %Tipo de linea
        Panel_rotado.LineWidth = 3; %Grosor

        
        delete(Panel_rotado);
        i=i+1
    t=datetime
    t
    end
    
    
end

%Definir Matrice de rotacion
Rz = [cosd(theta(n)) -sind(theta(n)) 0 ;
      sind(theta(n))  cosd(theta(n)) 0 ; 
            0                 0      1];

%Obtener puntos rotados en z
ar = Rz * [a(1);a(2);a(3)];
br = Rz * [b(1);b(2);b(3)];
cr = Rz * [c(1);c(2);c(3)];
dr = Rz * [d(1);d(2);d(3)];
er = Rz * [e(1);e(2);e(3)];
fr = Rz * [f(1);f(2);f(3)];
gr = Rz * [g(1);g(2);g(3)];
hr = Rz * [h(1);h(2);h(3)];

%Vertices panel rotado z
xr = [ar(1) br(1) cr(1) dr(1) ar(1) er(1) fr(1) gr(1) hr(1) er(1) hr(1) dr(1) hr(1) gr(1) cr(1) br(1) fr(1)];
yr = [ar(2) br(2) cr(2) dr(2) ar(2) er(2) fr(2) gr(2) hr(2) er(2) hr(2) dr(2) hr(2) gr(2) cr(2) br(2) fr(2)];
zr = [ar(3) br(3) cr(3) dr(3) ar(3) er(3) fr(3) gr(3) hr(3) er(3) hr(3) dr(3) hr(3) gr(3) cr(3) br(3) fr(3)];

%Dibujar panel rotado 
Panel_rotado = line(xr,yr,zr); %Dibujar panel
Panel_rotado.Color = 'blue'; %Color
%Panel_rotado.LineStyle = '--'; %Tipo de linea
Panel_rotado.LineWidth = 3; %Grosor
%}

%{
condensador=8;

for i=2:n/condensador
    if(i == 1)
       pause(t_pausaInicio) 
    end

    %Definir Matrice de rotacion
    Rz = [cosd(theta(i*condensador)) -sind(theta(i*condensador)) 0 ;
          sind(theta(i*condensador))  cosd(theta(i*condensador)) 0 ; 
                        0                 0  1];

    %Obtener puntos rotados en z
    ar = Rz * [a(1);a(2);a(3)];
    br = Rz * [b(1);b(2);b(3)];
    cr = Rz * [c(1);c(2);c(3)];
    dr = Rz * [d(1);d(2);d(3)];
    er = Rz * [e(1);e(2);e(3)];
    fr = Rz * [f(1);f(2);f(3)];
    gr = Rz * [g(1);g(2);g(3)];
    hr = Rz * [h(1);h(2);h(3)];

    %Vertices panel rotado z
    xr = [ar(1) br(1) cr(1) dr(1) ar(1) er(1) fr(1) gr(1) hr(1) er(1) hr(1) dr(1) hr(1) gr(1) cr(1) br(1) fr(1)];
    yr = [ar(2) br(2) cr(2) dr(2) ar(2) er(2) fr(2) gr(2) hr(2) er(2) hr(2) dr(2) hr(2) gr(2) cr(2) br(2) fr(2)];
    zr = [ar(3) br(3) cr(3) dr(3) ar(3) er(3) fr(3) gr(3) hr(3) er(3) hr(3) dr(3) hr(3) gr(3) cr(3) br(3) fr(3)];

    %Dibujar panel rotado 
    Panel_rotado = line(xr,yr,zr); %Dibujar panel
    Panel_rotado.Color = 'blue'; %Color
    %Panel_rotado.LineStyle = '--'; %Tipo de linea
    Panel_rotado.LineWidth = 3; %Grosor

    pause(t_pausa); %Pausa

    delete(Panel_rotado);
end

    %Definir Matrice de rotacion
    Rz = [cosd(theta(n)) -sind(theta(n)) 0 ;
          sind(theta(n))  cosd(theta(n)) 0 ; 
                        0                 0      1];

    %Obtener puntos rotados en z
    ar = Rz * [a(1);a(2);a(3)];
    br = Rz * [b(1);b(2);b(3)];
    cr = Rz * [c(1);c(2);c(3)];
    dr = Rz * [d(1);d(2);d(3)];
    er = Rz * [e(1);e(2);e(3)];
    fr = Rz * [f(1);f(2);f(3)];
    gr = Rz * [g(1);g(2);g(3)];
    hr = Rz * [h(1);h(2);h(3)];

    %Vertices panel rotado z
    xr = [ar(1) br(1) cr(1) dr(1) ar(1) er(1) fr(1) gr(1) hr(1) er(1) hr(1) dr(1) hr(1) gr(1) cr(1) br(1) fr(1)];
    yr = [ar(2) br(2) cr(2) dr(2) ar(2) er(2) fr(2) gr(2) hr(2) er(2) hr(2) dr(2) hr(2) gr(2) cr(2) br(2) fr(2)];
    zr = [ar(3) br(3) cr(3) dr(3) ar(3) er(3) fr(3) gr(3) hr(3) er(3) hr(3) dr(3) hr(3) gr(3) cr(3) br(3) fr(3)];

    %Dibujar panel rotado 
    Panel_rotado = line(xr,yr,zr); %Dibujar panel
    Panel_rotado.Color = 'blue'; %Color
    %Panel_rotado.LineStyle = '--'; %Tipo de linea
    Panel_rotado.LineWidth = 3; %Grosor

%}
%{
for i = 1:frames
    %Condicion para que se pause un poco al inicio
    if(i == 1)
       pause(t_pausaInicio) 
    end
    
    %Convercion para avanzar un poco cada ciclo
    div = (i/frames);
    thetan = theta*div;
    
    %Definir Matrice de rotacion
    Rz = [cosd(thetan) -sind(thetan) 0 ;
          sind(thetan)  cosd(thetan) 0 ; 
                0            0       1];
          
    %Obtener puntos rotados en z
    ar = Rz * [a(1);a(2);a(3)];
    br = Rz * [b(1);b(2);b(3)];
    cr = Rz * [c(1);c(2);c(3)];
    dr = Rz * [d(1);d(2);d(3)];
    er = Rz * [e(1);e(2);e(3)];
    fr = Rz * [f(1);f(2);f(3)];
    gr = Rz * [g(1);g(2);g(3)];
    hr = Rz * [h(1);h(2);h(3)];

    %Vertices panel rotado z
    xr = [ar(1) br(1) cr(1) dr(1) ar(1) er(1) fr(1) gr(1) hr(1) er(1) hr(1) dr(1) hr(1) gr(1) cr(1) br(1) fr(1)];
    yr = [ar(2) br(2) cr(2) dr(2) ar(2) er(2) fr(2) gr(2) hr(2) er(2) hr(2) dr(2) hr(2) gr(2) cr(2) br(2) fr(2)];
    zr = [ar(3) br(3) cr(3) dr(3) ar(3) er(3) fr(3) gr(3) hr(3) er(3) hr(3) dr(3) hr(3) gr(3) cr(3) br(3) fr(3)];

    %Dibujar panel rotado 
    Panel_rotado = line(xr,yr,zr); %Dibujar panel
    Panel_rotado.Color = 'blue'; %Color
    %Panel_rotado.LineStyle = '--'; %Tipo de linea
    Panel_rotado.LineWidth = 3; %Grosor
    
    pause(t_pausa); %Pausa
    
    %Condicion para borrar los paneles anteriores
    if(i < frames)
       delete(Panel_rotado);
    end

    end
%}
