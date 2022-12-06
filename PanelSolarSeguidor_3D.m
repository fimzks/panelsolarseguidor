%Limpiar todo
clear all; clc; clf;

%cargar variables
load('PanelSolarSeguidor.mat')

%inicializacion

    Simulation_Time=10;

   %---valores iniciales---%
    secs=0;     %segundos iniciales
    mn=0;       %minutos  iniciales
    hra=1;      %hora     inicial
    dia=1;      %dia      inicial
    mes=1;      %mes      inicial

    latitud=  -36.83078129680368;
    longitud= -73.0368870295651;
    timezone= -3;

    beta_i = pi/4;    % Ángulo de elevación [rad]


  %---parametros---%
   % panel total %
    M = 24.5;       % Masa [Kg]
    P_W = 0.996;    % Ancho de panel total [m]
    P_L = 2.18;     % Largo de panel total [m]
    dph = 0.35;     % Profundidad [m]
    Kd = 5;         % Constante de amortiguamiento [(N*m)/(rad/s)]

   %-paneles secundarios-%
    S_w=0.1;        % Ancho de panel secundario [m]
    S_l=0.2;        % Largo de panel secundario [m]

   % motor %  
    Kf = 1.8350;    % Constante de FEM trasera [V/(rad/s)]
    Kt = 1.8350;    % Constante de Torque [N*m/A]
    L = 0.1042;     % Inductancia [H]
    R = 3.05;       % Resistencia [Ohm]
    Kg= 20;         % Razón entre engranajes

  %---Voltaje---%
   G_retroceso = -20;

  %---Eficiencia---%
   selector_eficiencia=0;
   eficiencia_custom=10;

%iniciar simulacion
    sim('PanelSolarSeguidor_mk2_comp.slx',Simulation_Time);
    simout=sim("PanelSolarSeguidor_mk2_comp.slx");

%Definir angulos
 %alfa
    alfa = ans.logsout{1}.Values.Data;
    alfa_time = ans.logsout{1}.Values.Time;
    n = ans.logsout{1}.Values.Length;
 %beta
    beta=ans.logsout{2}.Values.Data;
    beta_time = ans.logsout{2}.Values.Time;
    nb = ans.logsout{2}.Values.Length;

%Vertices panel
    a=[P_W/2,-P_L/2,0];
    b=[P_W/2,P_L/2,0];
    d=[P_W/2,-P_L/2,dph];
    c=[P_W/2,P_L/2,dph];
    f=[-P_W/2,P_L/2,0];
    g=[-P_W/2,P_L/2,dph];
    h=[-P_W/2,-P_L/2,dph];
    e=[-P_W/2,-P_L/2,0];

    Ry= [cos(beta(1)) 0 sin(beta(1));
            0      1      0    ; 
        -sin(beta(1)) 0 cos(beta(1))];

%Obtener puntos rotados en y del panel
    a = Ry * [a(1);a(2);a(3)];
    b = Ry * [b(1);b(2);b(3)];
    c = Ry * [c(1);c(2);c(3)];
    d = Ry * [d(1);d(2);d(3)];
    e = Ry * [e(1);e(2);e(3)];
    f = Ry * [f(1);f(2);f(3)];
    g = Ry * [g(1);g(2);g(3)];
    h = Ry * [h(1);h(2);h(3)];

%Matriz de puntos del panel
    x = [a(1) b(1) c(1) d(1) a(1) e(1) f(1) g(1) h(1) e(1) h(1) d(1) h(1) g(1) c(1) b(1) f(1)];
    y = [a(2) b(2) c(2) d(2) a(2) e(2) f(2) g(2) h(2) e(2) h(2) d(2) h(2) g(2) c(2) b(2) f(2)];
    z = [a(3) b(3) c(3) d(3) a(3) e(3) f(3) g(3) h(3) e(3) h(3) d(3) h(3) g(3) c(3) b(3) f(3)];

%Vertices del cubo
    A=[-0.5,-0.5,-2];
    B=[0.5,-0.5,-2];
    C=[0.5,0.5,-2];
    D=[-0.5,0.5,-2];
    E=[-0.5,-0.5,-1];
    F=[0.5,-0.5,-1];
    G=[0.5,0.5,-1];
    H=[-0.5,0.5,-1];

%Matriz de puntos del cubo
    k = [A(1) B(1) C(1) D(1) A(1) E(1) F(1) G(1) H(1) E(1) H(1) D(1) H(1) G(1) C(1) B(1) F(1)];
    p = [A(2) B(2) C(2) D(2) A(2) E(2) F(2) G(2) H(2) E(2) H(2) D(2) H(2) G(2) C(2) B(2) F(2)];
    v = [A(3) B(3) C(3) D(3) A(3) E(3) F(3) G(3) H(3) E(3) H(3) D(3) H(3) G(3) C(3) B(3) F(3)];

%Dibujar cubo
    Cubo = line(k,p,v); %Dibujar cubo
    Cubo.Color = 'green'; %Color

%Cubo.LineStyle = '--'; %Tipo de linea
    Cubo.LineWidth = 3; %Grosor

%Vertices del palo
    A=[-0.25,-0.25,-1];
    B=[0.25,-0.25,-1];
    C=[0.25,0.25,-1];
    D=[-0.25,0.25,-1];
    E=[-0.25,-0.25,0];
    F=[0.25,-0.25,0];
    G=[0.25,0.25,0];
    H=[-0.25,0.25,0];

%Matriz de puntos del palo
    k = [A(1) B(1) C(1) D(1) A(1) E(1) F(1) G(1) H(1) E(1) H(1) D(1) H(1) G(1) C(1) B(1) F(1)];
    p = [A(2) B(2) C(2) D(2) A(2) E(2) F(2) G(2) H(2) E(2) H(2) D(2) H(2) G(2) C(2) B(2) F(2)];
    v = [A(3) B(3) C(3) D(3) A(3) E(3) F(3) G(3) H(3) E(3) H(3) D(3) H(3) G(3) C(3) B(3) F(3)];

%Dibujar palo
    Palo = line(k,p,v); %Dibujar palo
    Palo.Color = 'green'; %Color

%Palo.LineStyle = '--'; %Tipo de linea
    Palo.LineWidth = 3; %Grosor

%Formato plano 
    Xmin = -2; %Limite negativo x para el plano
    Xmax = 2; %Limite positivo x para el plano
    Ymin = -2; %Limite negativo y para el plano
    Ymax = 2; %Limite positivo y para el plano
    Zmin = -2; %Limite negativo y para el plano
    Zmax = 2; %Limite positivo y para el plano

    axis([Xmin,Xmax,Ymin,Ymax,Zmin,Zmax]); %Tamaño de plano carteciano
    axis square; %Cuadricula cuadrada
    grid on %Mostrar cuadricula
    view(3); %Vista isometrica
    xlabel('x[m]'); %Eje x
    ylabel('y[m]'); %Eje y
    zlabel('z[m]'); %Eje z

%------------------------Animacion----------------------%
for i=1:n
  %Definir Matrice de rotacion
    Rz = [cos(alfa(i)) -sin(alfa(i)) 0 ;
          sin(alfa(i))  cos(alfa(i)) 0 ; 
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
        pause(alfa_time(i+1)-alfa_time(i));
    end
    delete(Panel_rotado);
end

%Definir Matrice de rotacion
    Rz = [cos(alfa(n)) -sin(alfa(n)) 0 ;
          sin(alfa(n))  cos(alfa(n)) 0 ; 
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