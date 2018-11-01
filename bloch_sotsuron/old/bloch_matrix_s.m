%Dynamics of magnetization in secondary rotating coordinate system (i.e. while spin locking)
clear;
close all;

gamma = 2 * pi * 42.58e6; %rad/(s*T)
T1rho = 0.2;    %s
T2rho = 0.2;    %s
Bos = 1e-7; %T
Bsl = 1e-6; %spinlock magnetic field     %T
b_os_x0 = Bos/2; %Bosのt=0でのx成分 %T
b_os_z0 = 0; %Bosのt=0でのz成分 %T
omega_sl = gamma * Bsl;%angular frequency in secondary rotating coordinate system

timesize = 50e-3; %spinlock time
h =1e-5; %stepのこと
j = 1;

M_sl_x = zeros(1,int8(timesize/h)); %ここ要相談
M_sl_y = zeros(1,int8(timesize/h));
M_sl_z = zeros(1,int8(timesize/h));

R1rho = 2*pi/T1rho; %rad
R2rho = 2*pi/T2rho; %rad

omega_sl_x = gamma * b_os_x0;
omega_sl_y = gamma * Bsl;
omega_sl_z = gamma * b_os_z0;
d_omega_sl = omega_sl - omega_sl_y;

[M_x, M_y, M_z, i]=bloch_first(300e-3,150e-3,1e-5,0,5.8e-4,0.5);
M_sl_i = [M_x(i-1); M_y(i-1); M_z(i-1)]; %initial value of magnetization

MA = [-R2rho, omega_sl_z, d_omega_sl;...
      -omega_sl_z, -R1rho, omega_sl_x;...
      -d_omega_sl, -omega_sl_x, -R2rho];

    for t = 0:h:timesize
      M = expm(MA * t) * M_sl_i;
      M_sl_x(j) = M(1,1);
      M_sl_y(j) = M(2,1);
      M_sl_z(j) = M(3,1);
      j = j + 1;
    end

    figure;
    plot3(M_sl_x,M_sl_y,M_sl_z);
    xlabel('M_x');
    ylabel('M_y');
    zlabel('M_z');
    xlim([-1,1]);
    ylim([-1,1]);
    zlim([-1,1]);
    ax = gca;
    ax.FontName = 'Times New Roman';
    ax.FontSize = 16;
    grid on;
    grid minor;
