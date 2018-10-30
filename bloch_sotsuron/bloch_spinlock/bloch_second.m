%Dynamics of magnetization in secondary rotating coordinate system (i.e. while spin locking)
function [M_sl_x, M_sl_y, M_sl_z, j] = bloch_second(T1rho, T2rho, Bsl, Bos, sl_time, T1, T2, b_x0, b_y0, alpha_time, M_inf, M_x_i, M_y_i, M_z_i)

gamma = 2 * pi * 42.58e6; %rad/(s*T)
%T1rho = 0.2;    %s
%T2rho = 0.5;    %s
%Bsl = 0.3e-3; %spinlock magnetic field     %T
%Bos %oscillating magnetic field  %T
b_os_x0 = Bos/2; %x component of Bos in t=0 %T
b_os_z0 = 0; %z component of Bos in t=0 %T
omega_sl = gamma * Bsl;%angular frequency in secondary rotating coordinate system

%sl_time = 50e-3; %spinlock time
h =1e-5; %step
j = 1;

M_sl_x = zeros(1,int8(sl_time/h)); %ここ要相談
M_sl_y = zeros(1,int8(sl_time/h));
M_sl_z = zeros(1,int8(sl_time/h));

R1rho = 2*pi/T1rho; %rad
R2rho = 2*pi/T2rho; %rad

omega_sl_x = gamma * b_os_x0;
omega_sl_y = gamma * Bsl;
omega_sl_z = gamma * b_os_z0;
d_omega_sl = omega_sl - omega_sl_y;

[M_x, M_y, M_z, i]=bloch_first(T1, T2, b_x0, b_y0, alpha_time, M_inf, M_x_i, M_y_i, M_z_i);

M_sl_i = [M_x(i-1); M_y(i-1); M_z(i-1)]; %initial value of magnetization

MA = [-R2rho, omega_sl_z, d_omega_sl;...
      -omega_sl_z, -R1rho, omega_sl_x;...
      -d_omega_sl, -omega_sl_x, -R2rho];

    for t = 0:h:sl_time
      M = expm(MA * t) * M_sl_i;
      M_sl_x(j) = M(1,1);
      M_sl_y(j) = M(2,1);
      M_sl_z(j) = M(3,1);
      j = j + 1;
    end

end
