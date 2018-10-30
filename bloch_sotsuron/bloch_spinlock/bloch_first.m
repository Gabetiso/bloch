%Dynamics of magnetization in primary rotating coordinate system
function [M_x, M_y, M_z, i] = bloch_first(T1, T2, b_x0, b_y0, alpha_time, M_inf, M_x_i, M_y_i, M_z_i)

gamma = 2 * pi * 42.58e6; %rad/(s*T)
%T1 = 300e-3; %縦緩和時間;    %s
%T2 = 150e-3; %横緩和時間;    %s
B0 = 0.3; %静磁場         %T
%b_x0 = 1e-5; %B1のt=0でのx成分 %T
%b_y0 = 0; %B1のt=0でのy成分 %T
omega_0 = gamma * B0;%一次回転座標の角周波数(B1の角周波数)

%M_hは斉時解(ベクトル)
%M_sは特解(ベクトル)
%M_inf = 0.5;%熱平衡状態でのM_zの値(スカラー)
%alpha_time = 5.8e-4;
h =1e-5; %stepのこと
i = 1;

M_x = zeros(1,int8(alpha_time/h)); %ここ要相談
M_y = zeros(1,int8(alpha_time/h));
M_z = zeros(1,int8(alpha_time/h));

R1 = 2*pi/T1; %rad
R2 = 2*pi/T2; %rad

omega_x = gamma * b_x0;
omega_y = gamma * b_y0;
omega_z = gamma * B0;
d_omega = omega_0 - omega_z;

M_a = [0; 0; R1*M_inf];
M_i = [M_x_i; M_y_i; M_z_i]; %磁化の初期値(ベクトル)

I = eye(3);

MA = [-R2, -d_omega, -omega_y;...
      d_omega, -R2, omega_x;...
      omega_y, -omega_x, -R1];

M_s = -inv(MA) * M_a;
    for t = 0:h:alpha_time
      M = expm(MA * t) * M_i + (I - expm(MA * t)) * M_s;
      M_x(i) = M(1,1);
      M_y(i) = M(2,1);
      M_z(i) = M(3,1);
      i = i + 1;
    end

end
