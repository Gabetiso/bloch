clear;
close all;

%addpath('../bloch_spinlock'); %path of the using function

gamma = 2 * pi * 42.58e6;
FA = 90; %flip angle (degree)

%-------------------------------------------------------------------------------
%parameter of bloch_first
%-------------------------------------------------------------------------------
T1 = 884e-3;
T2 = 72e-3;
b_x0 = 1e-5;
b_y0 = 0;
alpha_time = ( FA/(gamma*b_x0) )/180*pi;
M_inf = 1;
M_x_i = 0;
M_y_i = 0;
M_z_i = 1;

%-------------------------------------------------------------------------------
%parameter of bloch_second
%-------------------------------------------------------------------------------
T1rho = 0.13;
T2rho = 0.5;
Bsl = 1e-6;
Bos = 1e-6;
sl_time = 50e-3;

%-------------------------------------------------------------------------------
%function
%-------------------------------------------------------------------------------
[M_x, M_y, M_z, i] = bloch_first( T1, T2, b_x0, b_y0, alpha_time, M_inf, M_x_i, M_y_i, M_z_i );
[M_sl_x, M_sl_y, M_sl_z, j] = bloch_second( T1rho, T2rho, Bsl, Bos, sl_time, T1, T2, b_x0, b_y0, alpha_time, M_inf, M_x_i, M_y_i, M_z_i );
[M_sl2_x, M_sl2_y, M_sl2_z, k] = bloch_first( T1, T2, -b_x0, b_y0, alpha_time, M_inf, M_sl_x(j-1), M_sl_y(j-1), M_sl_z(j-1) );

d_M_z = M_sl2_z(k-1)/M_z(1);
disp('The difference in M_z')
disp(d_M_z)

%-------------------------------------------------------------------------------
%figure
%-------------------------------------------------------------------------------
figure;
plot3(M_x,M_y,M_z);
hold on;
plot3(M_sl_x,M_sl_y,M_sl_z);
plot3(M_sl2_x,M_sl2_y,M_sl2_z);
hold off;
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
