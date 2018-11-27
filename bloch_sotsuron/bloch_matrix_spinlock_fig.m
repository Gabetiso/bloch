clear;
close all;

addpath('./bloch_spinlock'); %path of the using function(bloch_first and bloch_second)

gamma = 2 * pi * 42.58e6;
FA = deg2rad(90); %flip angle   %rad

%-------------------------------------------------------------------------------
%parameter of bloch_first
%-------------------------------------------------------------------------------
T1 = 884e-3;
T2 = 72e-3;
trf = 1e-3;  %given parameter
b_x0 = FA/(gamma*trf);
b_y0 = 0;
M_inf = 1;
M_i = [0; 0; 1];

%-------------------------------------------------------------------------------
%parameter of bloch_second
%-------------------------------------------------------------------------------
T1rho = 150e-3;
T2rho = 70e-3;
fsl = 100; %spin lock frequency   %Hz
fos = 100; %brain frequency   %Hz
Bsl = (fsl * 2 * pi)/gamma;
Bos = 80e-9;
omega_os = fos * 2 * pi;
tsl = 50e-3;

%-------------------------------------------------------------------------------
%function
%-------------------------------------------------------------------------------
[M] = bloch_first_fig( T1, T2, b_x0, b_y0, trf, M_inf, M_i );
[M_sl] = bloch_second_fig( T1rho, T2rho, Bsl, Bos, omega_os, tsl, M(:,end) );
[M_sl2] = bloch_first_fig( T1, T2, -b_x0, b_y0, trf, M_inf, M_sl(:,end) );

%-------------------------------------------------------------------------------
%figure
%-------------------------------------------------------------------------------
M_x = zeros(1,size(M,2));
M_y = zeros(1,size(M,2));
M_z = zeros(1,size(M,2));
M_sl_x = zeros(1,size(M_sl,2));
M_sl_y = zeros(1,size(M_sl,2));
M_sl_z = zeros(1,size(M_sl,2));
M_sl2_x = zeros(1,size(M_sl2,2));
M_sl2_y = zeros(1,size(M_sl2,2));
M_sl2_z = zeros(1,size(M_sl2,2));
for i = 1:size(M,2)
  M_x(i) = M(1,i);
  M_y(i) = M(2,i);
  M_z(i) = M(3,i);
end
for j = 1:size(M_sl,2)
  M_sl_x(j) = M_sl(1,j);
  M_sl_y(j) = M_sl(2,j);
  M_sl_z(j) = M_sl(3,j);
end
for k = 1:size(M_sl2,2)
  M_sl2_x(k) = M_sl2(1,k);
  M_sl2_y(k) = M_sl2(2,k);
  M_sl2_z(k) = M_sl2(3,k);
end

a = -1:1e-2:1;
b = zeros(size(a));

figure;
plot3(M_x,M_y,M_z,'b');
hold on;
plot3(M_sl_x,M_sl_y,M_sl_z,'r');
plot3(M_sl2_x,M_sl2_y,M_sl2_z,'b');
plot3(a,b,b,'k');
plot3(b,a,b,'k');
plot3(b,b,a,'k');
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
saveas(gcf,'./Result/spinlock','png');
