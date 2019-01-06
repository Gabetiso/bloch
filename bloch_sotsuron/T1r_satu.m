clear;
close all;

addpath('./bloch_spinlock'); %path of the using function(bloch_first and bloch_second)

gamma = 2 * pi * 42.58e6;
FA = deg2rad(90); %flip angle   %rad

%-------------------------------------------------------------------------------
%parameter of bloch_first
%-------------------------------------------------------------------------------
T1 = 121.7e-3;
T2 = 154.7e-3;
trf = 1e-3;  %given parameter
b_x0 = FA/(gamma*trf);
b_y0 = 0;
M_inf = 1;
M_i = [0; 0; 1];

%-------------------------------------------------------------------------------
%parameter of Spin Lock
%-------------------------------------------------------------------------------
T1r = linspace(10e-3,360e-3,6);
T2r = T2;
fsl = 100; %spin lock frequency   %Hz
fos = 100; %brain frequency   %Hz
omega_os = 2 * pi * fos;
Bsl = (fsl * 2 * pi)/gamma;
Bos = 600e-9;
tsl = linspace(0,1000e-3,1e3);

%-------------------------------------------------------------------------------
%parameter of Newton's methods
%-------------------------------------------------------------------------------
[M] = bloch_first( T1, T2, b_x0, b_y0, trf, M_inf, M_i );
R2r = 1/T2r;
b_os_x0 = Bos/2;
b_os_z0 = 0;
omega_sl = gamma * [b_os_x0; Bsl; b_os_z0];

%-------------------------------------------------------------------------------
%function
%-------------------------------------------------------------------------------
f = zeros(size(T1r,2),size(tsl,2));
for i = 1:size(T1r,2)
  R1r = 1/T1r(i);
  al = (R1r+R2r)/2;
  be = sqrt( omega_sl(1)^2 + omega_sl(3)^2 - (R1r-R2r)^2/4 );
  A = (-al+R2r)*M(2)-omega_sl(3)*M(1)+omega_sl(1)*M(3);
  for j = 1:size(tsl,2)
    f(i,j) = M(2)*exp(-tsl(j)/T1r(i))-exp(-al*tsl(j))*( M(2)*cos(be*tsl(j))+A/be*sin(be*tsl(j)) );
  end
end

figure;
plot(tsl*1e3,f);
legend('T_{1\rho}=10ms','T_{1\rho}=80ms','T_{1\rho}=150ms','T_{1\rho}=220ms','T_{1\rho}=290ms','T_{1\rho}=360ms');
xlabel('T_{sl}[ms]');
xlim([0,100]);
ax = gca;
ax.FontName = 'Times New Roman';
ax.FontSize = 16;
saveas(gcf,'./Result/T1r_satu','png');
