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
%parameter of Spin Lock
%-------------------------------------------------------------------------------
T1r = 150e-3;
T2r = 70e-3;
fsl = 100; %spin lock frequency   %Hz
fos = 100; %brain frequency   %Hz
omega_os = 2 * pi * fos;
Bsl = (fsl * 2 * pi)/gamma;
Bos = 8e-9;
tsl = linspace(0,1000e-3,1e3);

%-------------------------------------------------------------------------------
%parameter of Newton's methods
%-------------------------------------------------------------------------------
[M] = bloch_first( T1, T2, b_x0, b_y0, trf, M_inf, M_i );
R1r = 1/T1r;
R2r = 1/T2r;
b_os_x0 = Bos/2;
b_os_z0 = 0;
omega_sl = gamma * [b_os_x0; Bsl; b_os_z0];
al = (R1r+R2r)/2;
be = sqrt( omega_sl(1)^2 + omega_sl(3)^2 - (R1r-R2r)^2/4 );
A = (-al+R2r)*M(2)-omega_sl(3)*M(1)+omega_sl(1)*M(3);

%-------------------------------------------------------------------------------
%function
%-------------------------------------------------------------------------------
f = zeros(size(tsl));
df = zeros(size(tsl));
for i = 1:size(tsl,2)
  f(i) = M(2)*exp(-tsl(i)/T1r)-exp(-al*tsl(i))*( M(2)*cos(be*tsl(i))+A/be*sin(be*tsl(i)) );
  df(i) = -M(2)/T1r*exp(-tsl(i)/T1r)+exp(-al*tsl(i))*( al*M(2)*cos(be*tsl(i))+be*M(2)*sin(be*tsl(i))+A*al/be*sin(be*tsl(i))-A*cos(be*tsl(i)) );
end

figure;
plot(tsl*1e3,f,tsl*1e3,df);
xlabel('T_{sl}[ms]');
ylabel('difference of magnetization');
xlim([0,1000]);
ylim([-0.1,0.1]);
ax = gca;
ax.FontName = 'Times New Roman';
ax.FontSize = 16;
saveas(gcf,'Msl_tsl_opt','png');
