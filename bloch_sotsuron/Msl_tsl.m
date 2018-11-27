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
omega_os = 2 * pi * fos;
Bsl = (fsl * 2 * pi)/gamma;
Bos = 80e-9;

tsl = linspace(0,400e-3,1e3); %variable
M_sl_y = zeros( size(tsl) );

%-------------------------------------------------------------------------------
%function
%-------------------------------------------------------------------------------
[M] = bloch_first( T1, T2, b_x0, b_y0, trf, M_inf, M_i );

for i = 1:size(tsl,2)
  [M_sl] = bloch_second( T1rho, T2rho, Bsl, Bos, omega_os, tsl(i), M );
  %[M_sl2] = bloch_first( T1, T2, -b_x0, b_y0, trf, M_inf, M_sl );
  M_sl_y(i) = M(2) - M_sl(2);
end

figure;
plot(tsl*1e3,M_sl_y);
xlabel('T_{sl}[ms]');
ylabel('Difference of y component between M_{sl} and M');
xlim([0,400]);
ylim([0,1.1]);
ax = gca;
ax.FontName = 'Times New Roman';
ax.FontSize = 16;
saveas(gcf,'./Result/Msl_tsl','png');
