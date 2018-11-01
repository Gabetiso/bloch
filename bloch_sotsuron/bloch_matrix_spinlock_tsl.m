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
Bos = 160e-9;
%tsl = 50e-3;

tsl = linspace(0,500e-3,1e3); %variable
scr = zeros( size(tsl) );

%-------------------------------------------------------------------------------
%function
%-------------------------------------------------------------------------------
[M] = bloch_first( T1, T2, b_x0, b_y0, trf, M_inf, M_i );

for i = 1:size(tsl,2)
  [M_sl] = bloch_second( T1rho, T2rho, Bsl, Bos, omega_os, tsl(i), M );
  [M_sl2] = bloch_first( T1, T2, -b_x0, b_y0, trf, M_inf, M_sl );
  scr(i) = M_sl2(3);
end

figure;
plot(tsl*1e3,scr);
xlabel('T_{sl}[ms]');
ylabel('SCR');
xlim([0,500]);
ylim([-1,1]);
ax = gca;
ax.FontName = 'Times New Roman';
ax.FontSize = 16;
saveas(gcf,'tsl','png');