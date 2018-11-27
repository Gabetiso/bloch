clear;
close all;

addpath('./bloch_spinlock'); %path of the using function(bloch_first and bloch_second)

gamma = 2 * pi * 42.58e6;
FA = deg2rad(90); %flip angle   %rad

%-------------------------------------------------------------------------------
%parameter of bloch_first
%-------------------------------------------------------------------------------
T1 = 68e-3;
T2 = 81e-3;
trf = 1e-3;  %given parameter
b_x0 = FA/(gamma*trf);
b_y0 = 0;
M_inf = 1;
M_i = [0; 0; 1];

%-------------------------------------------------------------------------------
%parameter of bloch_second
%-------------------------------------------------------------------------------
T1r = 85e-3;
T2r = 60e-3;
fsl = 100; %spin lock frequency   %Hz
fos = 100; %brain frequency   %Hz
omega_os = 2 * pi * fos;
Bsl = (fsl * 2 * pi)/gamma;
Bos = 80e-9;
%tsl = 50e-3;

tsl = linspace(0,500e-3,200); %variable %1/fsl=10msの周期にするべき
scr = zeros( size(tsl) );
scr_t1r = zeros( size(tsl) );

%-------------------------------------------------------------------------------
%function
%-------------------------------------------------------------------------------
[M] = bloch_first( T1, T2, b_x0, b_y0, trf, M_inf, M_i );
[M_sl2off] = bloch_first( T1, T2, -b_x0, b_y0, trf, M_inf, M );

for i = 1:size(tsl,2)
  %SIRS
  [M_sl] = bloch_second( T1r, T2r, Bsl, Bos, omega_os, tsl(i), M );
  [M_sl2on] = bloch_first( T1, T2, -b_x0, b_y0, trf, M_inf, M_sl );

  %T1rho relaxation
  M_t1r = [0;M(2)*exp(-tsl(i)/T1r);0];
  [M_t1r2] = bloch_first( T1, T2, -b_x0, b_y0, trf, M_inf, M_t1r );

  %signal change ratio
  scr(i) = M_sl2on(3)/M_sl2off(3);
  scr_t1r(i) = M_t1r2(3)/M_sl2off(3);
end

figure;
plot(tsl*1e3,abs(scr));
hold on;
plot(tsl*1e3,abs(scr_t1r));
plot(tsl*1e3,abs(scr_t1r - scr));
hold off;
legend('SIRS','T_{1\rho} relaxation','T_{1\rho} relaxation - SIRS');
xlabel('T_{sl}[ms]');
ylabel('M_{on} / M_{ref}');
xlim([0,500]);
ylim([0,1]);
ax = gca;
ax.FontName = 'Times New Roman';
ax.FontSize = 16;
saveas(gcf,'./Result/tsl','png');
