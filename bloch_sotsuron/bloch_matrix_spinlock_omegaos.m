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
%fos = 100; %brain frequency   %Hz
Bsl = (fsl * 2 * pi)/gamma;
Bos = 160e-9;
tsl = 50e-3;

fos = linspace(0,200,1e3);
omega_os = 2*pi * fos;
scr = zeros( size(fos) );

%-------------------------------------------------------------------------------
%function
%-------------------------------------------------------------------------------
[M] = bloch_first( T1, T2, b_x0, b_y0, trf, M_inf, M_i );

for i = 1:size(fos,2)
  [M_sl] = bloch_second( T1rho, T2rho, Bsl, Bos, omega_os(i), tsl, M );
  [M_sl2] = bloch_first( T1, T2, -b_x0, b_y0, trf, M_inf, M_sl );
  scr(i) = M_sl2(3);
end

figure;
plot(fos,scr);
xlabel('\omega_{os}');
ylabel('SCR');
xlim([0,200]);
ylim([0,1]);
ax = gca;
ax.FontName = 'Times New Roman';
ax.FontSize = 16;
grid on;
grid minor;
