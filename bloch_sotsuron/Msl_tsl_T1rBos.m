clear;
close all;

addpath('./bloch_spinlock'); %path of the using function(bloch_first and bloch_second)
addpath('./Newton');

gamma = 2 * pi * 42.58e6;
FA = deg2rad(90); %flip angle   %rad

%-------------------------------------------------------------------------------
%parameter of bloch_first
%-------------------------------------------------------------------------------
T1 = 100e-3;
T2 = 80e-3;
trf = 1e-3;  %given parameter
b_x0 = FA/(gamma*trf);
b_y0 = 0;
M_inf = 1;
M_i = [0; 0; 1];

%-------------------------------------------------------------------------------
%parameter of Spin Lock
%-------------------------------------------------------------------------------
T1r = linspace(1e-3,500e-3,1e2);
T2r = 100e-3;
fsl = 100; %spin lock frequency   %Hz
fos = 100; %brain frequency   %Hz
omega_os = 2 * pi * fos;
Bsl = (fsl * 2 * pi)/gamma;
Bos = linspace(80e-9,560e-9,7);

%-------------------------------------------------------------------------------
%parameter of Newton's methods
%-------------------------------------------------------------------------------
[M] = bloch_first( T1, T2, b_x0, b_y0, trf, M_inf, M_i );
R2r = 1/T2r;
tsl = zeros(size(Bos,2),size(T1r,2));

for i = 1:size(Bos,2)
  b_os_x0 = Bos(i)/2;
  b_os_z0 = 0;
  omega_sl = gamma * [b_os_x0; Bsl; b_os_z0];
  for j = 1:size(T1r,2)
    R1r = 1/T1r(j);
    al = (R1r+R2r)/2;
    be = sqrt( omega_sl(1)^2 + omega_sl(3)^2 - (R1r-R2r)^2/4 );
    A = (-al+R2r)*M(2)-omega_sl(3)*M(1)+omega_sl(1)*M(3);
    %-------------------------------------------------------------------------------
    %Newton's methods
    %-------------------------------------------------------------------------------
    df = @(t) -M(2)/T1r(j)*exp(-t/T1r(j))+exp(-al*t)*( al*M(2)*cos(be*t)+be*M(2)*sin(be*t)+A*al/be*sin(be*t)-A*cos(be*t) );
    ddf = @(t) M(2)/T1r(j)^2*exp(-t/T1r(j))+exp(-al*t)*( -al^2*M(2)*cos(be*t)-2*al*be*M(2)*sin(be*t)-A*al^2/be*sin(be*t)...
    +2*al*A*cos(be*t)+A*be*sin(be*t)+be^2*M(2)*cos(be*t) ); %Derivative of df

    [ta] = newton_opt(df,ddf);
    [tsl_opt] = secondmin(ta);
    tsl(i,j) = tsl_opt;
  end
end

figure;
plot(T1r*1e3,tsl*1e3);
legend('B_{os}=80nT','B_{os}=160nT','B_{os}=240nT','B_{os}=320nT','B_{os}=400nT','B_{os}=480nT','B_{os}=560nT','Location','northwest');
xlabel('T_{1\rho}(ms)');
ylabel('Optimized T_{sl}(ms)');
xlim([0,200]);
ylim([0,50]);
ax = gca;
ax.FontName = 'Times New Roman';
ax.FontSize = 16;
saveas(gcf,'./Result/T1rbos_tsl','png');
