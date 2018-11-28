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
T1r = 100e-3;
T2r = 100e-3;
fsl = 100; %spin lock frequency   %Hz
fos = 100; %brain frequency   %Hz
omega_os = 2 * pi * fos;
Bsl = (fsl * 2 * pi)/gamma;
Bos = linspace(20e-9,600e-9,1e2);

%-------------------------------------------------------------------------------
%parameter of Newton's methods
%-------------------------------------------------------------------------------
[M] = bloch_first( T1, T2, b_x0, b_y0, trf, M_inf, M_i );
R1r = 1/T1r;
R2r = 1/T2r;
al = (R1r+R2r)/2;
%N_max = 100;
tsl = zeros(1,size(T1r,2));
%t = [1e-3,1e-2,4e-2,6e-2,8e-2,1e-1,1.5e-1];
%ta = zeros(size(t));

for i = 1:size(Bos,2)
  b_os_x0 = Bos(i)/2;
  b_os_z0 = 0;
  omega_sl = gamma * [b_os_x0; Bsl; b_os_z0];
  be = sqrt( omega_sl(1)^2 + omega_sl(3)^2 - (R1r-R2r)^2/4 );
  A = (-al+R2r)*M(2)-omega_sl(3)*M(1)+omega_sl(1)*M(3);
  %-------------------------------------------------------------------------------
  %Newton's methods
  %-------------------------------------------------------------------------------
  df = @(t) -M(2)/T1r*exp(-t/T1r)+exp(-al*t)*( al*M(2)*cos(be*t)+be*M(2)*sin(be*t)+A*al/be*sin(be*t)-A*cos(be*t) );
  ddf = @(t) M(2)/T1r^2*exp(-t/T1r)+exp(-al*t)*( -al^2*M(2)*cos(be*t)-2*al*be*M(2)*sin(be*t)-A*al^2/be*sin(be*t)...
  +2*al*A*cos(be*t)+A*be*sin(be*t)+be^2*M(2)*cos(be*t) ); %Derivative of df

  %  for j = 1:size(t,2)
  %    t_old = t(j);
  %
  %    for n = 1:N_max
  %      t_new = t_old - df(t_old)/ddf(t_old);
  %      if abs(t_new - t_old) < 1e-7  %Convergence determination
  %        break
  %      end
  %      t_old = t_new;
  %    end
  %
  %    if n==N_max
  %      ta(j) = 0;  %Convergence
  %    else
  %      if t_old < 0
  %        ta(j) = 0;  %Inappropriate
  %      elseif df(t_old - 1e-4)>0 & df(t_old + 1e-4)<0
  %        ta(j) = t_old;  %Ask for maximal value
  %      else
  %        ta(j) = 0;  %Minimal value
  %      end
  %    end
  %  end

  [ta] = newton_opt(df,ddf);
  [tsl_opt] = secondmin(ta);
  tsl(i) = tsl_opt;
end

figure;
plot(Bos*1e9,tsl*1e3);
xlabel('B_{os}(nT)');
ylabel('Optimized T_{sl}(ms)');
xlim([0,600]);
ylim([0,200]);
ax = gca;
ax.FontName = 'Times New Roman';
ax.FontSize = 16;
saveas(gcf,'./Result/Bos_tsl','png');
