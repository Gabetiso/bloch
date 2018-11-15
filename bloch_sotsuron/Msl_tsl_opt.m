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
Bos = 80e-9;

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
%Newton's methods
%-------------------------------------------------------------------------------
f = @(t) -M(2)/T1r*exp(-t/T1r)+exp(-al*t)*( al*M(2)*cos(be*t)+be*M(2)*sin(be*t)+A*al/be*sin(be*t)-A*cos(be*t) );%Function
df = @(t) M(2)/T1r^2*exp(-t/T1r)+exp(-al*t)*( -al^2*M(2)*cos(be*t)-2*al*be*M(2)*sin(be*t)-A*al^2/be*sin(be*t)...
          +2*al*A*cos(be*t)+A*be*sin(be*t)+be^2*M(2)*cos(be*t) ); %Derivative of f

t = [1e-4,1e-3,1e-2,1e-1,1.5e-1,2e-1,1];
N_max = 100;

for i = 1:size(t,2)
  t_old = t(i);

  for n = 1:N_max
    t_new = t_old - f(t_old)/df(t_old);
    if abs(t_new - t_old) < 1e-7  %Convergence determination
      break
    end
    t_old = t_new;
  end

  if n==N_max
    disp('Don''t convergent')
  else
    disp(t_old)
  end
end
