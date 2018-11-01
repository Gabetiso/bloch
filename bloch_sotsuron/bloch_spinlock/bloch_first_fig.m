function [M] = bloch_first_fig(T1, T2, b_x0, b_y0, trf, M_inf, M_i)
  %Dynamics of magnetization in primary rotating coordinate system
  gamma = 2 * pi * 42.58e6; %rad/(s*T)
  %T1 = 300e-3; %T1 time    %s
  %T2 = 150e-3; %T2 time    %s
  B0 = 0.3; %Static magnetic field    %T
  %b_x0 = 1e-5; %x component of B1 in t=0    %T
  %b_y0 = 0; %y component of B1 in t=0    %T
  omega_0 = gamma * B0;%Angular frequency in primary rotation coordinate
  %M_inf = 0.5; %M_z in thermal equilibrium state
  %trf = 5.8e-4;
  R1 = 1/T1; %1/s
  R2 = 1/T2; %1/s
  t = 0:1e-5:trf;
  M = zeros(3,size(t,2));

  omega = gamma * [b_x0; b_y0; B0] - [0;0;omega_0];

  M_a = [0; 0; R1*M_inf];
  %M_i   %Initial value of magnetization

  MA = [-R2, omega(3), -omega(2);...
  - omega(3),      -R2, omega(1);...
  omega(2),-omega(1), -R1];
  for i = 1:size(t,2)
    M(:,i) = expm(MA * t(i)) * M_i + ( eye(3) - expm(MA * t(i)) ) * -inv(MA) * M_a;
  end
end
