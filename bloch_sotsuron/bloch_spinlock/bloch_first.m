function [M] = bloch_first(T1, T2, b_x0, b_y0, trf, M_inf, M_i)
  %Dynamics of magnetization in primary rotating coordinate system
  %-----------------------------------------------------------------------------
  %detail of parameter
  %-----------------------------------------------------------------------------
  %T1 = 300e-3; %T1 time    %s
  %T2 = 150e-3; %T2 time    %s
  %b_x0 = 1e-5; %x component of B1 in t=0    %T
  %b_y0 = 0; %y component of B1 in t=0    %T
  %trf  %RF pulse time    %s
  %M_inf = 0.5; %M_z in thermal equilibrium state   %T
  %M_i  %Initial value of magnetization   %T
  %-----------------------------------------------------------------------------

  gamma = 2 * pi * 42.58e6; %rad/(s*T)
  B0 = 0.3; %Static magnetic field    %T
  omega_0 = gamma * B0; %Angular frequency in primary rotation coordinate   %rad

  R1 = 1/T1; %1/s
  R2 = 1/T2; %1/s

  omega = gamma * [b_x0; b_y0; B0] - [0;0;omega_0];

  M_a = [0; 0; R1*M_inf];

  MA = [-R2, omega(3), -omega(2);...
        -omega(3), -R2, omega(1);...
        omega(2), -omega(1), -R1];

  M = expm(MA * trf) * M_i + (eye(3) - expm(MA * trf) )* -inv(MA) * M_a;
end
