function [M_sl] = bloch_second(T1rho, T2rho, Bsl, Bos, omega_os, tsl, M_sl_i )
  %Dynamics of magnetization in secondary rotating coordinate system (i.e. while spin locking)
  %-----------------------------------------------------------------------------
  %detail of parameter
  %-----------------------------------------------------------------------------
  %T1rho = 0.2;    %s
  %T2rho = 0.5;    %s
  %Bsl = 0.3e-3; %Spinlock magnetic field     %T
  %Bos %Oscillating magnetic field    %T
  %omega_os = gamma * Bsl; %Angular frequency in secondary rotating coordinate system   %rad
  %tsl = 50e-3; %Spinlock time    %s
  %M_sl_i %Initial value of magnetization   %T
  %-----------------------------------------------------------------------------
  gamma = 2 * pi * 42.58e6; %rad/(s*T)

  b_os_x0 = Bos/2; %x component of Bos in t=0   %T
  b_os_z0 = 0; %z component of Bos in t=0   %T

  R1rho = 1/T1rho;    %1/s
  R2rho = 1/T2rho;    %1/s

  omega_sl = gamma * [b_os_x0; Bsl; b_os_z0] - [0;omega_os;0];    %rad

  MA = [-R2rho, omega_sl(3), - omega_sl(2);...
        -omega_sl(3),      -R1rho,   omega_sl(1);...
        omega_sl(2),-omega_sl(1), - R2rho];

  Mc2 = [cos(omega_os*tsl), 0, sin(omega_os*tsl);...
        0, 1, 0;...
        -sin(omega_os*tsl), 0, cos(omega_os*tsl)]; %Matrix for basis conversion

  M_sl = expm(MA * tsl) * M_sl_i;
  M_sl = Mc2 * M_sl; %basis conversion(2 -> 1)
end
