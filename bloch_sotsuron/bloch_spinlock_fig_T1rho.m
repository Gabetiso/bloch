clear;
close all;

addpath('./bloch_spinlock'); %path of the using function(bloch_first and bloch_second)

gamma = 2 * pi * 42.58e6;
FA = deg2rad(90); %flip angle   %rad

%-------------------------------------------------------------------------------
%parameter of bloch_first
%-------------------------------------------------------------------------------
T1 = 121.7e-3;
T2 = 154.7e-3;
trf = 1e-3;  %given parameter
b_x0 = FA/(gamma*trf);
b_y0 = 0;
M_inf = 1;
M_i = [0; 0; 1];

%-------------------------------------------------------------------------------
%parameter of bloch_second
%-------------------------------------------------------------------------------
T1rho = linspace(50e-3,200e-3,4);
T2rho = T2;
fsl = 100; %spin lock frequency   %Hz
fos = 100; %brain frequency   %Hz
Bsl = (fsl * 2 * pi)/gamma;
Bos = 600e-9;
omega_os = fos * 2 * pi;
tsl = 40e-3;
%-------------------------------------------------------------------------------
%function
%-------------------------------------------------------------------------------
[M] = bloch_first_fig( T1, T2, b_x0, b_y0, trf, M_inf, M_i );
M_sl=zeros(3,500,size(T1rho,2));
M_rev=zeros(3,100,size(T1rho,2));
for i = 1:size(T1rho,2)
  [M_sl(:,:,i)] = bloch_second_fig( T1rho(i), T2rho, Bsl, Bos, omega_os, tsl, M(:,end) );
  [M_rev(:,:,i)] = bloch_first_fig( T1, T2, -b_x0, b_y0, trf, M_inf, M_sl(:,end,i) );
end
%-------------------------------------------------------------------------------
%figure
%-------------------------------------------------------------------------------
%axis
a = -1:1e-2:1;
b = zeros(size(a));

figure;
for i = 1:size(T1rho,2)
  subplot(2,2,i);
  plot3(M(1,:),M(2,:),M(3,:),'b');
  hold on;
  plot3(M_sl(1,:,i),M_sl(2,:,i),M_sl(3,:,i),'r');
  plot3(M_rev(1,:,i),M_rev(2,:,i),M_rev(3,:,i),'b');
  plot3(a,b,b,'k',b,b,a,'k',b,a,b,'k'); %axis
  hold off;
  title(sprintf('T1rho=%d/ms',T1rho(i)*1e3));
  xlabel('M_x');
  ylabel('M_y');
  zlabel('M_z');
  xlim([-1,1]);
  ylim([-1,1]);
  zlim([-1,1]);
  ax = gca;
  ax.FontName = 'Times New Roman';
  ax.FontSize = 16;
  grid on;
  grid minor;
end
orient(gcf,'landscape');
print(gcf,'./Result/spinlockfigT1rho','-dpdf','-bestfit');
