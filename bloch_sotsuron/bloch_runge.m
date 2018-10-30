clear;
close all;

size = 1e3;
G = 42.6;%磁気回転比
B1 = 1;%x方向へのRF磁場

%回転座標系での磁化を定義
M_x = zeros(1,size);
M_x(1) = 0;
M_y = zeros(1,size);
M_y(1) = 0;
M_z = zeros(1,size);
M_z(1) = 1;
y1 = zeros(1,size);
y2 = zeros(1,size);
y3 = zeros(1,size);
y4 = zeros(1,size);
z1 = zeros(1,size);
z2 = zeros(1,size);
z3 = zeros(1,size);
z4 = zeros(1,size);

h = 0.001;
t = (0:h:size*h);

%まずはeular法
%for i = 1:size
%  M_x(i+1) = M_x(i);
%  M_y(i+1) = M_y(i) + h*G*B1*M_z(i);
%  M_z(i+1) = M_z(i) - h*G*B1*M_y(i);
%end

%runge-kutta法
for i = 1:size
  y1(i) = G*B1*M_z(i);
  z1(i) = -G*B1*M_y(i);
  y2(i) = G*B1*( M_z(i) + 0.5*h*z1(i) );
  z2(i) = -G*B1*( M_y(i) + 0.5*h*y1(i) );
  y3(i) = G*B1*( M_z(i) + 0.5*h*z2(i) );
  z3(i) = -G*B1*( M_y(i) + 0.5*h*y2(i) );
  y4(i) = G*B1*( M_z(i) + h*z3(i) );
  z4(i) = -G*B1*( M_y(i) + h*y3(i) );

  M_x(i+1) = M_x(i);
  M_y(i+1) = M_y(i) + h*(y1(i) + 2*y2(i) + 2*y3(i) + y4(i))/6;
  M_z(i+1) = M_z(i) + h*(z1(i) + 2*z2(i) + 2*z3(i) + z4(i))/6;
end


figure;

%plot(t,M_z);
plot3(M_x,M_y,M_z);
xlabel('M_x');
ylabel('M_y');
zlabel('M_z');
ax = gca;
ax.FontName = 'Times New Roman';
ax.FontSize = 16;
grid on;
grid minor;
