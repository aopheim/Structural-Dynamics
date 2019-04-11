%%%%%%%%%% function for building beam elements  %%%%%%%%%%

function [Kel,Mel] = beam(node1,node2,node3,mat),

%  node1=[x1 y1 z1];
%  node2=[x2 y2 z2];
%  node3=[x3 y3 z3]; for Ix direction
%  mat  =[E G m A Ix  Iy  Iz]

l = norm(node1-node2);          %Absoluttverdi
if l==0,
    disp('zero length')
    node1
    node2
    return
end
E = mat(1);
G = mat(2);
rho = mat(3);
A = mat(4);
J_x = mat(5);
I_y = mat(6);
I_z = mat(7);
r2  = (I_y+I_z)/A;

Kel_local  =  zeros(12,12);             %six degrees of freedom * 2. Matrix is in page 35 in lecture notes

Kel_local([1 7],[1 7]) = E*A/l*[ 1 -1;
                                -1  1];
Kel_local([2 6 8 12],[2 6 8 12]) = E*I_z/l^3* ...
  [12   6*l   -12   6*l;
   6*l  4*l^2 -6*l  2*l^2;
   -12  -6*l  12    -6*l;
   6*l  2*l^2 -6*l  4*l^2];

Kel_local([3 5 9 11],[3 5 9 11]) = E*I_y/l^3* ...
  [12   -6*l   -12   -6*l;
   -6*l  4*l^2 6*l  2*l^2;
   -12   6*l   12    6*l;
   -6*l  2*l^2 6*l  4*l^2];

Kel_local([4 10],[4 10]) = G*J_x/l*[ 1 -1;
                                    -1  1];
                                 
Mel_local  =  zeros(12,12);

Mel_local([1 7],[1 7]) = rho*A*l/6*[ 2 1;
                                     1 2];
Mel_local([2 6 8 12],[2 6 8 12]) = rho*A*l/420* ...
   [156   22*l   54   -13*l;
    22*l  4*l^2  13*l -3*l^2;
    54    13*l   156  -22*l;
   -13*l  -3*l^2 -22*l 4*l^2];

Mel_local([3 5 9 11],[3 5 9 11]) = rho*A*l/420* ...
   [156   -22*l   54   13*l;
    -22*l  4*l^2  -13*l -3*l^2;
    54    -13*l   156  22*l;
    13*l  -3*l^2  22*l 4*l^2];
Mel_local([4 10],[4 10]) = rho*A*l*r2/6*[ 2 1;
                                         1 2];
                                 

e_x = (node2-node1)/l;
d2  = node2-node1;
d3  = node3-node1;
e_y = cross(d3,d2)/norm(cross(d3,d2));
e_z = cross(e_x,e_y);           %cross product

R = [ e_x(1) e_x(2) e_x(3);
      e_y(1) e_y(2) e_y(3);
      e_z(1) e_z(2) e_z(3)];
   
T = zeros(12,12);
T([1:3],[1:3]) = R;
T([4:6],[4:6]) = R;
T([7:9],[7:9]) = R;
T([10:12],[10:12]) = R;

Kel = T'*Kel_local*T;
Mel = T'*Mel_local*T;
