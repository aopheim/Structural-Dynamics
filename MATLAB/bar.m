%%%%%%%%%% function for building bar elements  %%%%%%%%%%

function [Kel,Mel] = bar(node1,node2,mat),

%  node1=[x1 y1 z1];
%  node2=[x2 y2 z2];
%  mat  =[E G m A Ix  Iy  Iz]

l = norm(node1-node2);
E = mat(1);
A = mat(4);
rho = mat(3);

Kel_local  = E*A/l*[ 1 -1;
                    -1  1];
                 
Mel_local  = rho*A*l/6*[2 0 0 1 0 0;
 						0 2 0 0 1 0;
                        0 0 2 0 0 1;
                        1 0 0 2 0 0;
                        0 1 0 0 2 0;
                        0 0 1 0 0 2];
                     
e_x = (node2-node1)/l;
                     
T   = [e_x 0 0 0; 0 0 0 e_x];
Kel = T' * Kel_local *T;    % (4x4) matrix
Mel = Mel_local ;       % (6x6) matrix

