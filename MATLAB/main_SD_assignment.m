%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Finite Element assembly of a beam-bar 3D structure
%
%    Project part of the course Structural Dynamics, TUM MW 2136
%    D.R. 05.06.2017
%
%   This file contains the assembly of a general 3D mesh consisting of beam
%   and bar elements. The nodal, material and connectivity matrices are
%   defined in the 'Nodes', 'Mat' and 'Elements' matrices repsectively. The
%   example given here is a hanguar as described in a joined document.
%   The beam element matrices are given. The building of the K and M matrices, and
%   definining the geometry of the structure is done by Daniel Rixen.
%   Implementation of the inverse iteration and Newmark integration done by
%   Adrian Opheim, July 2017.
%
%%%%%%%%%%%%%%%%%%%%%

 clear
 close all

%%%%%%%%%% Geometric parameters %%%%%%%%%%%%%
%
% define here all the geometric quantities you need (diameters, lengths, thinkness, etc.)

L=80;
l1=4;
l2=5;
l3=5;
l4=4;
l5=6;
l6=4;
l7=10;
l8=3;
l9=9;
l10=10;

h=25;

H=25;

D1=0.20; d1=0.195;    % beams group 1 (lower arch)
D2=0.20; d2=0.195;    % beams group 2 (mid arch)
D3=0.25; d3=0.244;    % beams group 3 (high arch and span)
D4=0.05; d4=0.046;    % bars in arch
D5=0.40; d5=0.380;    % column beams
D6=0.08; d6=0.075;    % cross bars
D7=0.25; d7=0.24;    % cross beams


%  Material properties 
%
% define here the material constant. 

 E=2.1e11;
 nu=0.3;
 rho=7500;

%%%%%%% Computing characteristics %%%%%%%%%%%%%
%
% Here you compute the quantities needed for the 'Mat' array (see below). 
%

 G = E/(2*(1+nu));

 A1  = pi*(D1^2-d1^2)/4;
 A2  = pi*(D2^2-d2^2)/4;
 A3  = pi*(D3^2-d3^2)/4;
 A4  = pi*(D4^2-d4^2)/4;
 A5  = pi*(D5^2-d5^2)/4;
 A6  = pi*(D6^2-d6^2)/4;
 A7  = pi*(D7^2-d7^2)/4;
 
 I1  = pi*(D1^4-d1^4)/64; J1  = pi*(D1^4-d1^4)/32;
 I2  = pi*(D2^4-d2^4)/64; J2  = pi*(D2^4-d2^4)/32;
 I3  = pi*(D3^4-d3^4)/64; J3  = pi*(D3^4-d3^4)/32;
 I5  = pi*(D5^4-d5^4)/64; J5  = pi*(D5^4-d5^4)/32;
 I7  = pi*(D7^4-d7^4)/64; J7  = pi*(D7^4-d7^4)/32;

%%%%%%%  Nodes    %%%%%%%%
%  x    ,   y   ,   z  
%
% Create here the 'Nodes' matrix. It contains the nodal coordinates x,y,z.
%

Nodes(1,:)  = [ 0               0       H ];    % arch
Nodes(2,:)  = [ 0               0       H+l1];  
Nodes(3,:)  = [ l1              0       H+l1];
Nodes(4,:)  = [ l3              0       H+l1+l2];      
Nodes(5,:)  = [ l3+l1           0       H+l1+l2];
Nodes(6,:)  = [ l3+l5           0       H+l1+l2+l4];
Nodes(7,:)  = [ l3+l5+l1        0       H+l1+l2+l4];
Nodes(8,:)  = [ l3+l5+l7        0       H+l1+l2+l4+l6];
Nodes(9,:)  = [ l3+l5+l7+l9     0       H+l1+l2+l4+l6+l8];   %This is the node where the force is applied. 
Nodes(10,:) = [ l3+l5+l7+l9+l10 0       H+l1+l2+l4+l6+l8]; % This is the middle node of the arch
for i=21:29,
    Nodes(i,:) = Nodes(i-20,:); 
    Nodes(i,1) = L-Nodes(i,1);
end
Nodes(30,:)   = [0     0    0]; % ground nodes
Nodes(31,:)   = [L     0    0];
for i=1:31,
   Nodes(i+100,:) = Nodes(i,:) + [0 h 0];
end
Nodes(135,:) = [0 h/2 H];
Nodes(136,:) = [L h/2 H];

%%%%%%% Dummy node for beams %%%%%%%%%
%
% this node is use to orient the bam cross section in 3D (see lecture
% notes).
% Since in this example the cross section is axisymmetric any orientation
% is fine. This we choose a dummy node randonly, just making sure it is not
% colinear with the 2 nodes of a beam.

 Nodes(200,:) = [0 40 0]; 

%%%%%%%  Material and section properties    %%%%%
%
%  E    ,   G   , rho ,   A   ,  Ix  ,   Iy   , Iz
%  
% Matrix 'Mat' contains the material and section properties definition. Here below an example with
% 7 different materials. Materials 1, 2, 3, 5 and 7 are for beam elements, materials 4 and 6 are for bar elements.

 Mat(1,1:7)=  [ E   G   rho   A1  J1    I1     I1];  % beams group 1 (lower arch)
 Mat(2,:)  =  [ E   G   rho   A2  J2    I2     I2];  % beams group 2 (mid arch)
 Mat(3,:)  =  [ E   G   rho   A3  J3    I3     I3];  % beams group 3 (high arch and span)
 Mat(4,:)  =  [ E   G   rho   A4  0      0      0];  % bars in arch
 Mat(5,:)  =  [ E   G   rho   A5  J5    I5     I5];  % column beams
 Mat(6,:)  =  [ E   G   rho   A6  0      0      0];  % cross bars
 Mat(7,:)  =  [ E   G   rho   A7  J7    I7     I7];  % cross beams

%%%%%%% Elements   %%%%%%%%%
%
%  type(1=bar,2=beam), mat, node 1 , node 2, node 3(for beam)
%
%  The Matrix 'Elements' contains the connectivity (localization) informations for all the elements of the mesh,
%  together with the material and geometric related properties. 
%  For example, here below: element no.1 is a beam (type=2), is made of
%  material 7, connects node 1 and node 135 and has the dummy node 200
%  (necessary for a beam)

Elements(1,:)  = [2   1       1  2  200];  % arch
Elements(2,:)  = [2   1       1  3  200];
Elements(3,:)  = [2   1       2  3  200];
Elements(4,:)  = [2   2       2  4  200];  
Elements(5,:)  = [2   2       3  5  200]; 
Elements(6,:)  = [2   2       4  6  200]; 
Elements(7,:)  = [2   2       5  7  200]; 
Elements(8,:)  = [2   2       6  8  200];  
Elements(9,:)  = [2   2       7  8  200];
Elements(10,:) = [2   3       8  9  200];
Elements(11,:) = [1   4       3  4  200];
Elements(12,:) = [1   4       4  5  200];
Elements(13,:) = [1   4       5  6  200];
Elements(14,:) = [1   4       6  7  200];
Elements(15,:) = [2   3       9  10 200];   % middle span
Elements(16,:) = [2   3       10 29 200];  
% create right side of arch
for i=1:14,
    Elements(16+i,:) = Elements(i,:);
    Elements(16+i,3:4) = Elements(16+i,3:4) + [20 20];
end
% vertical columns
Elements(31,:) = [2   5       31 21  200];
Elements(32,:) = [2   5       30  1  200];
% duplicate the arch at y=h
for i=1:32,
    Elements(i+32,:) = Elements(i,:)+[0 0 100 100 0];
end
% add beams for crane rails
Elements(65,:) = [2   7       1   135 200];
Elements(66,:) = [2   7       135 101 200];
Elements(67,:) = [2   7       21  136 200];
Elements(68,:) = [2   7       136 121 200];
% add cross bar and beams
Elements(69,:) = [1   6       9  129  200];
Elements(70,:) = [1   6       29 109  200];
Elements(71,:) = [1   6       8  108  200];
Elements(72,:) = [1   6       28 128  200];
Elements(73,:) = [1   6       6  104  200];
Elements(74,:) = [1   6       26 124  200];
Elements(75,:) = [1   6       4  106  200];
Elements(76,:) = [1   6       24 126  200];

% plot the mesh

figure(1),plotmesh(Nodes,Elements)  

%%%%%%%%% build elementary matrices and assemble  %%%%%%

% build index table for dof: locnod(i,:) = list of degrees of freedom number attached to node i
%
Nodes_active = unique(Elements(:,3:4));   % find which node really are used in elements. Finds the unique elements. 
                                          % the degrees of freedom fixed on the ground will
                                          % be removed later.
Ndof = size(Nodes_active,1)*6;            % 6 dof per node
locnod(Nodes_active,1:6) = reshape([1:Ndof],6,Ndof/6)';        % ':  transpose of matrix.
Nele   = size(Elements,1);

K = sparse(Ndof,Ndof);
M = sparse(Ndof,Ndof);

for el = 1:Nele,

   type   = Elements(el,1);
   matel  = Mat(Elements(el,2),:);      %materials for the specific element
   node1  = Nodes(Elements(el,3),:);  % coordinates of node 1
   node2  = Nodes(Elements(el,4),:);  % coordinates of node 2
   node3  = Nodes(Elements(el,5),:);  % coordinates of node 3 (for beams)

   % Build element matrices
   if type==1,      % bar
       % write the function yourself!!
      [Kel,Mel] = bar(node1,node2,matel);   %Kel = (4x4), Mel = (6x6)
         dof=[locnod(Elements(el,3),[1 2 3]) ...
              locnod(Elements(el,4),[1 2 3]) ] ; % only the translational
   elseif type==2,
      [Kel,Mel] = beam(node1,node2,node3,matel);  
         dof=[locnod(Elements(el,3),:) ...              % ... means line break
              locnod(Elements(el,4),:) ] ; 
   end
   
   %% Assemble
   K(dof,dof) = K(dof,dof)+Kel;
   M(dof,dof) = M(dof,dof)+Mel;
   
   
      
   clear matel type Kel Mel dof
end


%% PART 1 of assignment. 
% Checking the model: Verify that, if no degrees of freedom are fixed,
% K * u_trans = 0
% where u_trans is a translation rigid body mode with unit amplitude.

% K: (264x264) matrix
% M: (264x264) matrix
% u_trans: (264x1) matrix

u_trans = zeros(size(K,1),1);
for i = 3:6:size(u_trans)
    u_trans(i) = 1;         % Setting u1 = u2 = 1 (one unit distance). 
end

null_check = full(K) * u_trans;

disp('Max. value of null_check vector: ')
disp(max(null_check))       %Gives 7.4506e-09, which has to be within reasonable limits. 



%% Checking the total mass of the structure by computing
% u_trans' * M * u_trans = m_total

m_total = u_trans' * M * u_trans;

%% Calculating the actual mass of the structure: m = rho * A * L
% We have 7 different bars and beams, all with different A and L. 
% Following the naming convention given in drawing: 

beam1 = 8 * Mat(1,3) * Mat(1,4) * l1;       %8 beam 1 elements with length l1. Mat(1,3): rho, Mat(1,4): A
beam1 = beam1 + 4 * Mat(1,3) * Mat(1,4) * sqrt((2 * l1)^2);   % 4 beam1 elements at 45 degree angle.

beam2 = 8 * Mat(2,3) * Mat(2,4) * sqrt(l2^2 + l3^2);        % 8 lower beam2 elements
beam2 = beam2 + 8 * Mat(2,3) * Mat(2,4) * sqrt(l4^2 + l5^2);    % 8 middle beam2 elements
beam2 = beam2 + 4 * Mat(2,3) * Mat(2,4) * sqrt(l6^2 + l7^2);    % 4 top-top beam2 elements
beam2 = beam2 + 4 * Mat(2,3) * Mat(2,4) * sqrt((l7 - l1)^2 + l6^2);     % 4 top-lower beam2 elements

beam3 = 4 * Mat(3,3) * Mat(3,4) * l10;          %4 topmost beam3 elements
beam3 = beam3 + 4 * Mat(3,3) * Mat(3,4) * sqrt(l8^2 + l9^2);    % 4 beam3 elements with an angle


bar4 = 8 * Mat(4,3) * Mat(4,4) * l1;            % 8 horizontal bar4 elements with length l1
bar4 = bar4 + 4 * Mat(4,3) * Mat(4,4) * sqrt((l3 - l1)^2 + l2^2);       % 4 lowermost bar4 elements
bar4 = bar4 + 4 * Mat(4,3) * Mat(4,4) * sqrt((l5 - l1)^2 + l4^2);       % 4 uppermost bar4 elements

beam5 = 4 * Mat(5,3) * Mat(5,4) * H;        % 4 beam5 elements (foundation blocks)

bar6 = 2 * Mat(6,3) * Mat(6,4) * h;         % 2 vertical bar6 elements
bar6 = bar6 + 4 * Mat(6, 3) * Mat(6,4) * sqrt((l5^2 + h^2));            %4 outermost cross bar6 elements
bar6 = bar6 + 2 * Mat(6, 3) * Mat(6,4) * sqrt((2 * l10)^2 + h^2);

beam7 = 2 * Mat(7,3) * Mat(7,4) * h;         % 2 vertical beam7 elements

m_total_calc = beam1 + beam2 + beam3 + bar4 + beam5 + bar6 + beam7;

disp('Absolute value of difference between u_trans^T * M * u_trans and calculated total mass: ')
disp(abs(m_total_calc - m_total))

mass_error_percent = ((abs(m_total_calc - m_total)) / m_total) * 100;           %Currently at 0.65%. Seems very reasonable. 




%% Apply boundary conditions
% fix dofs of nodes 30 31 130 131
%a. find the degrees of freedom that are fixed
dof_fix = [];
for n=[30 31 130 131], 
    dof_fix = [dof_fix  locnod(n,:) ];
end
Ndof_fix = size(dof_fix,2);
%b. remove these fixed dof from the list of dof and eliminate them in the matrices
dof_rem = setdiff([1:Ndof],dof_fix);% remaining degrees of freedom
Ndof_rem=Ndof-Ndof_fix;
K=K(dof_rem,dof_rem);   % from here on, K and M are related to the dof_remaining
M=M(dof_rem,dof_rem);

%% Compute eigensolutions
[X,omega2_eig] = eig(full(K),full(M));          % corresponds to K*X = M*X*Omega^2 (our original problem Ax = 0). 
                                           % [X, Omega2] = eig(K, M)  
                                           % <==>   K*X = M*X*Omega2
                                           % 
                                  
[omega2_eig,k] = sort(real(diag(omega2_eig)));      % extracting the real part of the diagonal entries of omega2_eig
X = real(X(:,k));
clear k

omega_eig = sqrt(omega2_eig)./(2*pi);            % Gives all the eigenfrequencies. Converted to Hertz.



%% PART 2: Dynamic analysis. 
% First task: Free vibration
%Compute the first 10 eigenmodes and eigenfrequencies of the system (do not forget
%to apply the boundary conditions from now on). For that, program and apply the
%inverse iteration technique with de
%ation. Verify your algorithm by comparing the
%eigensolutions obtained with the eigensolutions computed from the eig or eigs func-
%tion of Matlab. In particular verify that the accuracy of higher modes is deteriorating
%due to the successive application of de
%ation.

disp('Rank of K: ')
rank(full(K))
disp('Rank of M: ')
rank(full(M))

% As rank(K) = rank(M) = size(M) = size(K) = 240, meaning that M and K consist of 240
% linear independent entries, making them both non-singular.

M_constrained = full(M);        % Mass matrix, with degrees of freedom set (240 x 240)
K_constrained = full(K);        % Stiffness matrix, with degrees of freedom set (240 x 240)



%% Applying the inverse iteration technique 

no_eigenmodes = 10;
omega2_inv = zeros(no_eigenmodes, 1);
omega_inv = zeros(no_eigenmodes, 1);
z = ones(size(M,1), no_eigenmodes);         % z vector from inverse iteration, eq. (3.4.25)
x = ones(size(M,1), no_eigenmodes);         % eigenvectors. Size:(DOF x no_eigenmodes)
P = 1;              % Will be the projection operator. Must first be set to 1 for the first eigenmode. 
it_sum = 0;         %sum of all iterations. Used to find average number of iterations. 
for k = 1:1:no_eigenmodes
    [z(:,k), P, iterations] = inv_iter(K_constrained, M_constrained, x, k, P);
    
    x(:,k) = z(:,k);
    omega2_inv(k,1) = ( x(:,k).' * K * x(:,k) ) / ( x(:,k).' * M * x(:,k) );       % Rayleigh quotient (eq.3.4.15) [rad/s^2]
    omega_inv(k,1) = sqrt(omega2_inv(k,1)) ./ (2 * pi);                       % Converting from [rad/s] to [Hz]
    
    it_sum = it_sum + iterations;
end
it_avg = it_sum / no_eigenmodes;
fprintf('----Average iterations of inv_iter-------\n')
fprintf('\nit_avg =  %d\n', it_avg);
fprintf('-------------------------------\n')

% Accuracy measurement by comparing eigenfrequencies found by inverse
% iteration to the ones found by eigs in MATLAB:
omega_accuracy = zeros(no_eigenmodes,1);
eigenmode_accuracy = zeros(no_eigenmodes,1);
for i = 1:1:no_eigenmodes
   omega_accuracy(i,1) = (omega_eig(i,1) - omega_inv(i,1)) / omega_eig(i,1)  * 100;
   eigenmode_accuracy(i,1) = max(K * x(:,i) - (omega2_inv(i,1) * M * x(:, i)));
end




%% Time integration with implicit Newmark scheme

%Setting initial values:
q0 = 0;         % initial displacement at t = 0
q0_dot = 0;         % initial velocity at t = 0

h_vec = logspace(0, -2, 3);         
%h_vec = [8.957277139407033e-04, 5e-04, 1e-04];        %Used for computational analysis
cost_imp = zeros(size(h_vec,2), 3);         % Vector containing data for computational time. [h, t]

for i_h = 1:1:size(h_vec,2)
    tic
    h = h_vec(i_h);
    gamma_imp = 1/2;
    beta_imp = 1/4;
    C = zeros(size(M,1), size(M,2));

    [q_imp, q_dot_imp, q_dot_dot_imp, h, num_t_entries, t] = Newmark(gamma_imp, beta_imp, K, M, C, h, q0, q0_dot, omega_inv(4,1));
    
    cost_imp(i_h, 1) = h;
    cost_imp(i_h, 2) = num_t_entries;
    cost_imp(i_h, 3) = toc;
    
    
    %% Plotting:
    % Element 10 consists of node 8-9, and dof 43-54. Dof 51 is therefore u3,
    % corresponding to displacement in the z direction in the node where the
    % force is applied.
    
    %% Setting variables used in every plot: %%%%%%%%%%%%%%%%%%%%%%%%%%%
    plotStyle = {'-.','--','-', ':'};
    output_Path = 'C:\Users\adria\Documents\Dokumenter\TUM\Structural Dynamics\Assignment\figures';
    width = 175;                % width of exported plot figure
    height = 100;                   % height of exported plot figure
    lineSize = 2;                 %linesize of plotted lines (linewidth)
    titleSize = 30;
    labelSize = 22;
    valueSize = 100;         %Must be this high to get pdf export right. Originally at 20.
    legendInfo{i_h} = strcat('$h = $', num2str(h), 's');            % https://se.mathworks.com/matlabcentral/answers/29799-adding-a-legend-and-different-coloured-lines-in-a-for-loop
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure(2)
    plot(t, q_imp(51,:), 'k', 'Linestyle', plotStyle{i_h})
    title({'Displacement in Z direction of node subjected to', 'triangular force for varying step size'}, 'Fontsize', titleSize)
    xlabel('Time, $t$  $[s]$', 'FontSize', labelSize, 'Interpreter', 'Latex')
    ylabel('${q}$        $[m]$', 'FontSize', labelSize, 'Interpreter', 'Latex')
    set(gca,'FontSize', valueSize);
    set(gcf, 'PaperPosition', [0 0 width height])
    set(gcf, 'PaperSize', [width height])
    grid on
    grid minor
    hold on
    
    figure(3)
    plot(t, q_dot_imp(51,:), 'k', 'Linestyle', plotStyle{i_h})
    title({'Velocity in Z direction of node subjected to', 'triangular force for varying step size'}, 'Fontsize', titleSize)
    xlabel('Time, $t$  $[s]$', 'FontSize', labelSize, 'Interpreter', 'Latex')
    ylabel('$\dot{q}$       $[\frac{m}{s}]$', 'FontSize', labelSize, 'Interpreter', 'Latex')
    set(gca,'FontSize', valueSize);
    set(gcf, 'PaperPosition', [0 0 width height])
    set(gcf, 'PaperSize', [width height])
    grid on
    grid minor
    hold on
    
    figure(4)
    plot(t, q_dot_dot_imp(51,:), 'k', 'Linestyle', plotStyle{i_h}, 'Linewidth', 1.3)
    title({'Acceleration in Z direction of node subjected to', 'triangular force for varying step size'}, 'Fontsize', titleSize)
    xlabel('Time, $t$  $[s]$', 'FontSize', labelSize, 'Interpreter', 'Latex')
    ylabel('$\ddot{q}$     $[\frac{m}{s^2}]$', 'FontSize', labelSize, 'Interpreter', 'Latex')
    set(gca,'FontSize', valueSize);
    set(gcf, 'PaperPosition', [0 0 width height])
    set(gcf, 'PaperSize', [width height])
    grid on
    grid minor
    hold on
    
end
%Legends and exporting must be applied afterwards due to the for loop
figure(2)
legend(legendInfo, 'Location', 'southeast', 'Interpreter', 'latex')
saveas(figure(2), fullfile(output_Path, 'var_h_disp'), 'pdf')

figure(3)
legend(legendInfo, 'Location', 'southeast', 'Interpreter', 'latex')
saveas(figure(3), fullfile(output_Path, 'var_h_vel'), 'pdf')

figure(4)
legend(legendInfo, 'Location', 'southeast', 'Interpreter', 'latex')
saveas(figure(4), fullfile(output_Path, 'var_h_acc'), 'pdf')



%% Explicit time integration - central difference
% Defining all needed matrices
gamma_exp = 1/2;
beta_exp = 0;

h_stable = 2 / sqrt(max(omega2_eig));

% Part used for calculation costs:
%{
h_vec = [8.957277139407033e-04, 5e-04, 1e-04];        %REMOVE after computational analysis
cost_exp = zeros(size(h_vec,2), 3);         % Vector containing data for computational time. [h, t]

for i_h = 1:1:size(h_vec,2)
    tic
    h = h_vec(1,i_h)
    [q_exp_stable, q_dot_exp_stable, q_dot_dot_exp_stable, h_exp_stable, num_t_entries, t_exp_stable] = Newmark(gamma_exp, beta_exp, K, M, C, h, q0, q0_dot, omega_inv(4,1));
    cost_exp(i_h, 1) = h_exp_stable;
    cost_exp(i_h, 2) = num_t_entries;
    cost_exp(i_h, 3) = toc;
end
%}

[q_exp_stable, q_dot_exp_stable, q_dot_dot_exp_stable, h_exp_stable, num_t_entries, t_exp_stable] = Newmark(gamma_exp, beta_exp, K, M, C, h_stable, q0, q0_dot, omega_inv(4,1));
[q_exp_unstable, q_dot_exp_unstable, q_dot_dot_exp_unstable, h_exp_unstable, num_t_entries, t_exp_unstable] = Newmark(gamma_exp, beta_exp, K, M, C, 1.01*h_stable, q0, q0_dot, omega_inv(4,1));


figure(5)
subplot(2,1,1)
plot(t_exp_stable, q_exp_stable(51, :), 'k', 'Linewidth', lineSize)
title({'Stable explicit Newmark integration', 'with the central difference algorithm'}, 'Fontsize', titleSize)
xlabel('Time, $t$  $[s]$', 'FontSize', labelSize, 'Interpreter', 'Latex')
ylabel('${q}$       $[m]$', 'FontSize', labelSize, 'Interpreter', 'Latex')
h_legend = legend('$h = h_{stable}$' , 'Location', 'south');     % must be done like this... https://se.mathworks.com/matlabcentral/newsreader/view_thread/254118 
set(h_legend, 'Interpreter', 'Latex')
set(gca,'FontSize', valueSize);
set(gcf, 'PaperPosition', [0 0 width height])
set(gcf, 'PaperSize', [width height])
grid on
grid minor


subplot(2,1,2)
plot(t_exp_unstable, q_exp_unstable(51, :), 'k', 'Linewidth', lineSize)
title({'Unstable explicit Newmark integration', 'with the central difference algorithm'}, 'Fontsize', titleSize)
xlabel('Time, $t$  $[s]$', 'FontSize', labelSize, 'Interpreter', 'Latex')
ylabel('${q}$       $[m]$', 'FontSize', labelSize, 'Interpreter', 'Latex')
h_legend = legend('$h = 1.01 \cdot h_{stable}$', 'Location', 'south');
set(h_legend, 'Interpreter', 'Latex')
set(gca,'FontSize', valueSize);
set(gcf, 'PaperPosition', [0 0 width height])
set(gcf, 'PaperSize', [width height])
grid on
grid minor

saveas(figure(5), fullfile(output_Path, 'exp_Newmark'), 'pdf')



%% Plotting displacements of the hangar

figure(6)           % Plotting the computed displacement from the applied load
q_plot = zeros(Ndof,Ndof_rem);
q_plot(dof_rem, 1) = q_exp_stable(:,1000);
plotmesh_no_nodes(Nodes,Elements);
plotdis(Nodes,Elements,locnod,q_plot(:,1), 500);
set(gcf, 'PaperPosition', [0 0 width height])
set(gcf, 'PaperSize', [width height])
saveas(figure(6), fullfile(output_Path, 'force_displacement'), 'pdf')



%% Setting variables used in 10 eigenmodes plot: %%%%%%%%%%%%%%%%%%%%%%%%%%%
plotStyle = {'-.','--','-', ':'};
output_Path = 'C:\Users\adria\Documents\Dokumenter\TUM\Structural Dynamics\Assignment\figures';
width_eig = 100;                % width of exported plot figure
height_eig = 175;                   % height of exported plot figure
lineSize = 2;                 %linesize of plotted lines (linewidth)
titleSize = 30;
labelSize = 22;
valueSize = 100;         %Must be this high to get pdf export right. Originally at 20.
legendInfo{i_h} = strcat('$h = $', num2str(h), 's');            % https://se.mathworks.com/matlabcentral/answers/29799-adding-a-legend-and-different-coloured-lines-in-a-for-loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Xplot=zeros(Ndof,Ndof_rem); % express the modes in the non-fixed numbering used in locnod
Xplot(dof_rem,:)=X;

figure(7) % all 10 firstcomputed eigenmodes X
for i = 1:1:10
    subplot(5,2,i)
    plotdis(Nodes,Elements,locnod,Xplot(:,i),5);
end

set(gcf, 'PaperPosition', [0 0 width_eig height_eig])
set(gcf, 'PaperSize', [width_eig height_eig])
saveas(figure(7), fullfile(output_Path, 'eig_plot'), 'pdf')
% If this section only is run (mark and F9), the plot figures all get the
% same size, which it doesn't if you run the whole script.







