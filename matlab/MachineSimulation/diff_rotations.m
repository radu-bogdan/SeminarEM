%% This script solves the magnetostatic problem for different
% rotation angles of the rotor. To show an animation of the different
% configurations, type movie(anim) (see documentation for more details)
% into the comand window after executing this script. A data structure with
% 128 different angles is already stored in
% Animations/diff_angles_with_current.
%
% To show the animation, type e.g. movie(anim,1,4) into the command window
% (see also doc movie)

clear;

%% loading data
% -----------------------------------------------------------------------
% If the data in DataPreCalculation/motor.mat changed, you have to rerun
% DataPreCalculation/CreateData.m !!!!
% -----------------------------------------------------------------------
addpath(genpath('DataPreCalculation'));
addpath(genpath('Material'));
addpath(genpath('Plotting'));
addpath(genpath('Assemble'));
addpath(genpath('Linesearch'));
addpath(genpath('Meshing'));
load -mat motor_divided.mat % datastructure with double amount of points on the interface.
load -mat J3_vec.mat
load -mat M_vec.mat
load -mat material.mat
load -mat e_stator_outer.mat


%% settings
gamma = 10^5; % penalty for dirichlet boundary condition
Nitsche_pen = 1e7; %penalty factor for Nitsche-type mortaring
num_iteration = 50; %max number of iterations for newton method
eps = 10^-4; % tolerance for abortion criterium
J3_vec = 0*J3_vec; % with or without current (comment out: with current)
nRot = 32; % nRot rotations with angle 2pi/nRot


%% finding points of rotor
magnet_ind = find(contains(regions_2d,'magnet'));
rotorParts_ind = find(contains(regions_2d,'rotor'));
air_gap_ind = find(matches(regions_2d,'air_gap'));
shaft_ind = find(contains(regions_2d,'shaft'));
rotor_ind = [magnet_ind,rotorParts_ind,air_gap_ind,shaft_ind];

t_ind = ismember(t(4,:),rotor_ind);
t_rotor = t(:,t_ind);
p_rotor_ind = unique(t_rotor(1:3,:));
p_rotor_ref = p(:,p_rotor_ind);
M_vec_ref = M_vec;


%% 
% size of datastructures
np=size(p,2); nt=size(t,2); ne = size(e,2);
uh = zeros(np,1);        %initialize solution

% variable for creating animation afterwards
anim(nRot) = struct('cdata',[],'colormap',[]);
angle_ref = 2*pi/nRot;
for i = 1:nRot
    tic
    % rotating rotor
    angle = angle_ref*(i-1);
    Rot = [cos(angle), -sin(angle); sin(angle), cos(angle)];
    p(:,p_rotor_ind) = Rot*p_rotor_ref;
    
    % calculating the edges on the interface between rotor and stator
    eI = calcEdgeInterfaceCircle(p,eI1,eI2,t);
    
    M_vec = Rot*M_vec_ref;
    
    %plot mesh
    % figure
%     h = plotMesh(p,t,e,material);
%     set(h, 'LineStyle', 'none');
             
    
    %Caculate constant Matrices for iteration
    B = calcB(p,t,e);
    R_mat = Assemble_Robin_Matrix(p,e,B);      %calculate Robin Matrix
    
    [S1,S2] = Assemble_Interface_Matrices(p,eI,t,material.nu0,Nitsche_pen);
    % symmetric Nitsche-type mortaring
    S = S1 - S2 - S2';
%     % non-symmetric Nitsche-type mortaring
%     S = S1 - S2 + S2';
    
    [J_j,J_m] = Assemble_J_Vec(p,t,J3_vec,M_vec,B);
    J_vec = J_j - material.nu0*J_m;
    
    alpha = 1;
    vh = randn(np,1);
    for iter = 1:1:num_iteration
        %solve (N(n) + R)wh = J - F1(n) - F2(n)
    
        [duhx,duhy] = pdegrad(p,t,uh);
        grad_uh = [duhx;duhy];
        hessian = hessian_g(t(4,:),grad_uh,material);
        g_grad_val = grad_g(t(4,:),grad_uh,material);
    
        [F1_vec,F2_vec] = Assemble_F_Vec(p,t,e,uh,B,g_grad_val,M_vec,material);  
        N_mat = Assemble_Newton_Matrix(p,t,uh,B,hessian,material);
        Suh = S*uh;
        F = J_vec - gamma*F1_vec - F2_vec - Suh;
        A = N_mat + gamma*R_mat + S;
        wh = A\F;
        alpha = alpha*2;
        phi = @(alpha) Functional_Nitsche(t,uh + alpha*wh,gamma,J_j,J_m,F1_vec,Suh,material,B);
        phi_der0 = Functional_der0_Nitsche(wh,gamma,J_j,J_m,F1_vec,F2_vec,Suh,material);
    
        
        [alpha,iter_linesearch,energy] = ArmijoLineSearch(phi,phi_der0,alpha);
        uh = uh + alpha*wh;
        if iter == 1
            res0 = norm(F);
        end
        res = norm(F);
        if res < eps*res0
            break;
        end
    % alternative abortion criteria
    %     if max(abs(wh))/max(abs(uh)) < eps
    %         break;
    %     end
%     fprintf('Newton Iteration %d: energy=%d, alpha=%d\n',iter,energy,alpha);
    end
    fprintf('Summary:\n')
    toc;
    fprintf('Newton Iterations: %d\n',iter);
    
    [ux,uy] = pdegrad(p,t,uh);
    B = [ux;uy];
    t1 = t;
    t1(4,:) = t(4,:) + ones(1,nt);
    % figure
    % h = plotMesh(p,t,e,material);
    % set(h, 'LineStyle', 'none');
    % hold on
    % pdeplot(p,e,t1,'FlowData',[uy;-ux]);
    
    % B = [uy;-ux];
    % Babs = sqrt(B(1,:).^2 +  B(2,:).^2);
    % figure
    % pdesurf(p,t,Babs);
    % colorbar
    % title('Magnetic Flux |B| [T]')
    % colormap jet
    % view([-0.83 90.00])
    
    
    pdeplot(p,e,t,'XYData',uh,Contour="on",ColorMap="jet")
    % adding new frame to animation
    anim(i) = getframe(gcf);
end