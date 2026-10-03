%%
clear; clc;

syms t beta(t) theta(t) alpha(t) g
syms mM1 mM2 mDiscPully mHoleDisc mHoleMasses
syms lM1 lM2 lDisc lPully
syms rM1 rM2 rDisc rHoleMasses rHoleDisc
syms lsmall llarge pDisc pPully

%Rotation matricies
R10 = [cos(beta(t))  -sin(beta(t))  0;
       sin(beta(t))   cos(beta(t))  0;
       0              0             1];


R21 = [1  0             0;
       0  cos(theta(t))  -sin(theta(t));
       0  sin(theta(t))  cos(theta(t))];


R32 = [ cos(alpha(t))  0  sin(alpha(t));
        0              1  0;
       -sin(alpha(t))  0  cos(alpha(t))];

%full rotations to inertial frame
R20 = R10 * R21;           
R30 = R10 * R21 * R32;     

%Distance vectors
r_DB = [0; pDisc;  0];          
r_PB = [0; pPully; 0];          
r_DP = (r_DB + r_PB) / 2; %not using COM this is probably a source of some error   

r_SB = [0; lsmall; 0];          
r_LB = [0; llarge; 0];          

%Moving all vectors into the intertial frame
r_DP_0 = R30 * r_DP;
r_SB_0 = R20 * r_SB;
r_LB_0 = R20 * r_LB;


%angular velocities
w1   = [0; 0; diff(beta,  t)];  
w21  = [diff(theta, t); 0; 0];   
w32  = [0; diff(alpha, t); 0];   

w2_in_2 = R21.' * w1 + w21;              
w3_in_3 = R32.' * w2_in_2 + w32;         

%ang veloc in reference to inertial
w2_in_0 = R20 * w2_in_2
w3_in_0 = R30 * w3_in_3


%Inertia
%inertia m1 (large)
IM1 = diag([(1/4*mM1*rM1^2 + 1/12*mM1*lM1^2);
             (1/2*mM1*rM1^2);
             (1/4*mM1*rM1^2 + 1/12*mM1*lM1^2)]);

IM1_hole = diag([(1/4*mHoleMasses*rHoleMasses^2 + 1/12*mHoleMasses*lM1^2);
                  (1/2*mHoleMasses*rHoleMasses^2);
                  (1/4*mHoleMasses*rHoleMasses^2 + 1/12*mHoleMasses*lM1^2)]);

IM1_net = IM1 - IM1_hole;% hollow M1, at COM, in frame 2

%Inertia M2 (small)
IM2 = diag([(1/4*mM2*rM2^2 + 1/12*mM2*lM2^2);
             (1/2*mM2*rM2^2);
             (1/4*mM2*rM2^2 + 1/12*mM2*lM2^2)]);

IM2_hole = diag([(1/4*mHoleMasses*rHoleMasses^2 + 1/12*mHoleMasses*lM2^2);
                  (1/2*mHoleMasses*rHoleMasses^2);
                  (1/4*mHoleMasses*rHoleMasses^2 + 1/12*mHoleMasses*lM2^2)]);

IM2_net = IM2 - IM2_hole;% hollow M2, at COM, in frame 2

%Inertia Disc
IDP_full = diag([(1/4*mDiscPully*rDisc^2 + 1/12*mDiscPully*(lDisc+lPully)^2);
                  (1/2*mDiscPully*rDisc^2);
                  (1/4*mDiscPully*rDisc^2 + 1/12*mDiscPully*(lDisc+lPully)^2)]);

IDP_hole = diag([(1/4*mHoleDisc*rHoleDisc^2 + 1/12*mHoleDisc*(lDisc+lPully)^2);
                  (1/2*mHoleDisc*rHoleDisc^2);
                  (1/4*mHoleDisc*rHoleDisc^2 + 1/12*mHoleDisc*(lDisc+lPully)^2)]);


IDP_net = IDP_full - IDP_hole;% hollow disc-pulley, at COM, in frame 3

%setup parrallel axis
PAT = @(I, m, r) I + m * (dot(r,r)*eye(3) - r*r.');

%shift all to point B (PRM point)
IM1_B = PAT(IM1, mM1, r_LB);
IM1_hole_B = PAT(IM1_hole, mHoleMasses, r_LB);
ITM1  = IM1_B - IM1_hole_B;% net hollow M1 at B, frame 2


IM2_B = PAT(IM2, mM2, r_SB);
IM2_hole_B = PAT(IM2_hole, mHoleMasses, r_SB);
ITM2 = IM2_B - IM2_hole_B;% net hollow M2 at B, frame 2


IDP_B = PAT(IDP_full, mDiscPully, r_DP);
IDP_hole_B = PAT(IDP_hole, mHoleDisc,  r_DP);
ITDP = IDP_B - IDP_hole_B;% net hollow disc-pulley at B, frame 3

% ang momentum
HM1_2  = ITM1 * w2_in_2;     
HM2_2  = ITM2 * w2_in_2;    
HDP_3  = ITDP * w3_in_3;     
%HM1_2  = IM1_net * w2_in_2;     
%HM2_2  = IM2_net * w2_in_2;    
%HDP_3  = IDP_net * w3_in_3;     

%ang momentum for M derivation
HM1_2_dot = diff(HM1_2, t) + cross(w2_in_2, HM1_2);
HM2_2_dot = diff(HM2_2, t) + cross(w2_in_2, HM2_2);
HDP_3_dot = diff(HDP_3, t) + cross(w3_in_3, HDP_3);


HMass_2_dot = HM1_2_dot + HM2_2_dot;

%Linear velcoties in inertial frame
r_DP_0_dot  = diff(r_DP_0, t);
r_SB_0_dot  = diff(r_SB_0, t);
r_LB_0_dot  = diff(r_LB_0, t);

% Linear acceleration in inertial frame
r_DP_0_ddot = diff(r_DP_0_dot, t);
r_SB_0_ddot = diff(r_SB_0_dot, t);
r_LB_0_ddot = diff(r_LB_0_dot, t);

g_vec = [0; 0; -g];


F2 = mDiscPully * r_DP_0_ddot - R30.'*g_vec;


F_M1 = mM1 * g_vec;
F_M2 = mM2 * g_vec;


FPivot_0 = -(F2 + F_M1 + F_M2);

g_vec_3    = R30.' * g_vec; % F2 expressed in frame 3
M_F2_3     = cross(r_DP, mDiscPully*g_vec_3);% moment of F2 about B, frame 3
M_joint_3  = HDP_3_dot - M_F2_3;% joint moment, frame 3
M_joint_2  = R32 * M_joint_3;% joint moment, frame 2

F_M1_in_2  = R20.' * F_M1;% gravity on M1 in frame 2
F_M2_in_2  = R20.' * F_M2;% gravity on M2 in frame 2
M_grav_2   = cross(-r_LB, F_M1_in_2) + cross(-r_SB, F_M2_in_2);

%Equation of motion = momentum about point B 
EOM = HMass_2_dot + M_grav_2 + M_joint_2

syms theta_new theta_dot theta_ddot beta_new beta_dot beta_ddot alpha_new alpha_dot alpha_ddot

old_syms = [diff(theta,t,2), diff(theta,t), theta, diff(beta,t,2), diff(beta,t), beta, diff(alpha, t, 2), diff(alpha,t), alpha];
new_syms = [theta_ddot, theta_dot, theta_new, beta_ddot, beta_dot, beta_new, alpha_ddot, alpha_dot, alpha_new];

eqns = subs(EOM, old_syms, new_syms);   

% Extract mass matrix and forcing vector
[A_raw, rhs_raw] = equationsToMatrix(eqns, [beta_ddot, theta_ddot, alpha_ddot]);
b_raw = rhs_raw;  


%trying to make the equations less massive but this takes forever
% Simplify element-by-element (faster than simplify on full matrix)
%{
fprintf('Simplifying mass matrix A (3x3)...\n');
A_sym = sym(zeros(3,3));
for i = 1:3
    for j = 1:3
        A_sym(i,j) = simplify(A_raw(i,j), 'Steps', 50);
    end
end

fprintf('Simplifying forcing vector b (3x1)...\n');
b_sym = sym(zeros(3,1));
for i = 1:3
    b_sym(i) = simplify(b_raw(i), 'Steps', 50);
end

% Collect terms by velocity products for readable display
vel_syms = [beta_dot, theta_dot, alpha_dot];
A_disp = A_sym;  b_disp = b_sym;
for i = 1:3
    for j = 1:3
        A_disp(i,j) = collect(A_sym(i,j), vel_syms);
    end
    b_disp(i) = collect(b_sym(i), vel_syms);
end

disp('=== Mass matrix A(q) ==='); disp(A_disp)
disp('=== Forcing vector b(q, q_dot) ==='); disp(b_disp)
%}


Xdd_sym = simplify(A_raw \ b_raw); %solving for the three accelerations 

beta_dd_sol = simplify(Xdd_sym(1))
theta_dd_sol = simplify(Xdd_sym(2))
alpha_dd_sol = simplify(Xdd_sym(3))
%%
%masses
mM1_val = 0.9;
mM2_val = 0.03;
mDiscPully_val = 1.47;
mRod_val = 0.0;

%Lengths (m)
lM1_val = 0.0313;
lM2_val = 0.019;
lDisc_val = 0.0221;
lPully_val = 0.022;

%Radius (m)
rM1_val = 0.07/2;
rM2_val = 0.044/2;
rDisc_val = 0.254/2;
rPully_val = 0.0582/2;
rHoleMasses_val = 0.0128;
rHoleDisc_val = 0.0093;

pDisc_val = (12.8+2.3+(2.21/2))*0.01;
pPully_val = (12.8+(2.3/2))*0.01;
lsmall_val = (-28.2)*0.01;
llarge_val = (-31.0)*0.01;
g_var = 9.80665



VIDEO = 1;  
videoFile = 'Group_Gyro_Animation.mp4';

%Calculating mass of holes
volumeHoleDisc = (pi*rHoleDisc_val.^2*((lDisc_val+lPully_val)));
densityDisc = mDiscPully_val / (pi*rDisc_val.^2*lDisc_val);

densityMasses =  mM1_val / (pi*rM1_val.^2*abs(lM1_val));
volumeHoleMasses = (pi*rHoleMasses_val.^2*(abs(lM1_val)+abs(lM2_val))/2);

mHoleDisc_val = densityDisc * volumeHoleDisc
mHoleMasses_val = densityMasses * volumeHoleMasses


disp('=== Running Section ==='); 
param_old = [mM1 mM2 mDiscPully ...
             lM1 lM2 lDisc lPully ...
             rM1 rM2 rDisc rHoleDisc ...
             pDisc pPully lsmall llarge g ...
             mHoleDisc mHoleMasses rHoleMasses];

param_new = [mM1_val mM2_val mDiscPully_val ...
             lM1_val lM2_val lDisc_val lPully_val ...
             rM1_val rM2_val rDisc_val rHoleDisc_val ...
             pDisc_val pPully_val lsmall_val llarge_val g_var ...
             mHoleDisc_val mHoleMasses_val rHoleMasses_val];


beta_dd_sol  = simplify(subs(beta_dd_sol,  param_old, param_new));
theta_dd_sol = simplify(subs(theta_dd_sol, param_old, param_new));
alpha_dd_sol = simplify(subs(alpha_dd_sol, param_old, param_new));

symvar(beta_dd_sol)

beta0 = 0;
theta0 = pi/12;
alpha0 = 0;

beta_dot0 = 0;
theta_dot0 = 0.0;
alpha_dot0 = 500*2*pi/60; %500 was the measured RPM

X_init = [beta0; theta0; alpha0; beta_dot0; theta_dot0; alpha_dot0];

vars = {'beta_new','theta_new','alpha_new','beta_dot','theta_dot','alpha_dot'};
beta_dd_fun  = matlabFunction(beta_dd_sol,  'Vars', vars);
theta_dd_fun = matlabFunction(theta_dd_sol, 'Vars', vars);
alpha_dd_fun = matlabFunction(alpha_dd_sol, 'Vars', vars);


tspan = [0 90];
options = odeset('RelTol',1e-7,'AbsTol',1e-7);

sol = ode45(@(t,X)eom(t,X,beta_dd_fun,theta_dd_fun,alpha_dd_fun), tspan, X_init, options);

dt = 0.05;
t = tspan(1):dt:tspan(2);
X = deval(sol,t);

beta_dd_out  = zeros(1, length(t));
theta_dd_out = zeros(1, length(t));
alpha_dd_out = zeros(1, length(t));

for i = 1:length(t)
    b   = X(1,i);  th  = X(2,i);  al  = X(3,i);
    bd  = X(4,i);  thd = X(5,i);  ald = X(6,i);

    % Same singularity guard as in eom()
    if abs(cos(th)) < 1e-3
        th = th + sign(cos(th) + eps) * 1e-3;
    end

    beta_dd_out(i)  = beta_dd_fun( b, th, al, bd, thd, ald);
    theta_dd_out(i) = theta_dd_fun(b, th, al, bd, thd, ald);
    alpha_dd_out(i) = alpha_dd_fun(b, th, al, bd, thd, ald);
end

figure('Color','w');
plot(t,X,'LineWidth',1.2)
xlabel('time (s)')
ylabel('states')
h = legend('$\phi$','$\theta$','$\psi$', ...
           '$\dot{\phi}$','$\dot{\theta}$','$\dot{\psi}$');
set(h,'Interpreter','latex')
grid on
title('Gyroscope states')


if VIDEO
    fps = 15;
    MyVideo = VideoWriter(videoFile,'MPEG-4');
    MyVideo.FrameRate = fps;
    open(MyVideo);
end


fig = figure('Color','w');

for i = 1:length(t)
    clf(fig)
    hold on; grid on; axis equal;

    beta = X(1,i);
    theta = X(2,i);
    alpha = X(3,i);

    draw_gyro(beta,theta,alpha);

    title(sprintf('Gyroscope animation, t = %.2f s',t(i)))
    xlabel('X'); ylabel('Y'); zlabel('Z');
    view(35,25)
    axis([-0.6 0.6 -0.6 0.6 0 0.7])
    drawnow;

    if VIDEO
        frame = getframe(gcf);
        writeVideo(MyVideo,frame);
    else
        pause(dt)
    end
end

if VIDEO
    close(MyVideo)
    disp(['Saved video: ', videoFile])
end
accel  = zeros(3,length(t));

function Xdot = eom(~, X, f1, f2, f3)

    beta  = X(1);
    theta = X(2);
    alpha = X(3);

    beta_dot  = X(4);
    theta_dot = X(5);
    alpha_dot = X(6);

    if abs(cos(theta)) < 1e-3
        theta = theta + sign(cos(theta)+eps)*1e-3;
        
    end

    %theta = max(theta, 1e-4);
    %theta = min(theta, pi/2 - 1e-4);
   
    beta_dd  = f1(beta,theta,alpha,beta_dot,theta_dot,alpha_dot);
    theta_dd = f2(beta,theta,alpha,beta_dot,theta_dot,alpha_dot);
    alpha_dd = f3(beta,theta,alpha,beta_dot,theta_dot,alpha_dot);

    Xdot = [
        beta_dot;
        theta_dot;
        alpha_dot;
        beta_dd;
        theta_dd;
        alpha_dd
    ];

end

function draw_gyro(phi,theta,psi,p)

O = [0;0;0.35];

R01 = Rz(phi);
R12 = Rx(theta);
R23 = Ry(psi);

R02 = R01*R12;
R03 = R02*R23;

% Inertial axes
plot3([0 0.15],[0 0],[0 0],'r','LineWidth',2)
plot3([0 0],[0 0.15],[0 0],'g','LineWidth',2)
plot3([0 0],[0 0],[0 0.15],'b','LineWidth',2)

% Vertical support
plot3([0 0],[0 0],[-0.25 O(3)],'k','LineWidth',4)

% Bar along y2
sMin = -0.35;
sMax = 0.22;

barLocal = [0 0;
            sMin sMax;
            0 0];

barWorld = O + R02*barLocal;

plot3(barWorld(1,:),barWorld(2,:),barWorld(3,:), ...
      'Color',[0.6 0.35 0.1],'LineWidth',5)

% Rotor disk
rotorCenter = O + R02*[0;(0.128 + 0.020 + 0.023/2);0];
draw_disk(rotorCenter,R03,0.122,0.023,[0.75 0.75 0.75]);

% Big counterweight
bigCenter = O + R02*[0;-0.310;0];
draw_sphere(bigCenter,0.028,[0.2 0.3 0.8]);

% Small counterweight
smallCenter = O + R02*[0;-0.282;0];
draw_sphere(smallCenter,0.015,[0.4 0.7 0.9]);

end

%% ============================================================
% DRAW DISK
% ============================================================
function draw_disk(center,R,radius,thickness,color)

n = 60;
theta = linspace(0,2*pi,n);
rad = linspace(0,radius,2);
[Theta,Rad] = meshgrid(theta,rad);

% Disk spin axis is local y-axis.
% Disk face lies in local x-z plane.
X = Rad.*cos(Theta);
Z = Rad.*sin(Theta);
Yfront = thickness/2*ones(size(X));
Yback  = -thickness/2*ones(size(X));

P1 = R*[X(:)';Yfront(:)';Z(:)'];
P2 = R*[X(:)';Yback(:)'; Z(:)'];

X1 = reshape(P1(1,:),size(X)) + center(1);
Y1 = reshape(P1(2,:),size(X)) + center(2);
Z1 = reshape(P1(3,:),size(X)) + center(3);

X2 = reshape(P2(1,:),size(X)) + center(1);
Y2 = reshape(P2(2,:),size(X)) + center(2);
Z2 = reshape(P2(3,:),size(X)) + center(3);

if all(isfinite(X1(:))) && all(isfinite(Y1(:))) && all(isfinite(Z1(:)))
    surf(X1,Y1,Z1,'FaceColor',color,'EdgeColor','none','FaceAlpha',0.9);
end

if all(isfinite(X2(:))) && all(isfinite(Y2(:))) && all(isfinite(Z2(:)))
    surf(X2,Y2,Z2,'FaceColor',color,'EdgeColor','none','FaceAlpha',0.9);
end

% Stripe to show disk spin
stripe = [radius -radius;
          0       0;
          0       0];

stripeW = center + R*stripe;
plot3(stripeW(1,:),stripeW(2,:),stripeW(3,:),'k','LineWidth',2)

end

%% ============================================================
% DRAW SPHERE
% ============================================================
function draw_sphere(center,radius,color)

[x,y,z] = sphere(30);

Xs = radius*x + center(1);
Ys = radius*y + center(2);
Zs = radius*z + center(3);

if all(isfinite(Xs(:))) && all(isfinite(Ys(:))) && all(isfinite(Zs(:)))
    surf(Xs,Ys,Zs, ...
         'FaceColor',color, ...
         'EdgeColor','none', ...
         'FaceAlpha',0.95);
end

end

%% ============================================================
% ROTATION MATRICES
% ============================================================
function R = Rz(a)
R = [cos(a)  sin(a) 0;
     -sin(a)  cos(a) 0;
     0       0      1];
end

function R = Rx(a)
R = [1 0 0;
     0 cos(a) -sin(a);
     0 sin(a)  cos(a)];
end

function R = Ry(a)
R = [cos(a) 0 sin(a);
     0      1 0;
     -sin(a) 0 cos(a)];
end


%% MODEL BODY ANGULAR VELOCITY AND ACCELERATION AT SENSOR
close all;
tspan = [0 90];
dt = 0.05;
t = tspan(1):dt:tspan(2);

r_sensor = [0; -0.15; 0.05];  % sensor location in frame 2, 150mm along -y2

a_B       = zeros(3, length(t));
omega = zeros(3, length(t));

for i = 1:length(t)
    beta = X(1,i);
    theta  = X(2,i);

    R21_num = Rx(theta);
    R10_num = Rz(beta);
    R20_num = R10_num * R21_num;

    beta_d  = X(4,i);  
    theta_d = X(5,i);
    beta_dd = beta_dd_out(i);  
    theta_dd = theta_dd_out(i);

    %ang val frame 2
    w2 = [theta_d;
          beta_d*sin(theta);
          beta_d*cos(theta)];
    omega(:,i) = w2;

    %accel in frame 2
    alph2 = [theta_dd;
             beta_dd*sin(theta) + beta_d*theta_d*cos(theta);
             beta_dd*cos(theta) - beta_d*theta_d*sin(theta)];

    %move acceleration from PRM point to sensor mounting loc
    a_B(:,i) = cross(alph2, r_sensor) + cross(w2, cross(w2, r_sensor));

    a_B(:,i) = a_B(:,i) - R20_num.' * [0; 0; -9.81];%accout for gravity 
end 


T = readtable('ExportedData - Run 2.csv', 'VariableNamingRule', 'preserve');

AccelX = T.("Acceleration - x (m/s²) Run 1");
AccelY = T.("Acceleration - y (m/s²) Run 1");
AccelZ = T.("Acceleration - z (m/s²) Run 1");

AngVelX = T.("Angular Velocity - x (rad/s) Run 1");
AngVelY = T.("Angular Velocity - y (rad/s) Run 1");
AngVelZ = T.("Angular Velocity - z (rad/s) Run 1");


%MAgnitude
omega_mag = vecnorm(omega);
meas_mag  = sqrt(AngVelX.^2 + AngVelY.^2 + AngVelZ.^2);

figure(11);
plot(t, omega_mag); hold on;
plot(t, meas_mag(1:length(t)));
legend('model','measured');
grid on;


%X 
figure(4);
plot(t, omega(2,:), 'LineWidth', 1.5); hold on;
plot(t, AngVelX(1:length(omega(1,:)),:), 'LineWidth', 1.5);

xlabel('Time (s)');
ylabel('Angular Velocity');
legend('Calculated','Measured');
title('Measured Angular Velocity X Vs Calculated');
grid on;

%Y
figure(5);
plot(t, omega(1,:), 'LineWidth', 1.5); hold on;
plot(t, AngVelY(1:length(omega(1,:)),:), 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Angular Velocity');
legend('Calculated','Measured');
title('Measured Angular Velocity Y Vs Calculated');
grid on;

%Z
figure(6);
plot(t, -omega(3,:), 'LineWidth', 1.5); hold on;
plot(t, AngVelZ(1:length(omega(1,:)),:), 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Angular Velocity');
legend('Calculated','Measured');
title('Measured Angular Velocity Z Vs Calculated');
grid on;

%A Acel
meanX = mean(AccelX(1:length(omega(2,:)),:));
meanZ = mean(AccelZ(1:length(omega(2,:)),:));
figure(7);
plot(t, -a_B(2,:), 'LineWidth', 1.5); hold on;
plot(t, detrend(AccelX(1:length(omega(2,:)),:))+meanX, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Angular Acceleration (rad/s^2)');
legend('Calculated','Measured');
title('Measured Angular Acceleration X Vs Calculated');
grid on;

%Y Acel

figure(8);
plot(t, a_B(1,:), 'LineWidth', 1.5); hold on;
plot(t, AccelY(1:length(omega(2,:)),:), 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Angular Acceleration (rad/s^2)');
legend('Calculated','Measured');
title('Measured Angular Acceleration Y Vs Calculated');
grid on;

%Z ACELm
figure(9);
plot(t, a_B(3,:), 'LineWidth', 1.5); hold on;
plot(t, detrend(AccelZ(1:length(omega(2,:)),:))+meanZ, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Angular Acceleration (rad/s^2)');
legend('Calculated','Measured');
title('Measured Angular Acceleration Z Vs Calculated');
grid on;

figure(69)
hold on;
plot(t, X(3,:), 'LineWidth', 1.5); hold on;
plot(t, AngVelZ(1:length(omega(1,:)),:), 'LineWidth', 1.5);
plot(t, X(2,:), 'LineWidth', 1.5); hold on;
plot(t, AngVelY(1:length(omega(1,:)),:), 'LineWidth', 1.5);
plot(t, X(1,:), 'LineWidth', 1.5); hold on;
plot(t, AngVelX(1:length(omega(1,:)),:), 'LineWidth', 1.5);
legend('alpha','AngZ', 'theta', 'AngY', 'beta', "AngX");

hold off;

figure(98)
hold on;

plot(t, beta_dd_out, 'LineWidth', 1.5); hold on;

plot(t, AngVelZ(1:length(omega(1,:)),:), 'LineWidth', 1.5);
plot(t, theta_dd_out, 'LineWidth', 1.5); hold on;
plot(t, AngVelY(1:length(omega(1,:)),:), 'LineWidth', 1.5);

legend('betadot', 'AZ', 'thetadot', 'aY');
hold off;
