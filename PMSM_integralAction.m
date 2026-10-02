clear
clc
close all

% Load torque
% The primary IDA-PBC controller is designed for tau_L = 0.
tau_L = 0*1;

% Simulation length
Tsim = 10;

%% Define the plant
% Plant state:
% x1 = id
% x2 = iq
% x3 = omega

% Motor model and parameters from Petrovic 2001
plant.P   = 2;          % Number of pole pairs
plant.Phi = 0.199;      % Permanent-magnet flux constant [V s]
plant.Rs  = 1.5;        % Stator winding resistance [Ohm]
plant.J   = 1.08e-3;    % Rotor inertia [kg m^2]
plant.Ld  = 12e-3;      % d-axis inductance [H]
plant.Lq  = 6e-3;       % q-axis inductance [H]

% Physical input matrix in the coordinates [id; iq; omega]
%
% Ld*dot(id) = ... + vd
% Lq*dot(iq) = ... + vq
plant.g = [1/plant.Ld, 0;
           0, 1/plant.Lq;
           0, 0];

% Left annihilator of the input matrix
plant.gPerp = [0, 0, 1];

% Input distribution for the load torque
plant.gd = [0; 0; 1/plant.J];

% Open-loop physical energy
plant.H = @(x) 0.5*plant.Ld*x(1)^2 + 0.5*plant.Lq*x(2)^2 + 0.5*(plant.J/plant.P)*x(3)^2;

% Open-loop plant dynamics
%
% Ld*dot(id) = -Rs*id + omega*Lq*iq + vd
% Lq*dot(iq) = -Rs*iq - omega*Ld*id - omega*Phi + vq
% J*dot(omega) = P*((Ld-Lq)*id*iq + Phi*iq) - tau_L
%
plant.dx = @(x,u) [(-plant.Rs*x(1) + x(3)*plant.Lq*x(2) + u(1))/plant.Ld;
                    (-plant.Rs*x(2) - x(3)*plant.Ld*x(1) - x(3)*plant.Phi + u(2))/plant.Lq;
                    (plant.P*((plant.Ld-plant.Lq)*x(1)*x(2) + plant.Phi*x(2)) - tau_L)/plant.J];

%% Primary control design
% IDA-PBC design from Section 5.2 of
% Ortega and Garcia-Canseco (2004)

% IDA-PBC tuning parameters
ctrl.phi1 = 1;
ctrl.phi2 = 1;

ctrl.r1 = 20;
ctrl.r2 = 0.01;

% J12 is a free interconnection parameter
ctrl.J12 = 0;

% J23 must be negative
ctrl.J23 = -1;

% Desired equilibrium for zero load torque
ctrl.id_star = 0;
ctrl.iq_star = 0;
ctrl.omega_star = 31.42;       % [rad/s]

ctrl.x_star = [ctrl.id_star;
               ctrl.iq_star;
               ctrl.omega_star];

% Define function psi
% For tau_L = 0:
% psi'(iq) = -(P*Phi)/(J23*J)*iq
% Hence:
% psi(iq) = -(P*Phi)/(2*J23*J)*iq^2
ctrl.psi = @(iq) -(plant.P*plant.Phi)/(2*ctrl.J23*plant.J)*iq.^2;

% Desired interconnection term J13
ctrl.J13 = @(iq) plant.P/(ctrl.phi1*plant.J)*(plant.Lq-plant.Ld)*iq;

% For the primary control design we assume no load torque.
% The load disturbance will be handled by the secondary IA controller.
ctrl.tau_L_hat = 0;

% % % Primary control inputs
% d-axis physical voltage
ud = @(x) plant.Ld*((plant.Rs/plant.Ld)*x(1) - (plant.Lq/plant.Ld)*x(2)*x(3) ...
    - ctrl.r1*ctrl.phi1*x(1) + ctrl.J12/(ctrl.J23*plant.J)*(ctrl.tau_L_hat-plant.P*plant.Phi*x(2)) ...
    + plant.P*ctrl.phi2/(ctrl.phi1*plant.J)*(plant.Lq-plant.Ld)*x(2)*(x(3)-ctrl.omega_star));
% q-axis physical voltage
uq = @(x) plant.Lq*((plant.Rs/plant.Lq)*x(2) + (plant.Ld/plant.Lq)*x(1)*x(3) ...
    + (plant.Phi/plant.Lq)*x(3) - ctrl.J12*ctrl.phi1*x(1) ...
    - ctrl.r2/(ctrl.J23*plant.J)*(ctrl.tau_L_hat-plant.P*plant.Phi*x(2)) ...
    + ctrl.J23*ctrl.phi2*(x(3)-ctrl.omega_star));
% Total primary control input
ctrl.u = @(x) [ud(x); uq(x)];

% Closed-loop Hamiltonian from primary control design
ctrl.Hd = @(x) 0.5*ctrl.phi1*x(1).^2 + ctrl.psi(x(2)) + 0.5*ctrl.phi2*(x(3)-ctrl.omega_star).^2;

% Closed-loop Hamiltonian gradient
ctrl.dHd = @(x) [ctrl.phi1*x(1);
                -(plant.P*plant.Phi)/(ctrl.J23*plant.J)*x(2);
                ctrl.phi2*(x(3)-ctrl.omega_star)];

% Interconnection and damping structure after primary control design
ctrl.Fd = @(x) [ ...
    -ctrl.r1,             ctrl.J12,          ctrl.J13(x(2));
    -ctrl.J12,           -ctrl.r2,           ctrl.J23;
    -ctrl.J13(x(2)),     -ctrl.J23,          0
    ];

%% Integral action control design

% Integral-action gain
int.Ki = 1e2*eye(2);

% Signal to be integrated
int.Ao = [1, 0, 0;
          0, 0, ctrl.J23];

% Integral dynamics
% z1_dot = phi1*id
% z2_dot = J23*phi2*(omega-omega_star)
int.dz = @(x) int.Ao*ctrl.dHd(x);

% % % Matrices required by the IA controller
% Moore-Penrose left inverse of g
int.gDag = (plant.g.'*plant.g)\plant.g.';

% Free damping/interconnection terms
int.Lambda1 = zeros(2);
int.Lambda2 = @(x) int.gDag*(ctrl.Fd(x)+ctrl.Fd(x).')*plant.gPerp.'/(plant.gPerp*plant.gPerp.');

% Complete Lambda matrix
% Lambda = [Lambda1 Lambda2]*[g';gPerp]
int.Lambda = @(x) [int.Lambda1, int.Lambda2(x)]*[plant.g.'; plant.gPerp];

% gamma function
int.gamma = @(x) [0; x(2)];
int.dgamma_dx = [0, 0, 0;
                 0, 1, 0];

% Free F12 component
int.F12_1 = -int.gDag*int.Ao.';

% Integral control law
int.v = @(x,z) (-int.F12_1 + int.gDag*ctrl.Fd(x)*int.dgamma_dx.')*int.Ki*(int.gamma(x)-z) ...
                - int.Lambda(x)*(ctrl.dHd(x) + int.dgamma_dx.'*int.Ki*(int.gamma(x)-z));

% Integral energy
int.Hia = @(x,z) 0.5*(int.gamma(x)-z).'*int.Ki*(int.gamma(x)-z);

%% Simulate the system
% Initial conditions
id_0    = 0;       % [A]
iq_0    = 0;       % [A]
omega_0 = 0;       % [rad/s]
plant.x0 = [id_0;
            iq_0;
            omega_0];
int.z0 = [0;0];
sim.x0 = [plant.x0; int.z0];

% Full system ODE
sim.ode = @(t,x) [plant.dx(x(1:3), ctrl.u(x(1:3)) + int.v(x(1:3),x(4:5)));
                int.dz(x(1:3))];

% Simulate ODE
options = odeset('RelTol',1e-6);
[res.t,res.w] = ode23s(sim.ode,[0 Tsim],sim.x0,options);

%% Plot results
% Extract states
res.x = res.w(:,1:3);
res.z = res.w(:,4:5);
res.id    = res.x(:,1);
res.iq    = res.x(:,2);
res.omega = res.x(:,3);

% Evaluate energies and control inputs
N = length(res.t);
res.H    = zeros(N,1);
res.Hd   = zeros(N,1);
res.Hia  = zeros(N,1);
res.Hcl  = zeros(N,1);
res.u    = zeros(N,2);
for k = 1:N
    xk = res.x(k,:).';
    zk = res.z(k,:).';
    res.H(k)    = plant.H(xk);
    res.Hd(k)   = ctrl.Hd(xk);
    res.Hia(k)  = int.Hia(xk,zk);
    res.Hcl(k)  = res.Hd(k) + res.Hia(k);
    res.u(k,:) = (ctrl.u(xk) + int.v(xk,zk)).';
end

%%Plot states
figure
subplot(3,1,1)
plot(res.t,res.id,'LineWidth',1.5)
yline(ctrl.id_star,'--')
ylabel('$i_d$ [A]','Interpreter','latex')
grid on
subplot(3,1,2)
plot(res.t,res.iq,'LineWidth',1.5)
yline(ctrl.iq_star,'--')
ylabel('$i_q$ [A]','Interpreter','latex')
grid on
subplot(3,1,3)
plot(res.t,res.omega,'LineWidth',1.5)
yline(ctrl.omega_star,'--')
ylabel('$\omega$ [rad/s]','Interpreter','latex')
xlabel('Time [s]')
grid on

% Plot control voltages
figure
subplot(2,1,1)
plot(res.t,res.u(:,1),'LineWidth',1.5)
ylabel('$v_d$ [V]','Interpreter','latex')
grid on
subplot(2,1,2)
plot(res.t,res.u(:,2),'LineWidth',1.5)
ylabel('$v_q$ [V]','Interpreter','latex')
xlabel('Time [s]')
grid on

% Plot desired closed-loop Hamiltonian
figure
plot(res.t,res.Hcl,'LineWidth',1.5)
xlabel('Time [s]')
ylabel('$H_{cl}$','Interpreter','latex')
grid on