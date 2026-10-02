clear
clc
close all

% Constant ocean-current velocity expressed in the inertial frame.
% The primary IDA-PBC controller is designed assuming Vcn = 0.
Vcn = 0*[0.5; 0.5; 0];

% Desired vessel configuration
ctrl.etaStar = [1; 1; 0];

% Simulation length
Tsim = 60;

%% Define surface vessel model
% Model based on:
% "Manoeuvring Control of Underactuated Surface Vessels using Manifold
% Regulation for Port-Hamiltonian Systems", Renton et al.
sys.n = 3;

% Rotation from body-fixed to inertial coordinates
sys.J = @(eta) [cos(eta(3)), -sin(eta(3)), 0;
                sin(eta(3)),  cos(eta(3)), 0;
                0,            0,           1];

% Body-fixed inertia matrix
sys.Mb = diag([2; 1; 1]);

% Body-fixed Coriolis matrix
sys.Cb = @(pb) [0,     0,    -pb(2);
                0,     0,     pb(1);
                pb(2), -pb(1), 0];

% Body-fixed damping matrix
sys.Db = diag([1; 1; 1]);

% Relative body-fixed velocity
sys.nu = @(pb) sys.Mb\pb;

% Open-loop energy
sys.H = @(pb) 0.5*pb.'*(sys.Mb\pb);

% Open-loop vessel dynamics
% State coordinates used by the physical model are x_phys = [pb; eta],
% where pb is the body-fixed momentum.
%   dot(pb)  = -Cb(pb)*nu - Db(nu)*nu + tau
%   dot(eta) = J(eta)*nu + Vcn
sys.ode = @(pb,eta,tau) [-sys.Cb(pb)-sys.Db, -sys.J(eta).';
                        sys.J(eta), zeros(3)]*[sys.Mb\pb; zeros(3,1)] ...
                        + [zeros(3,1); Vcn] + [tau; zeros(3,1)];

%% Define transformed momentum coordinate representation
% Symbolic variables
syms eta1 eta2 eta3 real
syms pp1 pp2 pp3 real
eta_sym = [eta1; eta2; eta3];
p_sym   = [pp1; pp2; pp3];

% Transformed momentum  p = J(eta)^(-T) pb
sys.p = @(pb,eta) sys.J(eta)*pb;
sys.pb = @(p,eta) sys.J(eta).'*p;

% Transformed inverse inertia
% M^{-1}(eta) = J(eta) Mb^{-1} J(eta)'
ctrl.Mi = @(eta) sys.J(eta)*(sys.Mb\sys.J(eta).');

% Open-loop energy
sys.H_p = @(p, eta) 0.5*p.'*ctrl.Mi(eta)*p;

% Gradient of H in transformed momentum coordinates
sys.dH_p = matlabFunction(jacobian(sys.H_p(p_sym,eta_sym),[p_sym; eta_sym]).','vars',[{p_sym}, {eta_sym}]);

% C, D matrices in transformed momentum coordinates
sys.Dbar = @(eta) (sys.J(eta).')\sys.Db/sys.J(eta);
sys.dpb_deta = matlabFunction(jacobian((sys.J(eta_sym).')*p_sym,eta_sym),'vars',[{p_sym}, {eta_sym}]);
sys.Cbar = @(p,eta) -((sys.J(eta).')\sys.Cb(sys.pb(p,eta)))/sys.J(eta) + sys.dpb_deta(p,eta).'/sys.J(eta) - sys.J(eta).'\sys.dpb_deta(p,eta);

% Open-loop system represented in the transformed p coordinates
sys.ode_p = @(p,eta,tau) [sys.Cbar(p,eta)-sys.Dbar(eta), -eye(3);
                        eye(3), zeros(3)]*sys.dH_p(p,eta) ...
                        + [zeros(3,1); Vcn] + [-sys.J(eta).'; zeros(3)]*tau;
% % Open loop system dynamics transformed back to the pb coordinates. 
% % Only used for model validation
% sys.dp_deta = matlabFunction(jacobian((sys.J(eta_sym).')\p_sym,eta_sym),'vars',[{p_sym}, {eta_sym}]);
% sys.ode_pb = @(pb,eta,tau) [sys.J(eta).', sys.dpb_deta(pb,eta); zeros(3), eye(3)]*sys.ode_p(sys.p(pb,eta),eta,tau);

%% Primary IDA-PBC controller
% Primary control tuning gains
ctrl.Kp = eye(3);
ctrl.Kd = eye(3);

% Primary IDA-PBC input:
% tau = J'[-Kp(eta-etaStar) - Kd M^{-1}p + v]
% Since M^{-1}p = J Mb^{-1}pb, this becomes
% tau = J'[-Kp(eta-etaStar) - Kd J Mb^{-1}pb + v].
ctrl.tau = @(p,eta,v) sys.J(eta).'*(-ctrl.Kp*(eta-ctrl.etaStar) - ctrl.Kd*(ctrl.Mi(eta)*p) + v);

% Shaped Hamiltonian from primary control
ctrl.Hd = @(p,eta) 0.5*p.'*ctrl.Mi(eta)*p + 0.5*(eta-ctrl.etaStar).'*ctrl.Kp*(eta-ctrl.etaStar);

% Gradient of Hd in transformed momentum coordinates
ctrl.dHd = matlabFunction(jacobian(ctrl.Hd(p_sym,eta_sym),[p_sym; eta_sym]).','vars',[{p_sym}, {eta_sym}]);

% Input matrix for secondary control design
ctrl.g = [eye(3); zeros(3)];
ctrl.gPerp = [zeros(3), eye(3)];

% Interconnection and damping structure after primary control design
ctrl.Fd = @(p,eta) [ ...
    sys.Cbar(p,eta)-sys.Dbar(eta)-ctrl.Kd,   -eye(3);
    eye(3),                                     zeros(3)
    ];

%% Integral-action controller
% Tuning gain
int.Ki = 0.1*eye(1);

% Compute the matrix Theta
% Since
%   M^{-1} = J Mb^{-1} J'
% and J is orthogonal,
%   M^{-1/2} = J Mb^{-1/2} J'.
int.MbInvSqrt = diag(1./sqrt(diag(sys.Mb)));
int.Theta = @(eta) sys.J(eta)*int.MbInvSqrt*sys.J(eta).';

% Corresponding M^{1/2}
int.MbSqrt = diag(sqrt(diag(sys.Mb)));
int.ThetaInv = @(eta) sys.J(eta)*int.MbSqrt*sys.J(eta).';

% gamma(p,eta) = M^{-1/2}(eta)p
int.gamma = @(p,eta) int.Theta(eta)*p;
% Derivative of gamma with respect to eta
int.Gamma = matlabFunction(jacobian(int.gamma(p_sym,eta_sym),eta_sym),'vars',[{p_sym}, {eta_sym}]);
% Derivative of gamma with respect to full state
int.dgamma_dx = matlabFunction(jacobian(int.gamma(p_sym,eta_sym),[p_sym; eta_sym]),'vars',[{p_sym}, {eta_sym}]);
% ctrl.dgamma_dx = @(pb,eta) [ctrl.Theta(eta), ctrl.Gamma(pb,eta)];

% Prescribed integral output
% For the quadratic-Hamiltonian result in the paper,
%   L(eta) = -M^{-1/2}(eta)
% and hence
%   yo = -M^{-1/2}(eta) Kp (eta-etaStar).
% Therefore yo = 0 if and only if eta = etaStar.
int.yo = @(p,eta) -int.Theta(eta)*ctrl.Kp*(eta-ctrl.etaStar);

% % Integral-controller matrices
% For the surface-vessel example:
%   Lambda1 = 0
%   Lambda2 = 0
%   F12,1 = -M^{1/2} Gamma M^{-1/2}
%   F12,2 =  M^{-1/2}
%   F22   = 0
int.Lambda1 = zeros(3);
int.Lambda2 = zeros(3);
int.F12_1 = @(p,eta) -int.ThetaInv(eta)*int.Gamma(p,eta)*int.Theta(eta);
int.F12_2 = @(eta) int.Theta(eta);
int.F22   = zeros(3);
% Complete Lambda matrix
% Lambda = [Lambda1 Lambda2]*[g';gPerp]
int.Lambda = [int.Lambda1, int.Lambda2]*[ctrl.g.'; ctrl.gPerp];

% Moore-Penrose left inverse of g
int.gDag = (ctrl.g.'*ctrl.g)\ctrl.g.';

% Integral control law
int.v = @(p,eta,z) (-int.F12_1(p,eta) + int.gDag*ctrl.Fd(p,eta)*int.dgamma_dx(p,eta).')*int.Ki*(int.gamma(p,eta)-z) ...
                - int.Lambda*(ctrl.dHd(p,eta) + int.dgamma_dx(p,eta).'*int.Ki*(int.gamma(p,eta)-z));

% Closed-loop augmented energy
int.Hia = @(p,eta,z) 0.5*(int.gamma(p,eta)-z).'*int.Ki*(int.gamma(p,eta)-z);
int.Hcl = @(p,eta,z) ctrl.Hd(p,eta) + int.Hia(p,eta,z);

%% Run simulation
% Initial conditions
sim.eta0 = [1; 2; 0];
sim.pb0  = [-1; -1; 1];
sim.z0 = zeros(sys.n,1);
sim.x0 = [sim.pb0;
          sim.eta0;
          sim.z0];

% Closed-loop augmented dynamics:
v = zeros(3,1);
ode = @(t,x) [sys.ode(x(1:3),x(4:6),ctrl.tau(sys.p(x(1:3),x(4:6)),x(4:6),int.v(sys.p(x(1:3),x(4:6)),x(4:6),x(7:9))));
                int.yo(sys.p(x(1:3),x(4:6)),x(4:6))];

% Simulate system
options = odeset('RelTol',1e-6);
[res.t,res.x] = ode45(ode,[0 Tsim],sim.x0,options);

%% Plot results
% Unpack states
res.pb  = res.x(:,1:3);
res.eta = res.x(:,4:6);
res.z   = res.x(:,7:9);
N = length(res.t);

% Preallocate arrays
res.p     = zeros(N,3);
res.H   = zeros(N,1);
res.Hd  = zeros(N,1);
% Compute values of interest
for k = 1:N
    pb  = res.pb(k,:).';
    eta = res.eta(k,:).';
    z   = res.z(k,:).';

    res.p(k,:) = sys.p(pb,eta).';

    res.H(k)   = sys.H(pb);
    res.Hd(k)  = ctrl.Hd(res.p(k,:).',eta);
    res.Hcl(k) = int.Hcl(res.p(k,:).',eta,z);
end

% Plot vessel configuration
figure(1)
subplot(3,1,1)
plot(res.t,res.eta(:,1),'LineWidth',1.2)
hold on
yline(ctrl.etaStar(1),'--')
grid on
ylabel('\eta_1')
subplot(3,1,2)
plot(res.t,res.eta(:,2),'LineWidth',1.2)
hold on
yline(ctrl.etaStar(2),'--')
grid on
ylabel('\eta_2')
subplot(3,1,3)
plot(res.t,res.eta(:,3),'LineWidth',1.2)
hold on
yline(ctrl.etaStar(3),'--')
grid on
ylabel('\eta_3')
xlabel('Time (s)')

% Plot integral state
figure
plot(res.t,res.z,'LineWidth',1.2)
grid on
xlabel('Time (s)')
ylabel('z')
legend('z_1','z_2','z_3','Location','best')

% Plot transformed momentum
figure
plot(res.t,res.p,'LineWidth',1.2)
grid on
xlabel('Time (s)')
ylabel('p')
legend('p_1','p_2','p_3','Location','best')

% Plot energies
figure
subplot(3,1,1)
plot(res.t,res.H,'LineWidth',1.2)
grid on
ylabel('H')
subplot(3,1,2)
plot(res.t,res.Hd,'LineWidth',1.2)
grid on
ylabel('H_d')
subplot(3,1,3)
plot(res.t,res.Hcl,'LineWidth',1.2)
grid on
xlabel('Time (s)')
ylabel('H_{cl}')