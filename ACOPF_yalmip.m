%% AC Optimal Power Flow

% The branch flow limit is not considered.

% yalmip modeling, ipopt solving

% Dependency libraries: yalmip, ipopt, matpower

% Author: Qingju Luo
% Email: luoqingju@qq.com
% School of Electric Power Engineering, South China University of Technology
% Integrated Smart Energy System Optimal Operation and Control, ISESOOC

% The author's ability is limited, there will inevitably be errors and inadequacies, please criticize and correct!

% For learning and communication only!!!
% For learning and communication only!!!
% For learning and communication only!!!

%% YALMIP
% !!!!!!!!!! The yalmip (version: 20230622) call to ipopt can only use a quasi-Newton approximation to the Hessian !!!!!!!!!!
% !!!!!!!!!! The number of iterations is more than in the case of the exact Hessian !!!!!!!!!!
% !!!!!!!!!! When the system is large, it may not be solvable !!!!!!!!!!

clc
clear

define_constants;

% No flow limits
mpc = case14;
% mpc = case57;
% mpc = case_ieee30;
% mpc = case118;
% mpc = case300;

% With flow limits
% mpc = case9;
% mpc = case30;
% mpc = case2383wp;

%% Get Data
mpc = ext2int(mpc);

baseMVA = mpc.baseMVA;
bus = mpc.bus;
gen = mpc.gen;
branch = mpc.branch;
gencost = mpc.gencost;

nb = size(bus, 1); %% number of buses
nl = size(branch, 1); %% number of branches
ng = size(gen, 1); %% number of dispatchable injections

ref = find(bus(:, BUS_TYPE) == REF); % reference bus
Va_ref = bus(ref, VA) * (pi / 180); % reference bus voltage angle (radians)

Vmu = bus(:, VMAX); % maximum voltage magnitude (p.u.)
Vml = bus(:, VMIN); % minimum voltage magnitude (p.u.)

Pmax = gen(:, PMAX) / baseMVA; % maximum real power output (p.u.)
Pmin = gen(:, PMIN) / baseMVA; % minimum real power output (p.u.)

Qmax = gen(:, QMAX) / baseMVA; % maximum reactive power output (p.u.)
Qmin = gen(:, QMIN) / baseMVA; % minimum reactive power output (p.u.)

% find/prepare polynomial generator costs
Qpg = gencost(:, COST) * baseMVA^2; % quadratic
cpg = gencost(:, COST+1) * baseMVA; % linear
kpg = gencost(:, COST+2); % constant


Ybus = makeYbus(baseMVA, bus, branch); % admittance matrix

Cg = sparse(gen(:, GEN_BUS), (1:ng)', 1, nb, ng); %% connection matrix for generators & buses

Sd = (bus(:, PD) + 1j * bus(:, QD)) / baseMVA; % complex power demand (p.u.)

%% Define decision variables
disp('Start modeling...')

Va = sdpvar(nb, 1); % voltage angle (radians)
Vm = sdpvar(nb, 1); % voltage magnitude (p.u.)
Pg = sdpvar(ng, 1); % real power output (p.u.)
Qg = sdpvar(ng, 1); % reactive power output (p.u.)

%% Define constraints
con = [];
con = con + [Va(ref) == Va_ref]; %#ok<*NBRAK1>
con = con + [Vml <= Vm <= Vmu]; %#ok<*CHAIN>
con = con + [Pmin <= Pg <= Pmax];
con = con + [Qmin <= Qg <= Qmax];

V = Vm .* exp(1j*Va); % complex voltage
con = con + [V .* conj(Ybus * V) == Cg * (Pg + 1j * Qg) - Sd]; % complex power balance

%% Define the objective function
obj = Pg' * (Qpg .* Pg) + cpg' * Pg + sum(kpg); % objective function

%% Define initial value
assign(Va, Va_ref);
assign(Vm, (Vml + Vmu)/2);
assign(Pg, (Pmin + Pmax)/2);
assign(Qg, (Qmin + Qmax)/2);

% assign(Va, mpc.bus(:, VA) / 180 * pi);
% assign(Vm, mpc.bus(:, VM));
% assign(Pg, mpc.gen(:, PG) / baseMVA);
% assign(Qg, mpc.gen(:, QG) / baseMVA);

%% Solve
disp('Start solving...')
ops_yalmip = sdpsettings('solver', 'ipopt', 'verbose', 3, 'usex0', 1, 'savesolveroutput', 1);

ops_yalmip.ipopt.hessian_approximation
% hessian_approximation: Indicates what Hessian information is to be used.
% This determines which kind of information for the Hessian of the Lagrangian function is used by the algorithm.
% The default value for this string option is "exact".
% exact: Use second derivatives provided by the NLP.
% limited-memory: Perform a limited-memory quasi-Newton approximation


% ops_yalmip.ipopt.hessian_approximation = 'exact';
% *** Error using Ipopt Matlab interface: ***
% You must supply a callback function for computing the Hessian unless you decide to use a quasi-Newton approximation to the Hessian.



% !!!!!!!!!! The yalmip (version: 20230622) call to ipopt can only use a quasi-Newton approximation to the Hessian !!!!!!!!!!
% !!!!!!!!!! The number of iterations is more than in the case of the exact Hessian !!!!!!!!!!
% !!!!!!!!!! When the system is large, it may not be solvable !!!!!!!!!!



sol = solvesdp(con, obj, ops_yalmip);

if sol.problem ~= 0
    disp(sol.info)
    return
end

Va = value(Va);
Vm = value(Vm);
Pg = value(Pg);
Qg = value(Qg);

x_yalmip = [Va; Vm; Pg; Qg];
obj_yalmip = value(obj);

%% Compare with MATPOWER

mpc.branch(:, RATE_A) = 0; % The branch flow limit is not considered.

opt = ipopt_options();
opt.hessian_approximation = ops_yalmip.ipopt.hessian_approximation;
% opt.hessian_approximation = 'exact';
opt.max_iter = ops_yalmip.ipopt.max_iter;
mpopt = mpoption('opf.ac.solver', 'IPOPT', 'verbose', 2, 'out.all', 0, 'ipopt.opts', opt);
res = runopf(mpc, mpopt);
if res.success == 0
    disp('MATPOWER runopf failed !!!')
    return
end

fprintf(['\nThe error of the objective function is:', num2str(norm(obj_yalmip - res.f)), '\n'])
fprintf(['\nThe error of the decision variables is:', num2str(norm(x_yalmip - res.x)), '\n'])
