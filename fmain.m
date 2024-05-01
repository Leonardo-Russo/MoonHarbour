function fmain(workspace_path, verbose, misalignment_type)


% Define Options
global opt
opt = struct('name', "Options");
opt.saveplots = false;
opt.create_animation = false;
opt.show_progress = false;
opt.compute_target = false;
opt.additional_plots = false;
opt.showgui = false;
opt.N = 1000;                   % n° of points for the Interpolation
opt.RelTolODE = 1e-7;           % options for ode()
opt.AbsTolODE = 1e-6;

OptionsODE = odeset('RelTol', opt.RelTolODE, 'AbsTol', opt.AbsTolODE);

% Define Thrust Misalignment
misalignments = [];
misalignment = define_misalignment_error(misalignment_type);
misalignments = [misalignments; misalignment];

tic
%% Hyperparameters and Settings

% Define Global Variables
global DU TU MU Rm muM pbar log

% Define Constants
DU = 1738;                                              % km
TU = sqrt(1 / 4902.7779 * DU^3);                        % s
Rm = 1738.1 / DU;                                       % DU
muM = 4902.7779 * TU^2 / DU^3;                          % DU^3/TU^2
muE = 398600.4418 * TU^2 / DU^3;                        % DU^3/TU^2
muS = 132712440018 * TU^2 / DU^3;                       % DU^3/TU^2
deltaE = deg2rad(23.45);                                % rad
psiM = deg2rad(-81.7 + 360/18.6 * (5 + 4.5/12));        % rad
deltaM = deg2rad(1.5);                                  % rad
Dsol = 86400;                                           % s
sec2hrs = 1/3600;                                       % hrs
MU = 1/DU^2;

% Define Time Domain
date0 = datetime('2025-05-23 9:56:10');
datef = datetime('2025-05-23 21:56:10');    % @ periselenium
total_time = datef - date0;

% Define Control Limit
g0 = 9.80665;
u_limit = 5e-5*g0;

% Define the Chaser Initial Conditions ~ norm(rhodot_LVLH) = 1 m/s -> condition applied for the finetuning of the gain parameters
RHO0_LVLH = [1.5, 0, 0, -1e-6, -1e-3, -1e-3]';              % km, km/s
RHO0_LVLH = [RHO0_LVLH(1:3)/DU; RHO0_LVLH(4:6)/DU*TU];      % adim

% Define Desired Conditions for Docking
RHOf_LVLH = [-5e-3, 0, 0, 5e-5, 0, 0]';                     % km, km/s
RHOf_LVLH = [RHOf_LVLH(1:3)/DU; RHOf_LVLH(4:6)/DU*TU];      % adim

% Define Direct Approach Conditions
rhof_LVLH_DA = 5e-3/DU * RHO0_LVLH(1:3)/norm(RHO0_LVLH(1:3));
rhodotf_LVLH_DA = -5e-5/DU*TU * RHO0_LVLH(1:3)/norm(RHO0_LVLH(1:3));
RHOf_LVLH_DA = [rhof_LVLH_DA; rhodotf_LVLH_DA];

%% Propagate Reference Target Trajectory

if opt.compute_target

    % Interpolate the Ephemeris and Retrieve Target's Initial State
    [X0t_MCI, COE0t, MEE0t, EarthPPsMCI, DSGPPsMCI, SunPPsMCI, MoonPPsECI, time, t0, tf, opt.N] = ...
        EphemerisHandlerExp(deltaE, psiM, deltaM, opt.N, date0, datef);
    
    % Define the timespan for the propagation
    tspan_ref = linspace(t0, tf, opt.N);
    dt_ref = tspan_ref(2) - tspan_ref(1);
    
    % Perform the Propagation of the Target Trajectory
    if opt.show_progress
        pbar = waitbar(0, 'Performing the Target Trajectory Propagation');
    end
    [~, MEEt_ref] = ode113(@(t, MEE) DynamicalModelMEE(t, MEE, EarthPPsMCI, SunPPsMCI, ...
        muE, muS, MoonPPsECI, deltaE, psiM, deltaM, t0, tf), tspan_ref, MEE0t, OptionsODE);
    if opt.show_progress
        close(pbar);
    end
    
    % Conversion from MEE to COE
    COEt_ref = MEE2COE(MEEt_ref);
    
    % Conversion from COE to MCI
    Xt_MCI_ref = COE2rvPCI(COEt_ref, muM);
    
    % Interpolate the Angular Velocity of LVLH wrt MCI
    [~, omegadotPPsLVLH] = TargetHandler(Xt_MCI_ref, COEt_ref, MEEt_ref, tspan_ref, ...
        EarthPPsMCI, SunPPsMCI, MoonPPsECI, deltaE, psiM, deltaM, muE, muS);

    % save('Data/temp/Target Propagation.mat');

else

    load('Data/temp/Target Propagation.mat', 'X0t_MCI', 'COE0t', 'MEE0t', 'EarthPPsMCI', 'DSGPPsMCI', 'SunPPsMCI', 'MoonPPsECI', 'time', 't0', 'tf', 'tspan_ref', 'dt_ref', 'MEEt_ref', 'COEt_ref', 'Xt_MCI_ref', 'omegadotPPsLVLH');

end

%% Direct Approach: Chaser Backwards Propagation from Final Conditions

% Combine the Target and Chaser States into TC ~ [Target; Chaser] State
TC0 = [MEE0t; RHO0_LVLH];

% Define Natural Drift Propagation Settings
optODE_back = odeset('RelTol', 1e-9, 'AbsTol', 1e-9, 'Events', @(t, Y) driftstop_back(t, Y));

% Define the Backwards Drift tspan
tspan_back_DA = linspace(tf, t0 + (tf-t0)/2, opt.N);

% Define Final Conditions
MEEft = MEEt_ref(end, :)';
TCf_backdrift_DA = [MEEft; RHOf_LVLH_DA];

% Perform the Backpropagation of Target and Chaser Trajectories
[tspan_back_DA, TC_backdrift_DA, t_bwdstop_DA, ~, ~] = ode113(@(t, TC) NaturalRelativeMotion(t, TC, EarthPPsMCI, SunPPsMCI, ...
    muE, muS, MoonPPsECI, deltaE, psiM, deltaM, t0, tf, omegadotPPsLVLH, 0), tspan_back_DA, TCf_backdrift_DA, optODE_back);

% Check Successful Back Drift
if isempty(t_bwdstop_DA)
    warning('Could not complete Backwards Natural Drift from Direct Approach Conditions.');
end

% Retrieve the States and epoch for Final Drift
TC0_backdrift_DA = TC_backdrift_DA(end, :)';
t0_backdrift_DA = tspan_back_DA(end);


%% Direct Approach: Generate the Reference Chaser Trajectory for Hybrid Predictive Control

% Set the Via Points and Interpolate the Reference Trajectory
[RHOdPPsLVLH_DA, viapoints_DA, t_viapoints_DA] = ReferenceTrajectory(TC0, TC0_backdrift_DA, t0, t0_backdrift_DA, [0, 1]);

% Control Settings
prediction_interval_DA = 120;                  % s - one prediction every 2 minutes
prediction_dt_DA = prediction_interval_DA/TU;     % adim - prediction interval adimensionalized

prediction_maxlength_DA = fix(seconds(total_time)/prediction_interval_DA);        % max length of the predictive controls
check_times_DA = zeros(prediction_maxlength_DA, 1);           % gets all the times at which one should check in continue to saturate by future prediction

for k = 1 : prediction_maxlength_DA
    check_times_DA(k) = t0 + (k-1) * prediction_dt_DA;        % compute the check_times vector
end

prediction_delta_DA = (30*60)/TU;      % predictive propagation of 30 min forward
N_inner_DA = round(prediction_delta_DA/dt_ref) - 1;       % n° of points in the inner propagation


%% Direct Approach: Perform Chaser Rendezvous Manoeuvre

% Define Initial Conditions
x70 = 1;                % initial mass ratio
TCC0 = [TC0; x70];      % initial TargetControlledChaser State
event_odefun = 1;       % 1 means it'll stop at saturation

% Perform Hybrid-Predictive Feedback Control
if opt.show_progress
    pbar = waitbar(0, 'Performing the Direct Approach Rendezvous');
end
[tspan_ctrl_DA, TCC_ctrl_DA] = odeHamHPC(@(t, TCC) HybridPredictiveControl(t, TCC, EarthPPsMCI, SunPPsMCI, muM, ...
    muE, muS, MoonPPsECI, deltaE, psiM, deltaM, omegadotPPsLVLH, t0, tf, RHOdPPsLVLH_DA, DU, TU, check_times_DA, prediction_delta_DA, N_inner_DA, misalignment),...
    [t0, t0_backdrift_DA], TCC0, opt.N, @is_terminal_distance, event_odefun);

% Check if Terminal Conditions are reached
[~, terminal_flag] = is_terminal_distance(tspan_ctrl_DA(end), TCC_ctrl_DA(end, :));
if terminal_flag
    error('Reached Conditions under Saturation!')
end

% Define stopSaturationTime
stopSaturationTime = tspan_ctrl_DA(end);

% Retrive Final Kp
[~, ~, ~, ~, ~, ~, ~, ~, ~, ~, kp, ~] = HybridPredictiveControl(tspan_ctrl_DA(end), TCC_ctrl_DA(end, :), EarthPPsMCI, SunPPsMCI, muM, ...
    muE, muS, MoonPPsECI, deltaE, psiM, deltaM, omegadotPPsLVLH, t0, tf, RHOdPPsLVLH_DA, DU, TU, check_times_DA, prediction_delta_DA, N_inner_DA, misalignment);

% Perform Natural Feedback Control
[tspan_temp, TCC_temp] = odeHamHPC(@(t, TCC) NaturalFeedbackControl(t, TCC, EarthPPsMCI, SunPPsMCI, muM, ...
    muE, muS, MoonPPsECI, deltaE, psiM, deltaM, omegadotPPsLVLH, t0, tf, RHOdPPsLVLH_DA, kp, DU, TU, misalignment, opt.show_progress, 0), ...
    [tspan_ctrl_DA(end), t0_backdrift_DA], TCC_ctrl_DA(end, :), opt.N, @is_terminal_distance);

% Retrieve Final State Values
tspan_ctrl_DA = [tspan_ctrl_DA; tspan_temp];
TCC_ctrl_DA = [TCC_ctrl_DA; TCC_temp];
TCC_T0 = TCC_ctrl_DA(end, :);           % this will be the initial point for the Terminal Trajectory
t0_T = tspan_ctrl_DA(end);              % this will be the initial time for the Terminal Trajectory

if verbose
    % fprintf('Reached 50m Distance at: %02d hrs, %02d mins, %2.3f secs\n', time_elapsed(t0_T, t0, TU));
    disp(sprintf('Reached 100m Distance @ <strong>[%02d:%02d:%06.3f]</strong> (hh:mm:ss)', time_elapsed(t0_T, t0, TU)));
end

% % Show Original Trajectory
% figure('name', 'Original Trajectory');
% DrawTrajLVLH3D(TCC_ctrl_DA(:, 7:9)*DU);
% testPPs(RHOdPPsLVLH_DA(1:3), tspan_ctrl_DA);


%% Terminal Trajectory: Chaser Backwards Propagation

% Define Natural Drift Propagation Settings
optODE_back = odeset('RelTol', 1e-9, 'AbsTol', 1e-9, 'Events', @(t, Y) driftstop_back(t, Y));

% Define the Backwards Drift tspan
tspan_back = linspace(tf, t0 + (tf-t0)/2, opt.N);

% Define Final Conditions
MEEft = MEEt_ref(end, :)';
TCf_backdrift = [MEEft; RHOf_LVLH];

% Perform the Backpropagation of Target and Chaser Trajectories
[tspan_back, TC_backdrift, t_bwdstop, ~, ~] = ode113(@(t, TC) NaturalRelativeMotion(t, TC, EarthPPsMCI, SunPPsMCI, ...
    muE, muS, MoonPPsECI, deltaE, psiM, deltaM, t0, tf, omegadotPPsLVLH, 0), tspan_back, TCf_backdrift, optODE_back);

% Check Successful Back Drift
if isempty(t_bwdstop)
    warning('Could not complete Backwards Natural Drift with the given Initial Conditions.');
end

% Retrieve the States and epoch for Final Drift
TC0_backdrift = TC_backdrift(end, :)';
t0_backdrift = tspan_back(end);


%% Terminal Trajectory: Chaser Complete Docking Manoeuvre

% Initialize Regenerative Quantities
tspan_ctrl = [];
Y_ctrl = [];
RHOdPPsLVLH_T = [];
indices_ctrl = [1];
Q_N2C_AOCS_stack = [];
Qe_AOCS_stack = [];
u_rt_AOCS_stack = [];
xc_plot_stack = [];
yc_plot_stack = [];
zc_plot_stack = [];
u_AOCS_stack = [];
omega_cPPs_rt_stack = [];
omegadot_cPPs_rt_stack = [];
Q_N2C_PPs_rt_stack = [];
sign_qe0_0_stack = [];
omega_c_rt_stack = [];
omegadot_c_rt_stack = [];
Tc_AOCS_stack = [];
xb_AOCS_stack = [];
u_rt_norm_stack = [];
u_AOCS_norm_stack = [];
TCC_PPs_stack = [];


% Define Propagation Settings
dt_regen = 1*60/TU;             % renegerative propagation interval
dt_min = 0.5*60/TU;               % minimum propagation interval
prop_step = 1/TU;               % propagation time step
max_branches = 500;             % maximum n° of branches of the regenerative trajectory
optODE_rt = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);
optODE_AOCS = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);

% Define Global signvect Variable
global signvect

% Define Approach Options
ref_axis_tol = 1e-2;                % to make sure that initially xc_MCI is not aligned with c3_MCI
ref_stored = zeros(3, 1);           % stores the reference axes values at last time of previous branch

% Define Initial Conditions
TCC_rt0 = TCC_T0';
t0_rt = t0_T;

w_0 = zeros(3, 1);                      % rad/s
omegas_0 = zeros(4, 1);                 % rad/s
qb_0 = [0.1, 0.3, -0.5]';
qb0_0 = sqrt(1 - norm(qb_0)^2);
Xb0 = [qb0_0; qb_0];                    % body attitude wrt MCI

Y0_rt = [TCC_rt0; Xb0; w_0; omegas_0];     % AOCS extended State


% temporary stuff
debug = 0;
opt.initially_aligned = true;
omega_n = 0.05*TU;      % rad/s



for branch = 1 : max_branches


    % ----- Propagate Trajectory w/o Attitude ----- %

    TCC_rt0 = Y0_rt(1:13);

    misalignment = define_misalignment_error(misalignment_type);
    misalignments = [misalignments; misalignment];

    % Set the Via Points and Interpolate the Reference Terminal Trajectory
    if branch == 1
        [RHOdPPsLVLH_rt, viapoints_rt, t_viapoints_rt, rho1_0, t1_0, l_hat_0] = ReferenceTrajectory(TCC_rt0(1:12), TC0_backdrift, t0_rt, t0_backdrift, [1, 1]);
    else
        [RHOdPPsLVLH_rt, viapoints_rt, t_viapoints_rt] = ReferenceTrajectory(TCC_rt0(1:12), TC0_backdrift, t0_rt, t0_backdrift, [1, 1], NaN, t1_0, NaN);
    end

    % Set Final Propagation Time and Define timespan
    if length(t_viapoints_rt) >= 3
        if t0_rt + dt_min < t_viapoints_rt(2)
            tf_rt = t_viapoints_rt(2);
        else
            tf_rt = t0_backdrift;
        end
    else
        tf_rt = t0_backdrift;
    end

    % Define timespan
    tspan_rt = [t0_rt : prop_step : tf_rt]';

    % Propagate to Final Propagation Time
    [~, TCC_rt] = ode113(@(t, TCC) NaturalFeedbackControl(t, TCC, EarthPPsMCI, SunPPsMCI, muM, ...
        muE, muS, MoonPPsECI, deltaE, psiM, deltaM, omegadotPPsLVLH, t0, tf, RHOdPPsLVLH_rt, kp, DU, TU, misalignment, 0, 1), ...
        tspan_rt, TCC_rt0, optODE_rt); 

    % Interpolate Trajectory
    TCC_PPs = get_statePP(tspan_rt, TCC_rt);


    % ----- Retrieve Commanded Attitude ----- %

    % Commanded Attitude Frame will be defined by [xc, yc, zc] in MCI
    M_rt = length(tspan_rt);
    xc_MCI = zeros(M_rt, 3);
    yc_MCI = zeros(M_rt, 3);
    zc_MCI = zeros(M_rt, 3);
    R_N2C = zeros(3, 3, M_rt);
    Q_N2C = zeros(M_rt, 4);
    u_rt = zeros(M_rt, 3);
    u_rt_norm = zeros(M_rt, 1);
    
    MEEt_rt = TCC_rt(:, 1:6);       % retrieve Target State
    COEt_rt = MEE2COE(MEEt_rt);
    Xt_MCI_rt = COE2rvPCI(COEt_rt, muM);
    
    c3_MCI = [0, 0, 1]';
    
    for i = 1 : M_rt
    
        % Retrieve Thrust Acceleration in LVLH
        [~, ~, ~, ~, u_rt(i, :), ~, ~, ~, ~] = NaturalFeedbackControl(tspan_rt(i), TCC_rt(i, :), EarthPPsMCI, SunPPsMCI, muM, ...
            muE, muS, MoonPPsECI, deltaE, psiM, deltaM, omegadotPPsLVLH, t0, tf, RHOdPPsLVLH_rt, kp, DU, TU, misalignment, 0, 0);
        xc_LVLH = u_rt(i, :)' / norm(u_rt(i, :));        % normalize to find xc versor in LVLH
        u_rt_norm(i) = norm(u_rt(i, :));        % normalize to find xc versor in LVLH
    
        % Rotate from LVLH to MCI
        [R_LVLH2MCI, ~] = get_rotLVLH2MCI(Xt_MCI_rt(i, :)', tspan_rt(i), EarthPPsMCI, SunPPsMCI, MoonPPsECI, muE, muS, deltaE, psiM, deltaM);
        xc_MCI(i, :) = R_LVLH2MCI*xc_LVLH;

        % Choose Commanded Reference Frame approach
        if branch == 1 && i == 1
            if norm(cross(c3_MCI, xc_MCI(i, :)')) > ref_axis_tol
                ref3_MCI = c3_MCI;      % sets c3_MCI as the initial third axis reference
                ref_axis_flag = 1;
            else
                error('Second Commanded Attitude approach still needs to be defined.')
            end
        end
        if branch > 1 && i == 1
            if ref_axis_flag == 1
                ref3_MCI = ref_stored(:, 1);
            else
                error('Second Commanded Attitude approach still needs to be defined.')
            end

        end

        % Continue Computing Commanded Attitude
        yc_MCI(i, :) = cross(ref3_MCI, xc_MCI(i, :)')/norm(cross(ref3_MCI, xc_MCI(i, :)'));
        zc_MCI(i, :) = cross(xc_MCI(i, :)', yc_MCI(i, :)');
    
        R_N2C(:, :, i) = [xc_MCI(i, :)', yc_MCI(i, :)', zc_MCI(i, :)']';
        
        if i == 1
            if branch == 1
                signvect = [1 1 1 1]';
                [q0c, qc] = matrixToQuat(R_N2C(:, :, i));
            else
                signvect = sign(Q_N2C_AOCS_stack(end, :))';
                [q0c, qc] = matrixToQuat(R_N2C(:, :, i));
            end
        else
            [q0c, qc] = matrixToQuat(R_N2C(:, :, i)); 
        end

        Q_N2C(i, :) = [q0c, qc'];

    end

    % Interpolate R_N2C and compute its derivative
    R_N2C_PPs_rt = get_rotPPs(tspan_rt, R_N2C);
    Rdot_N2C_PPs_rt = fnder_rots(R_N2C_PPs_rt);

    % Interpolate Commanded Attitude Quaternions
    Q_N2C_PPs_rt = get_statePP(tspan_rt, Q_N2C);
    
    % Compute angular velocity of Commanded Attitude and its derivative
    omega_c_rt = zeros(M_rt, 3);
    for k = 1 : M_rt
        Rdot_N2C_rt = rotppsval(Rdot_N2C_PPs_rt, tspan_rt(k));
        omega_c_rt(k, :) = unskew(-Rdot_N2C_rt * R_N2C(:, :, k)');
    end
    omega_cPPs_rt = get_statePP(tspan_rt, omega_c_rt);
    omegadot_cPPs_rt = [fnder(omega_cPPs_rt(1), 1);
                      fnder(omega_cPPs_rt(2), 1);
                      fnder(omega_cPPs_rt(3), 1)];
    omegadot_c_rt = zeros(M_rt, 3);
    for k = 1 : M_rt
        omegadot_c_rt(k, :) = ppsval(omegadot_cPPs_rt, tspan_rt(k));
    end

    

    % ----- Propagate Trajectory AND Attitude ----- %
    
    % Define Time Domain
    t0_AOCS = t0_rt;
    tf_AOCS = min(t0_rt + dt_regen, tf_rt);
    % branch_squish_tol = 30/TU;  % 30 s
    branch_squish_tol = prop_step;
    if tf_rt - tf_AOCS < branch_squish_tol
        tf_AOCS = tf_rt;
    end
    tspan_AOCS = [t0_AOCS : prop_step : tf_AOCS]';
    
    qc0_0 = Q_N2C(1, 1);
    qc_0 = Q_N2C(1, 2:4)';
    
    % Align Body and Commanded at initial time
    if branch == 1 && opt.initially_aligned
        Y0_rt(14) = qc0_0;
        Y0_rt(15:17) = qc_0;
    end
    
    qb0_0 = Y0_rt(14);
    qb_0 = Y0_rt(15:17);
    sign_qe0_0 = sign(qc0_0*qb0_0 + qc_0'*qb_0);        % needed for short rotation
    
    % % Perform the Attitude Propagation - ode113
    % [tspan_AOCS, Y_AOCS] = ode113(@(t, Y) AOCS(t, Y, EarthPPsMCI, SunPPsMCI, muM, ...
    %     muE, muS, MoonPPsECI, deltaE, psiM, deltaM, omegadotPPsLVLH, t0, tf, RHOdPPsLVLH_rt, kp, omega_n, DU, TU, MU, TCC_PPs, omega_cPPs_rt, omegadot_cPPs_rt, Q_N2C_PPs_rt, sign_qe0_0, misalignment, opt.show_progress, 1), ...
    %     tspan_AOCS, Y0_rt, optODE_AOCS);
    
    % Perform the Attitude Propagation - odeHam
    [tspan_AOCS, Y_AOCS] = odeHamHPC(@(t, Y) AOCS(t, Y, EarthPPsMCI, SunPPsMCI, muM, ...
        muE, muS, MoonPPsECI, deltaE, psiM, deltaM, omegadotPPsLVLH, t0, tf, RHOdPPsLVLH_rt, kp, omega_n, DU, TU, MU, TCC_PPs, omega_cPPs_rt, omegadot_cPPs_rt, Q_N2C_PPs_rt, sign_qe0_0, misalignment, opt.show_progress, 0), ...
        [t0_AOCS, tf_AOCS], Y0_rt, length(tspan_AOCS)-1);
    
    % Retrieve Attitude Evolution
    Xb = Y_AOCS(:, 14:17);
    w = Y_AOCS(:, 18:20);
    omegas = Y_AOCS(:, 21:24);

    % In-Step Post Processing
    M_AOCS = length(tspan_AOCS);
    Q_N2C_AOCS = zeros(M_AOCS, 4);
    Qe_AOCS = zeros(M_AOCS, 4);
    u_AOCS = zeros(M_AOCS, 3);
    Tc_AOCS = zeros(M_AOCS, 3);
    xb_AOCS = zeros(M_AOCS, 3);
    u_AOCS_norm = zeros(M_AOCS, 1);

    for j = 1 : M_AOCS

        % AOCS Reconstruction
        [~, ~, ~, ~, u_AOCS(j, :), ~, ~, ~, ~, Tc_AOCS(j, :), xb_AOCS(j, :)] = AOCS(tspan_AOCS(j), Y_AOCS(j, :)', EarthPPsMCI, SunPPsMCI, muM, ...
        muE, muS, MoonPPsECI, deltaE, psiM, deltaM, omegadotPPsLVLH, t0, tf, RHOdPPsLVLH_rt, kp, omega_n, DU, TU, MU, TCC_PPs, omega_cPPs_rt, omegadot_cPPs_rt, Q_N2C_PPs_rt, sign_qe0_0, misalignment, 0, 1);

        % Attitude Reconstruction
        Q_N2C_AOCS(j, :) = ppsval(Q_N2C_PPs_rt, tspan_AOCS(j));
        qc0_AOCS = Q_N2C_AOCS(j, 1);
        qc_AOCS = Q_N2C_AOCS(j, 2:4)';
        qb0_AOCS = Xb(j, 1);
        qb_AOCS = Xb(j, 2:4)';
        Qe_AOCS(j, 1) = qc0_AOCS*qb0_AOCS + qc_AOCS' * qb_AOCS;
        Qe_AOCS(j, 2:4) = (-qc_AOCS*qb0_AOCS + qc0_AOCS*qb_AOCS - skew(qc_AOCS)*qb_AOCS)';

        u_AOCS_norm(j) = norm(u_AOCS(j, :));

    end

    % Store Commanded Attitude reference versor
    Q_N2C_final = Q_N2C_AOCS(end, :);
    R_N2C_final = q2C(Q_N2C_final(1), Q_N2C_final(2:4)');
    zc_ref = R_N2C_final(3, :)';
    if ref_axis_flag == 1
        ref_stored = zc_ref;    % store last reference axis value
    end

    % Additional Post-Processing
    % THIS ONLY WORKS SINCE THE TIMESTEP IS THE SAME
    u_rt_AOCS = u_rt(1:M_AOCS, :);
    xc_plot = squeeze(R_N2C(1, :, 1:M_AOCS))';
    yc_plot = squeeze(R_N2C(2, :, 1:M_AOCS))';
    zc_plot = squeeze(R_N2C(3, :, 1:M_AOCS))';


    % ----- Results Processing and Visualization ----- %    

    % Stack the Output
    tspan_ctrl = [tspan_ctrl; tspan_AOCS];
    Y_ctrl = [Y_ctrl; Y_AOCS];
    RHOdPPsLVLH_T = [RHOdPPsLVLH_T, RHOdPPsLVLH_rt];
    TCC_PPs_stack = [TCC_PPs_stack, TCC_PPs];
    indices_ctrl = [indices_ctrl; M_AOCS];
    Q_N2C_AOCS_stack = [Q_N2C_AOCS_stack; Q_N2C_AOCS];
    Qe_AOCS_stack = [Qe_AOCS_stack; Qe_AOCS];
    u_rt_AOCS_stack = [u_rt_AOCS_stack; u_rt_AOCS];
    xc_plot_stack = [xc_plot_stack; xc_plot];
    yc_plot_stack = [yc_plot_stack; yc_plot];
    zc_plot_stack = [zc_plot_stack; zc_plot];
    u_AOCS_stack = [u_AOCS_stack; u_AOCS];
    omega_cPPs_rt_stack = [omega_cPPs_rt_stack, omega_cPPs_rt];
    omegadot_cPPs_rt_stack = [omegadot_cPPs_rt_stack, omegadot_cPPs_rt];
    Q_N2C_PPs_rt_stack = [Q_N2C_PPs_rt_stack, Q_N2C_PPs_rt];
    sign_qe0_0_stack = [sign_qe0_0_stack; sign_qe0_0];
    omega_c_rt_stack = [omega_c_rt_stack; omega_c_rt(1:M_AOCS, :)];
    omegadot_c_rt_stack = [omegadot_c_rt_stack; omegadot_c_rt(1:M_AOCS, :)];
    Tc_AOCS_stack = [Tc_AOCS_stack; Tc_AOCS];
    xb_AOCS_stack = [xb_AOCS_stack; xb_AOCS];
    u_rt_norm_stack = [u_rt_norm_stack; u_rt_norm];
    u_AOCS_norm_stack = [u_AOCS_norm_stack; u_AOCS_norm];
    
    if verbose
        % fprintf('(%.0fVP) Branch %.0f: [%02d h, %02d m, %2.3f s] -> [%02d h, %02d m, %2.3f s]\n', length(t_viapoints_rt), branch, time_elapsed(t0_AOCS, t0, TU), time_elapsed(tf_AOCS, t0, TU));
        % fprintf('(%.0fVP) Branch %.0f: [%02d:%02d:%06.3f] -> [%02d:%02d:%06.3f]\n', length(t_viapoints_rt), branch, time_elapsed(t0_AOCS, t0, TU), time_elapsed(tf_AOCS, t0, TU));
        disp(sprintf('(%.0fVP) Branch %.0f: <strong>[%02d:%02d:%06.3f]</strong> -> <strong>[%02d:%02d:%06.3f]</strong>', length(t_viapoints_rt), branch, time_elapsed(t0_AOCS, t0, TU), time_elapsed(tf_AOCS, t0, TU)));
    end


    % In-Step Visualization
    if debug

        close all

        % fprintf('Propagating from t0 = %.4f min to tf = %.4f min\n', ([t0_AOCS, tf_AOCS]-t0)*TU/60)

        figure('Name', strcat("Branch ", string(branch), " - Natural Trajectory"))
        DrawTrajLVLH3D(TCC_rt(:, 7:9)*DU);

        figure('name', strcat("Branch ", string(branch), " - Trajectory"))
        for k = 1 : length(indices_ctrl) - 1
            bot_index = sum(indices_ctrl(1:k));
            top_index = sum(indices_ctrl(1:k+1))-1;
            testPPs(RHOdPPsLVLH_T(1:3, k), tspan_ctrl(bot_index:top_index));
        end
        DrawTrajLVLH3D(Y_ctrl(:, 7:9)*DU);

        figure('name', strcat("Branch ", string(branch), " - Body and Commanded Attitude"))
        subplot(1, 2, 1)
        plot((tspan_ctrl-t0)*TU*sec2hrs, Y_ctrl(:, 14), 'LineWidth', 1.5)
        hold on
        plot((tspan_ctrl-t0)*TU*sec2hrs, Y_ctrl(:, 15), 'LineWidth', 1.5)
        plot((tspan_ctrl-t0)*TU*sec2hrs, Y_ctrl(:, 16), 'LineWidth', 1.5)
        plot((tspan_ctrl-t0)*TU*sec2hrs, Y_ctrl(:, 17), 'LineWidth', 1.5)
        xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 12)
        ylabel('$q_i$', 'interpreter', 'latex', 'fontsize', 12)
        legend('q_{b0}', 'q_{b1}', 'q_{b2}', 'q_{b3}', 'fontsize', 10, 'location', 'best')
        grid on
        subplot(1, 2, 2)
        plot((tspan_ctrl-t0)*TU*sec2hrs, Q_N2C_AOCS_stack(:, 1), 'LineWidth', 1.5)
        hold on
        plot((tspan_ctrl-t0)*TU*sec2hrs, Q_N2C_AOCS_stack(:, 2), 'LineWidth', 1.5)
        plot((tspan_ctrl-t0)*TU*sec2hrs, Q_N2C_AOCS_stack(:, 3), 'LineWidth', 1.5)
        plot((tspan_ctrl-t0)*TU*sec2hrs, Q_N2C_AOCS_stack(:, 4), 'LineWidth', 1.5)
        xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 12)
        ylabel('$q_i$', 'interpreter', 'latex', 'fontsize', 12)
        legend('q_{c0}', 'q_{c1}', 'q_{c2}', 'q_{c3}', 'fontsize', 10, 'location', 'best')
        grid on  

        figure('name', strcat("Branch ", string(branch), " - Error Quaternions"))
        plot((tspan_ctrl-t0)*TU*sec2hrs, Qe_AOCS_stack(:, 1), 'LineWidth', 1.5)
        hold on
        plot((tspan_ctrl-t0)*TU*sec2hrs, Qe_AOCS_stack(:, 2), 'LineWidth', 1.5)
        plot((tspan_ctrl-t0)*TU*sec2hrs, Qe_AOCS_stack(:, 3), 'LineWidth', 1.5)
        plot((tspan_ctrl-t0)*TU*sec2hrs, Qe_AOCS_stack(:, 4), 'LineWidth', 1.5)
        xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 12)
        ylabel('$q_i$', 'interpreter', 'latex', 'fontsize', 12)
        legend('q_{e0}', 'q_{e1}', 'q_{e2}', 'q_{e3}', 'fontsize', 10, 'location', 'best')
        grid on

        figure('name', strcat("Branch ", string(branch), " - Natural vs AOCS Control Components"))
        subplot(1, 3, 1)
        plot((tspan_ctrl-t0)*TU*sec2hrs, u_rt_AOCS_stack(:, 1)*1000*DU/TU^2, 'LineWidth', 1.5);
        hold on
        plot((tspan_ctrl-t0)*TU*sec2hrs, u_AOCS_stack(:, 1)*1000*DU/TU^2, 'LineWidth', 1.5);
        xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 12)
        ylabel('$[m/s^2]$', 'interpreter', 'latex', 'fontsize', 12)
        legend('$u_{r, ideal}$', '$u_{r, AOCS}$', 'Location', 'best', 'Fontsize', 12, 'Interpreter','latex');
        grid on
        subplot(1, 3, 2)
        plot((tspan_ctrl-t0)*TU*sec2hrs, u_rt_AOCS_stack(:, 2)*1000*DU/TU^2, 'LineWidth', 1.5);
        hold on
        plot((tspan_ctrl-t0)*TU*sec2hrs, u_AOCS_stack(:, 2)*1000*DU/TU^2, 'LineWidth', 1.5);
        xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 12)
        ylabel('$[m/s^2]$', 'interpreter', 'latex', 'fontsize', 12)
        legend('$u_{\theta, ideal}$', '$u_{\theta, AOCS}$', 'Location', 'best', 'Fontsize', 12, 'Interpreter','latex');
        grid on
        subplot(1, 3, 3)
        plot((tspan_ctrl-t0)*TU*sec2hrs, u_rt_AOCS_stack(:, 3)*1000*DU/TU^2, 'LineWidth', 1.5);
        hold on
        plot((tspan_ctrl-t0)*TU*sec2hrs, u_AOCS_stack(:, 3)*1000*DU/TU^2, 'LineWidth', 1.5);
        xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 12)
        ylabel('$[m/s^2]$', 'interpreter', 'latex', 'fontsize', 12)
        legend('$u_{h, ideal}$', '$u_{h, AOCS}$', 'Location', 'best', 'Fontsize', 12, 'Interpreter','latex');
        grid on

        figure('name', strcat("Branch ", string(branch), " - omega_c"))
        plot((tspan_ctrl-t0)*TU*sec2hrs, omega_c_rt_stack(:, 1)/TU, 'LineWidth', 1.5)
        hold on
        plot((tspan_ctrl-t0)*TU*sec2hrs, omega_c_rt_stack(:, 2)/TU, 'LineWidth', 1.5)
        plot((tspan_ctrl-t0)*TU*sec2hrs, omega_c_rt_stack(:, 3)/TU, 'LineWidth', 1.5)
        xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 12)
        ylabel('$\omega_{c,i} \ [rad/s]$', 'interpreter', 'latex', 'fontsize', 12)
        legend('\omega_{c1}', '\omega_{c2}', '\omega_{c3}', 'fontsize', 10, 'location', 'best')
        grid on


        figure('name', strcat("Branch ", string(branch), " - omega_c"))
        plot((tspan_ctrl-t0)*TU*sec2hrs, omegadot_c_rt_stack(:, 1)/TU^2, 'LineWidth', 1.5)
        hold on
        plot((tspan_ctrl-t0)*TU*sec2hrs, omegadot_c_rt_stack(:, 2)/TU^2, 'LineWidth', 1.5)
        plot((tspan_ctrl-t0)*TU*sec2hrs, omegadot_c_rt_stack(:, 3)/TU^2, 'LineWidth', 1.5)
        xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 12)
        ylabel('$\dot{\omega}_{c,i} \ [rad/s^2]$', 'interpreter', 'latex', 'fontsize', 12)
        legend('\omega_{c1}', '\omega_{c2}', '\omega_{c3}', 'fontsize', 10, 'location', 'best')
        grid on


        figure('name', strcat("Branch ", string(branch), " - Tc"))
        plot((tspan_ctrl-t0)*TU*sec2hrs, Tc_AOCS_stack(:, 1)*1e6*DU^2/TU^2*MU, 'LineWidth', 1.5)
        hold on
        plot((tspan_ctrl-t0)*TU*sec2hrs, Tc_AOCS_stack(:, 2)*1e6*DU^2/TU^2*MU, 'LineWidth', 1.5)
        plot((tspan_ctrl-t0)*TU*sec2hrs, Tc_AOCS_stack(:, 3)*1e6*DU^2/TU^2*MU, 'LineWidth', 1.5)
        xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 12)
        ylabel('$T_{c,i} \ [Nm]$', 'interpreter', 'latex', 'fontsize', 12)
        legend('T_{c1}', 'T_{c2}', 'T_{c3}', 'fontsize', 10, 'location', 'best')
        grid on


        figure('name', strcat("Branch ", string(branch), " - xb_AOCS"))
        plot((tspan_ctrl-t0)*TU*sec2hrs, xb_AOCS_stack(:, 1), 'LineWidth', 1.5)
        hold on
        plot((tspan_ctrl-t0)*TU*sec2hrs, xb_AOCS_stack(:, 2), 'LineWidth', 1.5)
        plot((tspan_ctrl-t0)*TU*sec2hrs, xb_AOCS_stack(:, 3), 'LineWidth', 1.5)
        xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 12)
        ylabel('$x_{b,i}$', 'interpreter', 'latex', 'fontsize', 12)
        legend('x_{b1}', 'x_{b2}', 'x_{b3}', 'fontsize', 10, 'location', 'best')
        grid on
        

        if opt.additional_plots

            figure('name', strcat("Branch ", string(branch), " - xc"))
            plot((tspan_ctrl-t0)*TU*sec2hrs, xc_plot_stack(:, 1), 'LineWidth', 1.5)
            hold on
            plot((tspan_ctrl-t0)*TU*sec2hrs, xc_plot_stack(:, 2), 'LineWidth', 1.5)
            plot((tspan_ctrl-t0)*TU*sec2hrs, xc_plot_stack(:, 3), 'LineWidth', 1.5)
            xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 12)
            ylabel('$x_{c,i}$', 'interpreter', 'latex', 'fontsize', 12)
            legend('x_{c1}', 'x_{c2}', 'x_{c3}', 'fontsize', 10, 'location', 'best')
            grid on

            figure('name', strcat("Branch ", string(branch), " - yc"))
            plot((tspan_ctrl-t0)*TU*sec2hrs, yc_plot_stack(:, 1), 'LineWidth', 1.5)
            hold on
            plot((tspan_ctrl-t0)*TU*sec2hrs, yc_plot_stack(:, 2), 'LineWidth', 1.5)
            plot((tspan_ctrl-t0)*TU*sec2hrs, yc_plot_stack(:, 3), 'LineWidth', 1.5)
            xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 12)
            ylabel('$y_{c,i}$', 'interpreter', 'latex', 'fontsize', 12)
            legend('y_{c1}', 'y_{c2}', 'y_{c3}', 'fontsize', 10, 'location', 'best')
            grid on

            figure('name', strcat("Branch ", string(branch), " - zc"))
            plot((tspan_ctrl-t0)*TU*sec2hrs, zc_plot_stack(:, 1), 'LineWidth', 1.5)
            hold on
            plot((tspan_ctrl-t0)*TU*sec2hrs, zc_plot_stack(:, 2), 'LineWidth', 1.5)
            plot((tspan_ctrl-t0)*TU*sec2hrs, zc_plot_stack(:, 3), 'LineWidth', 1.5)
            xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 12)
            ylabel('$z_{c,i}$', 'interpreter', 'latex', 'fontsize', 12)
            legend('z_{c1}', 'z_{c2}', 'z_{c3}', 'fontsize', 10, 'location', 'best')
            grid on
    
            figure('name', strcat("Branch ", string(branch), " - Body and Wheels Angular Velocity"))
            subplot(1, 2, 1)
            plot((tspan_ctrl-t0)*TU*sec2hrs, Y_ctrl(:, 18)/TU, 'LineWidth', 1.5)
            hold on
            plot((tspan_ctrl-t0)*TU*sec2hrs, Y_ctrl(:, 19)/TU, 'LineWidth', 1.5)
            plot((tspan_ctrl-t0)*TU*sec2hrs, Y_ctrl(:, 20)/TU, 'LineWidth', 1.5)
            xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 12)
            ylabel('$\omega_i \, [rad/s]$', 'interpreter', 'latex', 'fontsize', 12)
            legend('\omega_{1}', '\omega_{2}', '\omega_{3}', 'fontsize', 10, 'location', 'best')
            grid on
            subplot(1, 2, 2)
            plot((tspan_ctrl-t0)*TU*sec2hrs, Y_ctrl(:, 21)/TU, 'LineWidth', 1.5)
            hold on
            plot((tspan_ctrl-t0)*TU*sec2hrs, Y_ctrl(:, 22)/TU, 'LineWidth', 1.5)
            plot((tspan_ctrl-t0)*TU*sec2hrs, Y_ctrl(:, 23)/TU, 'LineWidth', 1.5)
            plot((tspan_ctrl-t0)*TU*sec2hrs, Y_ctrl(:, 24)/TU, 'LineWidth', 1.5)
            xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 12)
            ylabel('$\omega_{si} \, [rad/s]$', 'interpreter', 'latex', 'fontsize', 12)
            legend('\omega_{s1}', '\omega_{s2}', '\omega_{s3}', '\omega_{s4}', 'fontsize', 10, 'location', 'best')
            grid on

            figure('name', strcat("Branch ", string(branch), " - Thrust Norms"))
            plot((tspan_ctrl-t0)*TU*sec2hrs, u_rt_norm_stack*1000*DU/TU^2, 'LineWidth', 1.5)
            hold on
            plot((tspan_ctrl-t0)*TU*sec2hrs, u_AOCS_norm_stack*1000*DU/TU^2, 'LineWidth', 1.5)
            xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 12)
            ylabel('$|u| \, [m/s^2]$', 'interpreter', 'latex', 'fontsize', 12)
            legend('traj', 'aocs', 'fontsize', 10, 'location', 'best')
            grid on

        end
        
        % pause

    end    

    % Set next Initial Conditions
    t0_rt = tf_AOCS;
    Y0_rt = Y_AOCS(end, :)';
    Yf_rt = Y_AOCS(end, :)';
    if tf_AOCS >= t0_backdrift
        break
    end

end

% Store Final Mass Ratio Value
TCC_ctrl = Y_ctrl(:, 1:13);
x7f = Y_ctrl(end, 13);


%% Final Natural Drift

% Define Forward Drift Propagation Parameters
tspan_drift = linspace(tspan_ctrl(end), tf + abs(tf-tspan_ctrl(end)), opt.N);                   % define the Forward Drift tspan
TC0_drift = TCC_ctrl(end, 1:12);                                                                % define Initial Conditions
optODE_fwd = odeset('RelTol', 1e-9, 'AbsTol', 1e-9, 'Events', @(t, Y) driftstop_fwd(t, Y));     % define options

% Perform the Final Natural Drift of the Chaser
[tspan_drift, TC_drift, t_fwdstop, ~, ~] = ode113(@(t, TC) NaturalRelativeMotion(t, TC, EarthPPsMCI, SunPPsMCI, ...
    muE, muS, MoonPPsECI, deltaE, psiM, deltaM, t0, tf, omegadotPPsLVLH, opt.show_progress), tspan_drift, TC0_drift, optODE_fwd);

if opt.show_progress
    close(pbar)
end

% Check Successful Forward Drift
if isempty(t_fwdstop)
    warning('Could not complete Forward Natural Drift with the given Initial Conditions.');
end


% Stack the States
TCC_drift = [TC_drift, x7f*ones(length(tspan_drift), 1)];
TCC = [TCC_ctrl_DA; TCC_ctrl; TCC_drift];
tspan = [tspan_ctrl_DA; tspan_ctrl; tspan_drift];


%% Post-Processing

% Retrieve States from Propagation
MEEt = TCC(:, 1:6);
RHO_LVLH = TCC(:, 7:12);
x7 = TCC(:, 13);

% Retrieve Target MCI State
COEt = MEE2COE(MEEt);
Xt_MCI = COE2rvPCI(COEt, muM);

M = length(tspan);

M_ctrl_DA = length(tspan_ctrl_DA);
M_ctrl = length(tspan_ctrl);
M_drift = length(tspan_drift);

% Initialize Post-Processing Variables
RHO_MCI = zeros(M, 6);
RHOd_LVLH = zeros(M_ctrl_DA + M_ctrl, 9);
Xc_MCI = zeros(M, 6);
dist = zeros(M, 1);
vel = zeros(M, 1);
u = zeros(M, 3);
u_norms = zeros(M, 1);
f_norms = zeros(M, 1);
kp_store = zeros(M, 1);
qe0 = zeros(M_ctrl, 1);
qe = zeros(M_ctrl, 3);
q0_LVLHt2MCI = zeros(M_ctrl, 1);
q_LVLHt2MCI = zeros(M_ctrl, 3);
Tc = zeros(M_ctrl, 3);

acc = zeros(M, 3);


% Perform Post-Processing
if opt.show_progress
    pbar = waitbar(0, 'Performing the Final Post-Processing');
end

for i = 1 : size(RHO_LVLH, 1)
    
    % Compute MCI States
    RHO_MCI(i, :) = rhoLVLH2MCI(RHO_LVLH(i, :)', Xt_MCI(i, :)', tspan(i), EarthPPsMCI, SunPPsMCI, MoonPPsECI, muE, muS, deltaE, psiM, deltaM)';
    Xc_MCI(i, :) = Xt_MCI(i, :) + RHO_MCI(i, :);
    
    % Compute the Distance and Velocity norms
    dist(i) = norm(RHO_MCI(i, 1:3));
    vel(i) = norm(RHO_MCI(i, 4:6));

    % Direct Approach: Hybrid-Predictive Feedback Control Law
    if i <= M_ctrl_DA
        
        [dY, ~, ~, ~, u(i, :), ~, ~, ~, ~, f_norms(i), kp_store(i), ~] = ...
            HybridPredictiveControlPostProcessing(tspan_ctrl_DA(i), TCC_ctrl_DA(i, :), EarthPPsMCI, SunPPsMCI, muM, ...
            muE, muS, MoonPPsECI, deltaE, psiM, deltaM, omegadotPPsLVLH, t0, tf, RHOdPPsLVLH_DA, DU, TU, stopSaturationTime, misalignment);
        RHOd_LVLH(i, :) = ppsval(RHOdPPsLVLH_DA, tspan(i));
        u_norms(i) = norm(u(i, :));    
    
    % Terminal Trajectory: Natural Feedback Control Law
    elseif i > M_ctrl_DA && i <= M_ctrl_DA + M_ctrl

        s = i - M_ctrl_DA;          % auxiliary index

        % [dY, ~, ~, ~, u(i, :), ~, ~, ~, f_norms(i)] = ...       % compute control values
        %     NaturalFeedbackControl(tspan_ctrl(s), TCC_ctrl(s, :), EarthPPsMCI, SunPPsMCI, muM, ...
        %     muE, muS, MoonPPsECI, deltaE, psiM, deltaM, omegadotPPsLVLH, t0, tf, RHOdPPsLVLH_T, kp, DU, TU, misalignment, 0, 0);
        
        for k = 1 : length(indices_ctrl) - 1        % retrieve reference trajectory
            bot_index = sum(indices_ctrl(1:k));
            top_index = sum(indices_ctrl(1:k+1))-1;
            if s >= bot_index && s <= top_index
                RHOd_LVLH(i, :) = ppsval(RHOdPPsLVLH_T(:, k), tspan(i));
                break
            end
        end
        
        % AOCS Reconstruction
        [~, ~, ~, ~, u(i, :), ~, ~, ~, f_norms(i), Tc(s, :), ~] = ...
            AOCS(tspan_ctrl(s), Y_ctrl(s, :)', EarthPPsMCI, SunPPsMCI, muM, ...
            muE, muS, MoonPPsECI, deltaE, psiM, deltaM, omegadotPPsLVLH, t0, tf, RHOdPPsLVLH_T(:, k), kp, omega_n, DU, TU, MU, TCC_PPs_stack(:, k), omega_cPPs_rt_stack(:, k), omegadot_cPPs_rt_stack(:, k), Q_N2C_PPs_rt_stack(:, k), sign_qe0_0_stack(k), misalignment, 0, 1);
            
        kp_store(i) = kp;
        u_norms(i) = norm(u(i, :));
        
        qc0_post = Q_N2C_AOCS_stack(s, 1);          % compute error quaternions
        qc_post = Q_N2C_AOCS_stack(s, 2:4)';
        qb0_post = Y_ctrl(s, 14);
        qb_post = Y_ctrl(s, 15:17)';
        qe0(s) = qc0_post*qb0_post + qc_post' * qb_post;
        qe(s, :) = (-qc_post*qb0_post + qc0_post*qb_post - skew(qc_post)*qb_post)';

        R_LVLHt2MCI_post = get_rotLVLH2MCI(Xt_MCI(i, :)', tspan(i), EarthPPsMCI, SunPPsMCI, MoonPPsECI, muE, muS, deltaE, psiM, deltaM);
        % [q0_LVLHt2MCI(s), q_LVLHt2MCI(s, :)] = C2q(R_LVLHt2MCI_post);
        [q0_LVLHt2MCI(s), q_LVLHt2MCI(s, :)] = C2q(R_LVLHt2MCI_post'); % I though it was the opposite but it feels like it's working better with the transpose

    % Terminal Trajectory: Natural Drift
    elseif i > M_ctrl_DA + M_ctrl
        
        s = i - M_ctrl_DA - M_ctrl;     % auxiliary index
        [dY, ~, ~, ~, f] = NaturalRelativeMotion(tspan_drift(s), TC_drift(s,:)', EarthPPsMCI, SunPPsMCI, ...
            muE, muS, MoonPPsECI, deltaE, psiM, deltaM, t0, tf, omegadotPPsLVLH, 0);
        f_norms(i) = norm(f);

    end

    acc(i, :) = dY(10:12);
    
    % Update Progress Bar
    if opt.show_progress
        waitbarMessage = sprintf('Post-Processing Progress: %.2f%%\n', i/M*100);
        waitbar(i/M, pbar, waitbarMessage);      % update the waitbar
    end

end

if opt.show_progress
    close(pbar)
end


% Create .json file for rendering
Q_LVLH2MCI = [q0_LVLHt2MCI, q_LVLHt2MCI];
QB_LVLH = quatmultiply(Y_ctrl(:, 14:17), Q_LVLH2MCI);
renderdata = [Y_ctrl(:, 7:9)*DU*1e3, QB_LVLH];

desState = [RHOf_LVLH(1:3)'*DU*1e3, RHOf_LVLH(4:6)'*DU/TU*1e3];
finalState = [TC_drift(end, 7:9)*DU*1e3, TC_drift(end, 10:12)*DU/TU*1e3];
deltaState = desState-finalState;

runtime = toc;

save(workspace_path);


end

