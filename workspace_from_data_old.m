function [RHO_LVLH, M_ctrl_DA, M_ctrl, M_drift, DU, TU, RHOd_LVLH, dist, vel, ...
          renderdata, TCC, Xt_MCI, RHO_MCI, u, u_norms, f_norms, kp_store, ...
          qe0, qe, Tc, Ta, betas, gammas, acc, deltaState, tspan, tspan_ctrl, ...
          Y_ctrl, t0, tf, failure_times, misalignments, Y_drift, Q_N2C_drift, ...
          qe0_drift, qe_drift, Tc_drift, Ta_drift] = workspace_from_data_old(data)

RHO_LVLH = data.RHO_LVLH;
M_ctrl_DA = data.M_ctrl_DA;
M_ctrl = data.M_ctrl;
M_drift = data.M_drift;
DU = data.DU;
TU = data.TU;
RHOd_LVLH = data.RHOd_LVLH;
dist = data.dist;
vel = data.vel;
renderdata = data.renderdata;
TCC = data.TCC;
Xt_MCI = data.Xt_MCI;
RHO_MCI = data.RHO_MCI;
u = data.u;
u_norms = data.u_norms;
f_norms = data.f_norms;
kp_store = data.kp_store;
qe0 = data.qe0;
qe = data.qe;
Tc = data.Tc;
Ta = data.Ta;
betas = data.betas;
gammas = data.gammas;
acc = data.acc;
deltaState = data.deltaState;
tspan = data.tspan;
tspan_ctrl = data.tspan_ctrl;
Y_ctrl = data.Y_ctrl;
t0 = data.t0;
tf = data.tf;
failure_times = data.failure_times;
misalignments = data.misalignments;
Y_drift = data.Y_drift;
Q_N2C_drift = data.Q_N2C_drift;
qe0_drift = data.qe0_drift;
qe_drift = data.qe_drift;
Tc_drift = data.Tc_drift;
Ta_drift = data.Ta_drift;


end
