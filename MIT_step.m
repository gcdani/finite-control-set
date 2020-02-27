function [MIT_out, MIT_states] = MIT_step(MIT_param, MIT_sim, MIT_init)

% Aspectos contrutivos
Vrms = MIT_param.Vrms;
f = MIT_param.f;
rs = MIT_param.rs; rr = MIT_param.rr; J = MIT_param.J;
P = MIT_param.P;
Lls = MIT_param.Lls; Llr = MIT_param.Llr; 
Lm = MIT_param.Lm;
Power = MIT_param.Power;

% Aspectos da simulacao
Tl = MIT_sim.Tl;
dt = MIT_sim.dt;
sample_time = MIT_sim.sample_time;
reference_frame = MIT_sim.reference_frame;

% Initial states
iqd = MIT_init.iqd;
lamb = MIT_init.lamb;
wm = MIT_init.wm;
theta_r = MIT_init.theta_r;
theta = MIT_init.theta;
w_matrix = MIT_init.w_matrix;
init_time = MIT_init.init_time;

% Calculos extras
Ls = Lm + Lls; Lr = Lm + Llr;
wb = 2 * pi * f;
Ns = 120 * f / P;
Power = Power / 1.341 * 1e3;
Tr = Power / ( Ns * pi / 30 );
wr = P * 0.5 * wm;

% Initializing variables
dlamb = zeros(4,1);
Te = 0;
idx = 0;
vdr = 0;
vqr = 0;

L_matrix = [ 
    Ls 0 Lm 0;
    0 Ls 0 Lm;
    Lm 0 Lr 0;
    0 Lm 0 Lr];

r_matrix = diag([rs; rs; rr; rr]);

tsim = init_time + sample_time;

if reference_frame == 0
    w = 0;
elseif reference_frame == 1
    w = wr;
elseif reference_frame == 2
    w = wb;
end
    
w_matrix(3,4) =  - w + wr;
w_matrix(4,3) =  - w_matrix(3,4);
w_matrix(1,2) =  - w;
w_matrix(2,1) =  - w_matrix(1,2);

MIT_out.iabcs = zeros(3, round(sample_time / dt));
MIT_out.iabcr = zeros(3, round(sample_time / dt));
MIT_out.vabcs = zeros(3,round(sample_time / dt));
MIT_out.i_qd = zeros(4,round(sample_time / dt));
MIT_out.time = zeros(1,round(sample_time / dt));
MIT_out.speed = zeros(1,round(sample_time / dt));
MIT_out.torque = zeros(1,round(sample_time / dt));
MIT_out.load_torque = zeros(1,round(sample_time / dt));
MIT_out.v_qd = zeros(2, round(sample_time / dt));

AlphaBeta2ABC = [1 0; -0.5 sqrt(3)/2; -0.5 -sqrt(3)/2];

for t = init_time:dt:tsim
    
    idx = idx + 1;
    
%     vas = sqrt(2/3) * Vrms * cos(wb * t);
%     vbs = sqrt(2/3) * Vrms * cos(wb * t - 2 * pi / 3);
%     vcs = sqrt(2/3) * Vrms * cos(wb * t + 2 * pi / 3);
%     
%     Ks = [  cos(theta) cos(theta - 2 * pi / 3) cos(theta + 2 * pi / 3);
%             sin(theta) sin(theta - 2 * pi / 3) sin(theta + 2 * pi / 3)];
%     Ks = (2/3) .* Ks;
%     
%     vqd = Ks * [vas;vbs;vcs];
%     vqd(3) = 0;
%     vqd(4) = 0;

    
    vqd = MIT_sim.vqd;
    
    dlamb = vqd - r_matrix * iqd + w_matrix * lamb;
    
    lamb = lamb + dt .* dlamb;
    
    iqd = L_matrix \ lamb;
    
    Te = 0.75 * P * Lm * (iqd(1) * iqd(4) - iqd(2) * iqd(3));
    
%     if t >= 0.8 && t <= 1
%         Tl = 10;
%     end
%     
    dtheta_r = dt * wr;
    theta_r = dtheta_r + theta_r;
    dwm = dt * (Te - Tl) / J;
    wm = dwm + wm;
    wr = P * 0.5 * wm;
    
    if reference_frame == 0
        w = 0;
    elseif reference_frame == 1
        w = wr;
    elseif reference_frame == 2
        w = wb;
    end
    
    theta = theta + w * dt;
    
    w_matrix(3,4) =  - w + wr;
    w_matrix(4,3) =  - w_matrix(3,4);
    w_matrix(1,2) =  - w;
    w_matrix(2,1) =  - w_matrix(1,2);
    
    invKs = [cos(theta) sin(theta); 
             cos(theta - 2 * pi / 3) sin(theta - 2 * pi / 3);
             cos(theta + 2 * pi / 3) sin(theta + 2 * pi / 3)];
    
    beta = theta - theta_r;
         
    invKr = [cos(beta) sin(beta);
             cos(beta - 2 * pi / 3) sin(beta - 2 * pi / 3);
             cos(beta + 2 * pi / 3) sin(beta + 2 * pi / 3)];
    
    MIT_out.iabcs(:,idx) = invKs * iqd(1:2);
    MIT_out.vabcs(:,idx) = invKs * vqd(1:2);
    MIT_out.iabcr(:,idx) = invKr * iqd(3:4);
    MIT_out.i_qd(:,idx) = iqd;
    MIT_out.v_qd(:,idx) = vqd(1:2);
    MIT_out.time(idx) = t;
    MIT_out.speed(idx) = wm * 30 / pi;
    MIT_out.torque(idx) = Te;
    MIT_out.load_torque(idx) = Tl;
    
end

% Initial states
MIT_states.iqd = iqd;
MIT_states.lamb = lamb;
MIT_states.wm = wm;
MIT_states.theta_r = theta_r;
MIT_states.theta = theta;
MIT_states.w_matrix = w_matrix;
MIT_states.init_time = tsim;
% 
% AlphaBeta2dq = [-sin(wb * t) cos(wb * t); cos(wb * t) sin(wb * t)]; 
% ialpha_beta = AlphaBeta2dq \ iqd(1:2);
T = (2 / 3) .* [1 -0.5 -0.5; 0 sqrt(3)/2 -sqrt(3)/2];
ialpha_beta = T * MIT_out.iabcs(:,end);

MIT_states.is = ialpha_beta(1) + 1i * ialpha_beta(2);