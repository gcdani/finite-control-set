function MIT_output = MIT_sim_step( MIT_input, MIT_motor, MIT_Initial_State )

% MIT_input.reference_frame: 0 = stationary, 1 = rotor, 2 = synchronous
rs = MIT_motor.rs; rr = MIT_motor.rr; J = MIT_motor.J; P = MIT_motor.Poles;
Lls = MIT_motor.Lls; Llr = MIT_motor.Llr; Lm = MIT_motor.Lm;

Ls = Lm + Lls; Lr = Lm + Llr;
wb = 2 * pi * MIT_input.freq;
Tl = 0;
tsim = MIT_input.tsim;
Power = MIT_motor.Power / 1.341 * 1e3;

Ns = 120 * MIT_input.freq / P;
Tr = Power / ( Ns * pi / 30 );

theta = MIT_Initial_State.theta;
Te = 0;
wr = MIT_Initial_State.wr;
theta_r = MIT_Initial_State.theta_r;
wm = MIT_Initial_State.wm;
dt = MIT_input.dt;
idx = 0;
vdr = 0;
vqr = 0;

L_matrix = [
    Ls 0 Lm 0;
    0 Ls 0 Lm;
    Lm 0 Lr 0;
    0 Lm 0 Lr];

iqd = MIT_Initial_State.iqd; % zeros(4,1);
lamb = MIT_Initial_State.lamb; % zeros(4,1);
dlamb = zeros(4,1);
r_matrix = [rs; rs; rr; rr] .* eye(4);

if MIT_input.reference_frame == 0 || MIT_input.reference_frame == 1
    w = 0;
elseif MIT_input.reference_frame == 2
    w = wb;
end

w_matrix = [
    0 -w 0 0;
    w  0 0 0;
    0  0 0 -w;
    0  0 w 0];

w_matrix(3,4) =  - w + wr;
w_matrix(4,3) =  - w_matrix(3,4);
w_matrix(1,2) =  - w;
w_matrix(2,1) =  - w_matrix(1,2);

iabcs = zeros(3, round(tsim / dt));
iabcr = zeros(3, round(tsim / dt));
vabcs = zeros(3,round(tsim / dt));
i_qd = zeros(4,round(tsim / dt));
time = zeros(1,round(tsim / dt));
speed = zeros(1,round(tsim / dt));
torque = zeros(1,round(tsim / dt));
load_torque = zeros(1,round(tsim / dt));
v_qd = zeros(2, round(tsim / dt));

for t = 0:dt:tsim
    
    idx = idx + 1;
    
    Tl = load_torque_create(MIT_input.Tl * Tr, MIT_input.Tl_start, MIT_input.Tl_end, t);
    
    voltage = voltage_input( MIT_input.voltage, wb, t, theta, MIT_input.voltage_option );
    
    vqd = voltage.vqd;
    
    dlamb = vqd - r_matrix * iqd + w_matrix * lamb;
    
    lamb = lamb + dt .* dlamb;
    
    iqd = L_matrix \ lamb;
    
    Te = 0.75 * P * Lm * (iqd(1) * iqd(4) - iqd(2) * iqd(3));
    
    dtheta_r = dt * wr;
    dwm = dt * (Te - Tl) / J;
    theta_r = dtheta_r + theta_r;
    wm = dwm + wm;
    wr = P * 0.5 * wm;
    
    if MIT_input.reference_frame == 0
        w = 0;
    elseif MIT_input.reference_frame == 1
        w = wr;
    elseif MIT_input.reference_frame == 2
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
    
    iabcs(:,idx) = invKs * iqd(1:2);
    vabcs(:,idx) = invKs * vqd(1:2);
    iabcr(:,idx) = invKr * iqd(3:4);
    i_qd(:,idx) = iqd;
    v_qd(:,idx) = vqd(1:2);
    time(idx) = t;
    speed(idx) = wm * 30 / pi;
    torque(idx) = Te;
    load_torque(idx) = Tl;
    
end

MIT_output.iabcs = iabcs;
MIT_output.iabcr = iabcr;
MIT_output.vabcs = vabcs;
MIT_output.iqd = i_qd;
MIT_output.vqd = v_qd;
MIT_output.Time = time;
MIT_output.Speed = speed;
MIT_output.Torque = torque;
MIT_output.LoadTorque = load_torque;

end

