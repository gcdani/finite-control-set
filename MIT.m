clear; clc;

reference_frame = 2; % 0 = stationary, 1 = rotor, 2 = synchronous

Vrms = 220;
f = 60;
wb = 2 * pi * f;
rs = 0.531; rr = 0.408; J = 0.1; P = 4;
Lls = 2.52e-3; Llr = Lls; Lm = 84.7e-3;
Ls = Lm + Lls; Lr = Lm + Llr;
Tl = 0;
tsim = 1;
Power = 5;
Power = Power / 1.341 * 1e3;

Ns = 120 * f / P;
Tr = Power / ( Ns * pi / 30 );

theta = 0;
Te = 0;
wr = 0;
theta_r = 0;
wm = 0;
dt = 1e-5;
idx = 0;
vdr = 0;
vqr = 0;

L_matrix = [ 
    Ls 0 Lm 0;
    0 Ls 0 Lm;
    Lm 0 Lr 0;
    0 Lm 0 Lr];

iqd = zeros(4,1);
dlamb = zeros(4,1);
lamb = zeros(4,1);
r_matrix = diag([rs; rs; rr; rr]);

if reference_frame == 0 || reference_frame == 1
    w = 0;
elseif reference_frame == 2
    w = wb;
end
    
w_matrix = [
    0 -w 0 0;
    w  0 0 0;
    0  0 0 -w;
    0  0 w 0];

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
    
    if t >= 0.7 && t <= 0.8
        Tl = 0.83 * Tr;
    else
        Tl = 0;
    end
    
    vas = sqrt(2/3) * Vrms * cos(wb * t);
    vbs = sqrt(2/3) * Vrms * cos(wb * t - 2 * pi / 3);
    vcs = sqrt(2/3) * Vrms * cos(wb * t + 2 * pi / 3);
    
    Ks = [  cos(theta) cos(theta - 2 * pi / 3) cos(theta + 2 * pi / 3);
            sin(theta) sin(theta - 2 * pi / 3) sin(theta + 2 * pi / 3)];
    Ks = (2/3) .* Ks;
    
    vqd = Ks * [vas;vbs;vcs];
    vqd(3) = 0;
    vqd(4) = 0;
    
    dlamb = vqd - r_matrix * iqd + w_matrix * lamb;
    
    lamb = lamb + dt .* dlamb;
    
    iqd = L_matrix \ lamb;
    
    Te = 0.75 * P * Lm * (iqd(1) * iqd(4) - iqd(2) * iqd(3));
    
    dtheta_r = dt * wr;
    dwm = dt * (Te - Tl) / J;
    theta_r = dtheta_r + theta_r;
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

figure;
subplot(8,1,1);
plot(time,v_qd(2,:));
title('Vds');
ylabel('[V]');

subplot(8,1,2);
plot(time,v_qd(1,:));
title('Vqs');
ylabel('[V]');

subplot(8,1,3);
plot(time,i_qd(2,:));
title('Ids');
ylabel('[A]');

subplot(8,1,4);
plot(time,i_qd(1,:));
title('Iqs');
ylabel('[A]');

subplot(8,1,5);
plot(time,i_qd(4,:));
title('Idr');
ylabel('[A]');

subplot(8,1,6);
plot(time,i_qd(3,:));
title('Iqr');
ylabel('[A]');

subplot(8,1,7);
plot(time,speed);
title('Speed');
ylabel('[RPM]');

subplot(8,1,8);
plot(time,torque);
title('Torque');
xlabel('time [s]');
ylabel('[N.m]');