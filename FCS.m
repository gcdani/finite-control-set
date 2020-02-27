function [output] = FCS( Te_ref, lamb_s_ref, MIT_param, Ts, rho, MIT_init, flux_est )

g_opt = inf;

w = 2 * pi * MIT_param.f;
t = MIT_init.init_time;
wt = w * t;
vdc = MIT_param.Vrms * sqrt(3 / 2);
wr = MIT_param.P * 0.5 * MIT_init.wm;
we = 0.5 * MIT_param.P * MIT_init.wm;
is = MIT_init.iqd(1) - 1i * MIT_init.iqd(2);
Lr = MIT_param.Lm + MIT_param.Llr;
Ls = MIT_param.Lm + MIT_param.Lls;
Lm = MIT_param.Lm;
rs = MIT_param.rs;
rr = MIT_param.rr;
sigma = (1 - Lm * Lm / (Ls*Lr));
L_sigma = Ls * sigma;

AlphaBeta2dq = [-sin(w * t) cos(w * t); cos(w * t) sin(w * t)]; 
Tc = (2/3) .* [0 -sqrt(3)/2 sqrt(3)/2; 1 -0.5 -0.5];
Sabc = [0 1 1 0 0 0 1 1;
        0 0 1 1 1 0 0 1;
        0 0 0 0 1 1 1 1];
kr = Lm / Lr;
R_sigma = MIT_param.rs + MIT_param.rr * kr * kr;
tau_sigma = L_sigma / R_sigma;
tau_r = Lr / MIT_param.rr;

 theta = MIT_init.theta;
 
 Ks = [  cos(theta) cos(theta - 2 * pi / 3) cos(theta + 2 * pi / 3);
            sin(theta) sin(theta - 2 * pi / 3) sin(theta + 2 * pi / 3)];
 Ks = (2/3) .* Ks;
 AlphaBeta2ABC = [1 0; -0.5 sqrt(3)/2; -0.5 -sqrt(3)/2];

for j = 0:6
    
    if j == 0 || j == 7
        vd = 0;
        vq = 0;
    elseif j == 1
        vd = 2 * vdc / 3;
        vq = 0;
    elseif j == 2
        vd = vdc / 3;
        vq = vdc / sqrt(3);
    elseif j == 3
        vd = - vdc / 3;
        vq = vdc / sqrt(3);
    elseif j == 4
        vd = - 2 * vdc / 3;
        vq = 0;
    elseif j == 5
        vd = - vdc / 3;
        vq = - vdc / sqrt(3);
    elseif j == 6
        vd = vdc / 3;
        vq = - vdc / sqrt(3);
    end
    
    vs = vd + 1i * vq;
    vqd = [vd; vq];
    
    flux_r = (flux_est(2) + Ts * ( rr * Lm * is / Lr - ( rr / Lr - 1i * we ) * flux_est(2) ));
    flux_s = (Lm * flux_est(2) / Lr + sigma * Ls * is);
    Te_est = 0.75 * MIT_param.P * imag( conj(flux_s) * is );
    
    flux_s_p = flux_s + Ts * vs - Ts * rs * is;
    is_p = (1 - Ts / tau_sigma) * is + (Ts / tau_sigma + Ts) * ( ( ( kr / tau_r - 1i*kr*we ) * flux_r + vs ) / R_sigma );
    Te_p =  0.75 * MIT_param.P * imag( conj(flux_s) * is_p );
    
    gt = abs( Te_ref - Te_p );
    glamb = abs( lamb_s_ref - flux_s );
    g = gt + rho * glamb;
    
    if t <= 0.4 && j == 0
        g = inf;
    end
    
    
    if g < g_opt
        g_opt = g;
        output.gopt = g;
        output.jopt = j;
        output.vs = vs;
        output.vqd = [vqd; 0; 0];
        output.Te_p = Te_p;
        output.lamb_s_p = flux_s;
        output.is_p = is;
        output.Te_e = Te_est;
        output.flux = [flux_s; flux_r];
    end 
    
end
