clear; clc;

% Model
MIT_param.Vrms = 220;
MIT_param.f = 60;
MIT_param.rs = 0.531; MIT_param.rr = 0.408; MIT_param.J = 0.1;
MIT_param.P = 4;
MIT_param.Lls = 2.52e-3; MIT_param.Llr = 2.52e-3;
MIT_param.Lm = 84.7e-3;
MIT_param.Power = 5;

% Simulation parameters
MIT_sim.Tl = 0;
MIT_sim.dt = 1e-5;
MIT_sim.sample_time = 5e-5;
MIT_sim.reference_frame = 2; % 0 = stationary, 1 = rotor, 2 = synchronous
MIT_sim.vs = 220;
MIT_sim.vqd = 0;

% Initial states
MIT_init.iqd = zeros(4,1);
MIT_init.lamb = zeros(4,1);
MIT_init.wm = 0;
MIT_init.theta_r = 0;
MIT_init.theta = 0;
MIT_init.w_matrix = zeros(4);
MIT_init.init_time = 0;
MIT_init.is = 0;

tsim = 1;
tsim = tsim / MIT_sim.sample_time;
lamb_est = [0, 0];
rho = 0;
Te_ref = 0;
lamb_s_ref = 0.55;
tsim = 1;
tsim = tsim / MIT_sim.sample_time;

for k = 0
    
    MIT_sim.vs = 220;
    MIT_sim.vqd = 0;
    
    % Initial states
    MIT_init.iqd = zeros(4,1);
    MIT_init.lamb = zeros(4,1);
    MIT_init.wm = 0;
    MIT_init.theta_r = 0;
    MIT_init.theta = 0;
    MIT_init.w_matrix = zeros(4);
    MIT_init.init_time = 0;
    MIT_init.is = 0;
    speed = [];
    time = [];
    jopt = [];
    ls = [];
    wref = 188.4956;
    control.integral = 0;
    tref = [];
    te_error = [];
    te_est = [];
    torque = [];
    is = [];
    vs = [];

    
    for idx = 0:tsim
        
        [T_ref, control] = speed_PI(wref - MIT_init.wm, MIT_sim.sample_time, control);
        tref = [tref T_ref];
        % Estimation of lamb_s and lamb_r
        % [lamb_s, lamb_r] = flux_estimation( MIT_sim.vs, MIT_init.is, MIT_param, MIT_sim.sample_time, lamb_est );
        % ls = [ls lamb_s];
        
        % Best switch
        if k == 0
            fsc = FCS( 0.75, lamb_s_ref, MIT_param, MIT_sim.sample_time, rho, MIT_init, lamb_est);
            lamb_est = fsc.flux;
            MIT_sim.vqd = (fsc.vqd);
            MIT_sim.vs = (fsc.vs);
            jopt = [jopt fsc.jopt];
        else
            MIT_sim.vqd = voltage_input( MIT_param.Vrms, 2*pi*60, MIT_init.init_time, MIT_init.theta, 0 );
        end
        
        % Simulate until next sample
        [out, states] = MIT_step(MIT_param, MIT_sim, MIT_init);
        speed = [speed out.speed];
        time = [time out.time];
        torque = [torque out.torque];
        
        MIT_init.iqd = states.iqd;
        MIT_init.lamb = states.lamb;
        MIT_init.wm = states.wm;
        MIT_init.theta_r = states.theta_r;
        MIT_init.theta = states.theta;
        MIT_init.w_matrix = states.w_matrix;
        MIT_init.init_time = states.init_time;
        MIT_init.is = (states.is);
        is = [is MIT_init.is];
        if k == 0
            vs = [vs fsc.vs];
            Te_error = fsc.Te_p - out.torque(end);
            te_error = [te_error Te_error];
            te_est = [te_est fsc.Te_e];
        end
        
    end
    
    hold on
    plot(time, speed);
end
