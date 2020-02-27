function  [lamb_s, lamb_r] = flux_estimation( vs, is, MIT_param, Ts, lamb_est )
% Ts = sampling time

    dlamb_s = Ts * ( vs - MIT_param.rs * is );
    lamb_s = lamb_est(1) + dlamb_s;
    
    Ls = MIT_param.Lm + MIT_param.Lls; Lr = MIT_param.Lm + MIT_param.Llr;
    lamb_r = Lr * lamb_est(1) / MIT_param.Lm + is * ( MIT_param.Lm - Lr * Ls / MIT_param.Lm);


end

