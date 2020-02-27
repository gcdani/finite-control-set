clear; clc;

MIT_motor.rs = 0.531; 
MIT_motor.rr = 0.408;
MIT_motor.J = 0.1; 
MIT_motor.Poles = 4;
MIT_motor.Lls = 2.52e-3; 
MIT_motor.Llr = 2.52e-3; 
MIT_motor.Lm = 84.7e-3;
MIT_motor.Power = 5;

MIT_input.voltage = 220;
MIT_input.freq = 60;
MIT_input.tsim = 1;
MIT_input.dt = 1e-5;
MIT_input.reference_frame = 2;
MIT_input.Tl = 0.83;
MIT_input.Tl_start = 1;
MIT_input.Tl_end = 1.2;
MIT_input.voltage_option = 0; % 0 = no pwm, 1 = VSI, 2 = FCS

MIT_output = MIT_sim(MIT_input, MIT_motor);

plot(MIT_output.Time,MIT_output.Speed);
title('Vds');
ylabel('[V]');

% hold on
% MIT_input.voltage_option = 0;
% MIT_input.voltage = 220;
% MIT_output = MIT_sim(MIT_input, MIT_motor);
% 
% plot(MIT_output.Time,MIT_output.Torque);