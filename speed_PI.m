function [torque, control] = speed_PI(error, Ts, control)

Kp = 8000;                       % proportional constant
Ki = 40;                        % proportional constant
integral = control.integral;

Ns = 120 * 60 / 4;
Power = 5 / 1.341 * 1e3;
Tr = Power / ( Ns * pi / 30 );

% CONTROL LOOP
integral = integral + error*Ts;     % integral term
torque = Kp*error + Ki*integral; % action of control
if torque >= Tr
    torque = Tr;
elseif torque <= -Tr
    torque = -Tr;
end
control.integral = integral;

end