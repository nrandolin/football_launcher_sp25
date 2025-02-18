% Hand calcs for force on hub/roll pin
% 420 SS / 1070-1905 Steel 5/32 roll pin has max double shear force of 2750
% lb
f_ball = 475; % N
r_wheel = 0.15; % M
motor_shaft_dia = 0.008; %mm

t_outer_wheel = 475*.15; % N-M

f_roll_pin = t_outer_wheel/0.008; % N

f_roll_pin_lbf = f_spring_pin * .225; % lbf

