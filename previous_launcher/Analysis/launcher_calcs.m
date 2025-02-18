clear all
% football parameters
polar_moment = 0.0501*0.1; % polar moment of inertia: kg*m^2
a0 = 0.089; 
b0 = 0.089;
c0 = 0.1395;
a1 = 0.085;
b1 = 0.085;
c1 = 0.136;
leather_density = 1000;
mass = leather_density*(4/3)*pi*(a0*b0*c0 - a1*b1*c1)*0.1
m_ball = 0.425243;
r_ball = 0.084;

% launcher parameters
wheel_speed_rpm = 1000;
omega_w1 = wheel_speed_rpm*2*pi()/60 % wheel 1 angular velocity: rad/s
omega_w2 = wheel_speed_rpm*2*pi()/60; % wheel 2 angular velocity: rad/s
r_wheel = 0.125; % radius of wheels: meters
angle_w1 = 22.5; % angle between wheel and horizontal: deg
angle_w2 = -22.5; % angle between wheel and horizontal: deg
efficiency = 0.85; % efficiency due to slip/friction
angle_launch = 25; % angle that ball is launched: deg
cont_time = 0.05; % contact time: sec


v_w1 = omega_w1 * r_wheel; % wheel 1 linear velocity: m/s
v_w2 = omega_w2 * r_wheel; % wheel 2 linear velocity: m/s
v_ball = sqrt(((v_w1*cosd(angle_w1) + v_w2*cosd(angle_w2))*efficiency)^2 ...
    + ((v_w1*sind(angle_w1) + v_w2*sind(angle_w2))*efficiency)^2) % velocity of ball: m/s
dist_traveled = (v_ball^2*sind(2*angle_launch))/9.81 % distance traveled: m
v_ball_1 = sqrt((v_ball*cosd(angle_w1))^2 ...
    + (v_ball*sind(angle_w1))^2); % velocity imparted on ball from wheel 1: m/s
v_ball_2 = sqrt((v_ball*cosd(angle_w2))^2 ...
    + (v_ball*sind(angle_w2))^2); % velocity imparted on ball from wheel 2: m/s
force_1 = mass*v_ball_1/cont_time*sind(angle_w1); % y force from wheel 1 on ball: N
force_2 = mass*v_ball_2/cont_time*sind(angle_w2); % y force from wheel 2 on ball: N
omega_ball = (cont_time*(abs(force_1)*r_wheel ...
    + abs(force_2)*r_wheel))/polar_moment; % angular velocity of ball: rad/s
rpm = omega_ball/(2*pi)*60 % rpm of ball


force_applied_ball = mass*v_ball_1/cont_time;



%% plot trajectory (w air resistance)
a= pi*r_ball^2;
rho = 1.296 ; % air density kg/m^3
C_d = 0.2;
v_traj = v_ball; % m/s
cannon_angle = angle_launch; % DEGREES
g = 9.8; % m/s^2
x = []; x(1) = 0;
y = []; y(1) = 0;
v_x = v_traj*cosd(cannon_angle);
v_y = v_traj*sind(cannon_angle);
drag_const = 1/2*rho*C_d*a ;
delta_t = 0.01; % s
t = []; t(1) = delta_t;
i = 1;
while y(i) >= 0
    % calculate acceleration and velocity
    a_x = -(drag_const/m_ball)*v_traj*v_x;
    a_y = -g-(drag_const/m_ball)*v_traj*v_y;
    v_x = v_x + a_x*delta_t ;
    v_y = v_y + a_y*delta_t ;
    v = sqrt(v_x^2 + v_y^2);
    % calculate position
    x(i+1) = x(i) + v_x*delta_t + a_x*delta_t^2/2;
    y(i+1) = y(i) + v_y*delta_t + a_y*delta_t^2/2;
   
    % calculate loop values
    t(i+1) = t(i) + delta_t;
    i = i+1; 
end
figure()
plot(x,y)
axis("equal")
title("football Trajecotry")
xlabel("vertical travel (m)")
ylabel("distance traveled (m)")

ball_distance = max(x)


%% calcualte ideal angle
angles = linspace(15,60,45);
ball_dist = zeros(1,45);
for j = 1:45
    cannon_angle = angles(j); % DEGREES
    g = 9.8; % m/s^2
    x = []; x(1) = 0;
    y = []; y(1) = 0;
    v_x = v_traj*cosd(cannon_angle);
    v_y = v_traj*sind(cannon_angle);
    drag_const = 1/2*rho*C_d*a ;
    delta_t = 0.01; % s
    t = []; t(1) = delta_t;
    i = 1;
    while y(i) >= 0
        % calculate acceleration and velocity
        a_x = -(drag_const/m_ball)*v_traj*v_x;
        a_y = -g-(drag_const/m_ball)*v_traj*v_y;
        v_x = v_x + a_x*delta_t ;
        v_y = v_y + a_y*delta_t ;
        v = sqrt(v_x^2 + v_y^2);
        % calculate position
        x(i+1) = x(i) + v_x*delta_t + a_x*delta_t^2/2;
        y(i+1) = y(i) + v_y*delta_t + a_y*delta_t^2/2;
       
        % calculate loop values
        t(i+1) = t(i) + delta_t;
        i = i+1; 
    end
    ball_dist(j) = max(x);
end
[max_dist, index] = max(ball_dist);
ideal_angle = angles(index);

figure()
plot(angles, ball_dist)
title("ball distance v.s. launch angle")
xlabel("launcher angle (degrees)")
ylabel("distance traveled (m")

%% spin study

angles = linspace(0,50,100);
ball_dist = zeros(1,100);
spin_rate = zeros(1,100);
for j = 1:100
    cannon_angle = angles(j); % DEGREES

    wheel_speed_rpm = 1500;
    omega_w1 = wheel_speed_rpm*2*pi()/60; % wheel 1 angular velocity: rad/s
    omega_w2 = wheel_speed_rpm*2*pi()/60; % wheel 2 angular velocity: rad/s
    r_wheel = 0.125; % radius of wheels: meters
    angle_w1 = angles(j); % angle between wheel and horizontal: deg
    angle_w2 = -angles(j); % angle between wheel and horizontal: deg
    efficiency = 0.95; % efficiency due to slip/friction
    angle_launch = 25; % angle that ball is launched: deg
    cont_time = 0.05; % contact time: sec
    
    
    v_w1 = omega_w1 * r_wheel; % wheel 1 linear velocity: m/s
    v_w2 = omega_w2 * r_wheel; % wheel 2 linear velocity: m/s
    v_ball = sqrt(((v_w1*cosd(angle_w1) + v_w2*cosd(angle_w2))*efficiency)^2 ...
        + ((v_w1*sind(angle_w1) + v_w2*sind(angle_w2))*efficiency)^2) ;% velocity of ball: m/s
    ball_dist(j) = (v_ball^2*sind(2*angle_launch))/9.81; % distance traveled: m
    v_ball_1 = sqrt((v_ball*cosd(angle_w1))^2 ...
        + (v_ball*sind(angle_w1))^2); % velocity imparted on ball from wheel 1: m/s
    v_ball_2 = sqrt((v_ball*cosd(angle_w2))^2 ...
        + (v_ball*sind(angle_w2))^2); % velocity imparted on ball from wheel 2: m/s
    force_1 = mass*v_ball_1/cont_time*sind(angle_w1); % y force from wheel 1 on ball: N
    force_2 = mass*v_ball_2/cont_time*sind(angle_w2); % y force from wheel 2 on ball: N
    omega_ball = (cont_time*(abs(force_1)*r_wheel ...
        + abs(force_2)*r_wheel))/polar_moment; % angular velocity of ball: rad/s
    spin_rate(j) = omega_ball/(2*pi)*60*1.3 ;% rpm of ball


    force_applied_ball = mass*v_ball_1/cont_time;
end


figure()
plot(angles, ball_dist)
title("ball distance v.s. wheel angle")
xlabel("wheel angle (degrees)")
ylabel("distance traveled (m)")

figure()
plot(angles, spin_rate)
title("spin rate v.s. wheel angle")
xlabel("wheel angle (m)")
ylabel("spin rate (rpm)")

