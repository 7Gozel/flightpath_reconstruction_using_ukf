function xdot = xdot_less_bias(ts, x, u_input, xNoise)

% State variables
u     = x(1);
v     = x(2);
w     = x(3);
phi   = x(4);
theta = x(5);
psi   = x(6);
x_e   = x(7);
y_e   = x(8);
h_e   = x(9);

%Biases
delta_a_x   = x(10);
delta_a_y   = x(11);
delta_a_z   = x(12);
delta_p     = x(13);
delta_q     = x(14);
delta_r     = x(15);

W_x_e       = 0;
W_y_e       = 0;
W_h_e       = 0;

% Input variables
a_x = u_input(1);
a_y = u_input(2);
a_z = u_input(3);
p = u_input(4);
q = u_input(5);
r = u_input(6);

% Constants
g = 9.805416; %9.81

% Derivative of state variables
u_dot = -(q-delta_q)*w + (r-delta_r)*v - g*sin(theta)          + (a_x-delta_a_x);
v_dot = -(r-delta_r)*u + (p-delta_p)*w + g*cos(theta)*sin(phi) + (a_y-delta_a_y);
w_dot = -(p-delta_p)*v + (q-delta_q)*u + g*cos(theta)*cos(phi) + (a_z-delta_a_z);

phi_dot   = (p-delta_p) + (q-delta_q)*sin(phi)*tan(theta) + (r-delta_r)*cos(phi)*tan(theta);
theta_dot =               (q-delta_q)*cos(phi)            - (r-delta_r)*sin(phi);
psi_dot   =               (q-delta_q)*sin(phi)*sec(theta) + (r-delta_r)*cos(phi)*sec(theta);

L_EB = [ cos(theta)*cos(psi), sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi), cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi);
         cos(theta)*sin(psi), sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi), cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi);
         sin(theta),         -sin(phi)*cos(theta),                          -cos(phi)*cos(theta)                        ]; 

velocity_Earth = L_EB * [u; v; w] + [W_x_e; W_y_e; W_h_e];

x_e_dot = velocity_Earth(1);
y_e_dot = velocity_Earth(2);
h_e_dot = velocity_Earth(3);

% State derivatives
xdot(1) = u_dot     + xNoise(1);
xdot(2) = v_dot     + xNoise(2);
xdot(3) = w_dot     + xNoise(3);
xdot(4) = phi_dot   + xNoise(4);
xdot(5) = theta_dot + xNoise(5);
xdot(6) = psi_dot   + xNoise(6);
xdot(7) = x_e_dot   + xNoise(7);
xdot(8) = y_e_dot   + xNoise(8);
xdot(9) = h_e_dot   + xNoise(9);

% Derivative of the biases are assumed to be zero
xdot(10:length(x)) = 0; 

% xdot must be a column vector
xdot = xdot';

return
% end of function
