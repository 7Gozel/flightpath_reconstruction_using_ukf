function y = y_obs_less_bias(ts, x, yNoise)

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

V     = sqrt(u^2+v^2+w^2);
alpha = atan(w/u);
beta  = asin(v/V);

% Output values
y(1) = phi     + yNoise(1);
y(2) = theta   + yNoise(2);
y(3) = psi     + yNoise(3);
y(4) = V       + yNoise(4);
y(5) = alpha   + yNoise(5);
y(6) = beta    + yNoise(6);
y(7) = x_e     + yNoise(7);
y(8) = y_e     + yNoise(8);
y(9) = h_e     + yNoise(9);
  
% y must be a column vector
y = y';

return
% end of function
