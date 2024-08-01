% Reference:

% Wan, Eric, and Rudolph Van Der Merwe. The Unscented Kalman Filter. 2004.

% Chapter 7: Recursive Parameter Estimation 
% "Flight Vehicle System Identification - A Time Domain Methodology"
% Second Edition
% by Ravindra V. Jategaonkar
% published by AIAA, Reston, VA 20191, USAA

% Otávio, Bruno & Teixeira, Soares & Torres, Leonardo & Henriques, 
% Paulo & Oliveira, Iscold & Aguirre, Luis. (2005). 
% Flight path reconstruction using the unscented Kalman filter algorithm. 

clc
clear
% Read the CSV file
filename = 'example.csv';
data = readmatrix(filename);
% data = [Time	a_x fts^2	a_y fts^2	a_z fts^2	p degs	q degs	r degs	
% Roll deg	Pitch deg	Heading deg	Velocity kt	Alpha deg	Beta deg
% Latitude	Longitude	Altitude ft]

% Convert units
data(:,2:4) = data(:,2:4)*0.3048; % ft/s to m/s
data(:,5:10) = data(:,5:10)*pi/180; % deg to rad
data(:,11) = data(:,11)*0.514444; % kt to m/s
data(:,12:13) = data(:,12:13)*pi/180; %deg to rad
data(:,16) = data(:,16)*0.3048; % ft to m

origin = [data(1,14), data(1,15), data(1,16)];
[xEast,yNorth] = latlon2local(data(:,14),data(:,15),data(:,16),origin);
data(:,14:15) = [xEast,yNorth];

v_i = data(1,11);
a_i = data(1,12);
b_i = data(1,13);

uvw_i = [v_i*cos(a_i)*cos(b_i) v_i*sin(b_i) v_i*sin(a_i)*cos(b_i);];

% Initial guess for the biases
param = [1 1 1 1 1 1]*0.001; 

% Measured data
Z = data(:,8:16);

x_i  = [uvw_i, data(1, 8:10), data(1, 14:16) param];
% Step 2: Extract the time column and the rest of the data columns
data_time = data(:, 1); % First column is time
data_time = data_time-data_time(1);

u_input = data(:,2:7);
xNoise = zeros(9, 1);

alpha = 1;
beta = 2;
nx = 9;
nxp = length(param)+nx;
ny = nx;
L = nxp+nx+ny; % states+param(aug)+process noise+measurement noise
kappa = 3-L;
lamda=alpha^2*(L+kappa)-L;
gamma=sqrt(L+lamda);
nsp = 2*L+1;

% Compute weights for 2*L+1 sigma points (wm: for mean and wc for covariance)
Wm(1,1)     = lamda/(L+lamda);
Wm(2:nsp,1) = ones(2*L,1)/(2*(L+lamda));
Wc(1,1)     = Wm(1,1) + 1 - alpha^2 + beta;
Wc(2:nsp,1) = Wm(2:nsp,1);


rr   = [0.0002  0.002  0.0002  0.0002...
        0.0000000001  0.0000000001 1 1  1]/1000';

qq(1:nx,1)     = [0.005  0.004  0.002 0.004 0.004  0.004  0.1  0.1 0.1]*10';
qq(nx+1:nxp,1) = [zeros(nxp-nx,1)];
    
% Initial state propagation error covariance matrix - pcov(nxp): only diagonal terms
pa0(1:nx,1)     = 1;     % p0 for system states
pa0(nx+1:nxp,1) = 1;     % p0 for system biases

% State propagation error covariance matrix - pcov: (nxp X nxp)
pcov = diag(pa0);
% Process noise covariance matrix - Qnoise: (nx x nx)
Qnoise = diag(qq(1:nx));
% Measurement noise covariance matrix - Rnoise: (ny x ny)
Rnoise = diag(rr);

Pa = blkdiag(pcov, Qnoise, Rnoise);

Psqrtm = gamma*chol(Pa)';
xhat = x_i';
xa = [xhat; zeros(nx,1); zeros(ny,1)];

Xsp    = [zeros(L,1) -Psqrtm Psqrtm] + repmat(xa,1,nsp);
xtSigPt(1:nxp,:) = Xsp(1:nxp,:);                 % for the first data point


for k=1:length(data)
    disp(k)
    u1c = u_input(k,:)'; 
    dt = 0.1;
    % Augmented Covariance matrix and augmented state vector
    Pa = blkdiag(pcov, Qnoise, Rnoise);
    xa(:,1) = [xhat; zeros(nx,1); zeros(ny,1)];

    % Generate (2Na+1) sigma points:
    Psqrtm  = gamma*chol(Pa)';
    Xsp    = [zeros(L,1) -Psqrtm Psqrtm] + repmat(xa,1,nsp);
    
    %----------------------------------------------------------------------
    % Prediction step (Time update): Propagate sigma points through state functions
    if k > 1
    u1p = u_input(k-1,:)';
        
    dt = (data(k, 1)-data(k-1, 1));
        for ip=1:nsp
           xNoise(1:nx,1) = Xsp(nxp+1:nxp+nx,ip);
           xSP = ruku4(dt , Xsp(1:nxp,ip), xNoise, dt , u1p, u1c);
           xtSigPt(:,ip) = xSP;       
        end
    end
    
    %----------------------------------------------------------------------
    % System states as weighted sum
    xtilde = xtSigPt*Wm; 

    % Compute covariances of predicted states:
    xtSPdiff = xtSigPt - kron(xtilde,ones(1,nsp));
    Pxx      = zeros(nxp,nxp);
    for ip=1:nsp
        Pxx = Pxx + Wc(ip) * (xtSPdiff(:,ip)*xtSPdiff(:,ip)');
    end
    
    %----------------------------------------------------------------------
    % Correction step (Measurement update):
    
    % Compute model outputs; 
    ytSigPt= zeros(ny,nsp);
    for ip=1:nsp
        yNoise(1:ny,1) = Xsp(nxp+nx+1:L,ip);
        ySP = y_obs_less_bias(dt, xtSigPt(1:nxp,ip), yNoise);
        ytSigPt(:,ip)  = ySP;  
    end

    % Model output as weighted sum 
    ytilde = ytSigPt*Wm;
    
    % Compute covariances and cross correlations, Pyy and Pxy;
    ytSPdiff = ytSigPt - kron(ytilde,ones(1,nsp));
    Pxy = zeros(nxp,ny);
    Pyy = zeros(ny,ny);
    for ip=1:nsp
        Pyy = Pyy + Wc(ip) * ytSPdiff(:,ip) * ytSPdiff(:,ip)'; 
        Pxy = Pxy + Wc(ip) * xtSPdiff(:,ip) * ytSPdiff(:,ip)';
    end
         
    % Kalman gain matrix Kgain
    Kgain = Pxy/Pyy;  

    % Corrected states
    z    = Z(k,1:ny)';                           % Measured outputs
    res  = z - ytilde;                           % Output Residual
    xhat = xtilde + Kgain*res;                   % State Update

    % Covariance update P(k) = P(k-) - KPyy*K';
    pcov = Pxx - Kgain * Pyy * Kgain';
   
    % Save estimated states for plotting
    y_estimated(k,1:ny)  = ytilde';
    xhSave(k,1:nxp) = xhat';

end

close all

state_vector_names = {'\phi_m', '\theta_m', '\psi_m', 'V_{TAS,m}', '\alpha_m', '\beta_m', 'x_{E,m}', 'y_{E,m}', 'H_m'};

% Function to clean the state vector name for file saving
clean_name = @(name) regexprep(name, '[^a-zA-Z0-9]', '_');

% Ensure the plots folder exists
if ~exist('plots', 'dir')
   mkdir('plots');
end

for i = 1:9
    figure;
    plot( data_time, Z(:,i), 'b', data_time, y_estimated(:,i), 'r--');
    title(['State Vector Component: ', state_vector_names{i}]);
    xlabel('Time');
    ylabel(state_vector_names{i});
    legend('Measured','Estimated');
    grid on;
    
    % Clean the name for the filename
    cleaned_name = clean_name(state_vector_names{i});
    filename = ['plots/plot_', cleaned_name, '.png'];
    
    % Save the figure
    saveas(gcf, filename);
end

figure('Position', [100, 100, 800, 600]);
plot3(Z(:,7),Z(:,8),Z(:,9),"b")

hold on
grid on
ylabel("Y [meters]")
zlabel('Altitude [meters]');
xlabel('X [meters]');
plot3(y_estimated(:,7),y_estimated(:,8),y_estimated(:,9),"r--","LineWidth", 2)
lgd = legend('Measured', 'Estimated');
lgd.Location = 'northoutside';
lgd.Orientation = 'horizontal';
lgd.Position = [0.5, 0.95, 0, 0];
filename = "plots/3D_plot.png";
    
% Save the figure
saveas(gcf, filename);

% titles = {'\Delta a_x fts^2', '\Delta a_y fts^2', '\Delta a_z fts^2', '\Delta p degs', '\Delta q degs', '\Delta r degs'};
% 
% for i = 10:15
%     figure;
%     plot(data_time, xhSave(:,i))
%     grid on
%     xlabel('Time');
%     ylabel(titles{i-9});
%     legend('Bias');
% end

function xa = ruku4(ts, xa, xNoise, dt, u1, u2)

dt2  = 0.5*dt;

xadt  = xdot_less_bias(ts, xa, u1, xNoise);
rk1 = xadt*dt2;

u12 = (u1 + u2)/2;    % average of u(k) and u(k+1)
xa1 = xa + rk1;
xadt  = xdot_less_bias(ts, xa1, u12, xNoise);
rk2 = dt2*xadt;

xa1 = xa + rk2;
xadt  = xdot_less_bias(ts, xa1, u12, xNoise);
rk3 = dt2*xadt;

xa1 = xa + 2*rk3;
xadt  = xdot_less_bias(ts, xa1, u2, xNoise);
rk4 = dt2*xadt;

xa  = xa + (rk1 + 2.0*(rk2+rk3) + rk4)/3;
end
