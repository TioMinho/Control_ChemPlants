function [x_hat, x] = luenberger(sys, ssPoints, L, u, t, x_0)
%LUENBERGER Simulates the system given an input signal and a observer gain
%   sys - State-Space system object
%   u   - Input signal vector (#elements = numel(t))
%   t   - Time vector
%   x_0 - Initial states
    
    % == Auxiliary variables for the simulation ==
    deltaX = t(2) - t(1); % Size of the time step
    nx = sys.sizeX;       % Dimension of the states
    xe = ssPoints(1:nx);       % State operation values for the linear model
    ue = ssPoints(nx+1:end);   % Input operation values for the linear model
    
    % == Retrives the correspondent State-Space matrices for the model selected index ==
    A = sys.ss_model.A_l(xe, ue);   B = sys.ss_model.B_l(xe, ue); 
    C = sys.ss_model.C;             D = sys.ss_model.D; 
    
    A_cl = A - L*C;
    
    % Simulation variables
    x       = zeros(nx, numel(t)); 
    x_hat   = zeros(nx, numel(t));

    % Initial States
    x(:,1)      = x_0' + 2; 
    x_hat(:,1)  = x_0';
%     y_aux       = x_0;

    % == Simulation for each timestep (i) on time signal (t) ==
    for i = 1:numel(t)-1
        % Simulates the real physical plant using the non-linear model
        n_resp = expm( A * (t(i+1) - t(1)) ) * x(:, 1);
        f_resp = 0;
        for j =1:i+1
            f_resp = f_resp + deltaX * expm( A * (t(i+1) - t(j)) ) * ( B * u(:, j) );
        end
        
        x(:, i+1) = n_resp + f_resp;

%         [~, y_aux] = odeSolver(sys.model, t(i:i+1), u(:,i)+ue, y_aux(end, :), 100);
%         x(:, i+1) = y_aux(end, :);

        % Simulates the linear model using the Lagrange
        % formulas
        %       x_k = e^{A (t - t_0)} x_0 + \int{ e^{A (t - \tau)} B u_\tau + L( Y_\tau - C x_\tau ) }
        n_resp = expm( A_cl * (t(i+1) - t(1)) ) * x_hat(:, 1);
        f_resp = 0;
        for j =1:i+1
            f_resp = f_resp + deltaX * expm( A_cl * (t(i+1) - t(j)) ) * ( B*u(:, j) + L*(x(:, j)) );
        end

        x_hat(:, i+1) = n_resp + f_resp;
    end
end

