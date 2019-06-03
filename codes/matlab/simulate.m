function [ tout, yout, xout, uout ] = simulate( varargin )
% SIMULATE Simulates a system in open-loop or closed-loop conditions.
%   [TOUT, YOUT, XOUT, UOUT] = SIMULATE(MODEL,T,U,X_0) simulates the free 
%   response of a system represented by the object MODEL, for a time T 
%   and initial state X_0, given an input signal U. Returns the time TOUT, 
%   the model response  YOUT, the system real response XOUT and the 
%   input signal UOUT. 
%   
%   SIMULATE(MODEL,T,R,X_0,'control',CONTROLLER) simulates the closed-loop 
%   of the system MODEL with the controller defined by CONTROLLER and
%   reference signal given by R.
%   
%   SIMULATE(...,'noise',[w; z]) adds white Gaussian noise to represent
%   model and measurement uncertainty.
%   
%   SIMULATE(...,'disturbance',W) adds the disturbance signal W to the
%   appropriate model free parameters.

    % == FREE RESPONSE SIMULATION ==
    if length(varargin) == 4
        model = varargin{1}; t = varargin{2}; u = varargin{3}; x_0 = varargin{4}; 

        tout = t; uout = u;
        yout = zeros(numel(x_0), numel(t));
        yout(:,1) = x_0;
        for i = 1:numel(t)-1
            [~, y_aux] = odeSolver(model, t(i:i+1), u(:,i), yout(:,i), 100); 
            yout(:,i+1) = y_aux(end, :);
        end

        xout = yout;

    % == CONTROLLER SIMULATION ==
    elseif length(varargin) > 4
        % == STORING POSSIBLE PARAMETERS ==
        i = 1;
        while(5+i*2 < length(varargin))
            if(strcmpi(varargin{5+2*i}, "noise"))
                noise = varargin{6+2*i};
            elseif(strcmpi(varargin{5+2*i}, "disturbance"))
                W = varargin{6+2*i}; W = W';
            end
            i = i + 1;
        end

        % == CLOSED-LOOP SIMULATION ==
        model = varargin{1}; t = varargin{2}; r = varargin{3}; x_0 = varargin{4}; 
        controller = varargin{6};

        % == Auxiliary variables for the simulation ==
        deltaT = t(2) - t(1); % Size of the time step
        nx = model.sizeX;     % Dimension of the states
        nu = model.sizeU;     % Dimension of the inputs
        ny = model.sizeY;     % Dimension of the outputs

        ue = controller.oper.ue;         % Input operation values for the linear model
        xe = controller.oper.xe;         % State operation values for the linear model
        
        % == Retrives the correspondent State-Space matrices for the model selected index ==
        A = model.ss_model.A_l(xe, ue);  
        B = model.ss_model.B_l(xe, ue); 
        C = model.ss_model.C;                      
        D = model.ss_model.D;

        % == Separates the process and measurement noises ==
        if exist('noise', 'var')
            w = noise(1:nx, :)';
            z = noise(nx+1:end, :)';
        else
            w = zeros(nx,numel(t));
            z = zeros(ny,numel(t));
        end

        % Simulation variables
        u        = zeros(nu, numel(t)); 
        x        = zeros(nx, numel(t)); 
        x_hat    = zeros(nx, numel(t));

        % Re-scale the reference signal to the linearized system output
        r = r - C * xe;

        % Initial States
        x(:,1)     = x_0; 
        x_hat(:,1) = x_0' - xe;

        % == LINEAR QUADRATIC REGULATOR ==
        if(strcmpi(controller.type, 'lqr'))
            % Calculate the LQR gain matrices
            K = lqr_(A, B, controller.Q, controller.R, 0:deltaT:controller.N); 

            % == Simulation for each timestep (i) on time signal (t) ==
            for i = 1:numel(t)-1
                % Calculates the input signal
                u(:,i) = r(:,i) - K{i}*(x(:, i) - xe);

                % Simulates the actual plant (represented by the nonlinear model)
                if exist('W', 'var')
                    [~, y_aux] = odeSolver(model.model_W, t(i:i+1), [u(:,i) + ue; W(:,i)], x(:,i), 100);
                else
                    [~, y_aux] = odeSolver(model.model, t(i:i+1), u(:,i) + ue, x(:,i), 100);
                end
                
                x(:, i+1) = y_aux(end, :);

                % Defines the closed-loop matrix for this instant
                A_cl = (A - B * K{i});

                % Simulates the linear model using the Lagrange formula
                sys = ss(A_cl, B, eye(nx), zeros(nx,nu));
                aux = lsim(sys, r(:,i:i+1), t(i:i+1), x_hat(:,i));
                x_hat(:, i+1) = aux(end, :);
            end

        % == LINEAR QUADRATIC REGULATOR WITH INTEGRAL ACTION ==
        elseif(strcmpi(controller.type, 'lqri'))
            % Augments the state auxiliary vector in time and updates the initial conditions
            x_hat = zeros(nx+ny, numel(t));
            x_hat(:,1) = [x_0' - xe; r(:,1) - C*(x_0' - xe)];

            x_a = zeros(ny, numel(t));
            x_a(:,1) = r(:,1) - C*(x_0' - xe);

            % Augment the state and input matrices
            A_i = [A, zeros(nx, ny); 
                   -C, zeros(ny)];
            B_i = [B; 
                   zeros(ny, nu)];
            C_i = [C, zeros(ny, ny)];
            D_i = D;

            % Calculate the LQRI gain matrix
            K = lqr_(A_i, B_i, controller.Q, controller.R, 0:deltaT:controller.N);

            % Calculate the closed-loop augmented matrix
            B_cl = [zeros(nx, ny); eye(ny)];

            % == Simulation for each timestep (i) on time signal (t) ==
            for i = 1:numel(t)-1
                % Calculates the input signal
                u(:,i) = - K{i} * [x(:, i)-xe; x_a(:,i)];

                % Simulates the actual plant (represented by the nonlinear model)
                if exist('W', 'var')
                    [~, y_aux] = odeSolver(model.model_W, t(i:i+1), [u(:,i) + ue; W(:,i)], x(:,i), 100);
                else
                    [~, y_aux] = odeSolver(model.model, t(i:i+1), u(:,i) + ue, x(:,i), 100);
                end
                
                x(:, i+1) = y_aux(end, :);
                x_a(:, i+1) = x_a(:, i) + (r(:,i+1) - C*(x(:,i+1) - xe));

                % Defines the closed-loop matrix for this instant
                A_cl = [A-B*K{i}(:,1:nx) -B*K{i}(:,nx+1:end); 
                            -C                 zeros(ny)];

                % Simulates the linear model
                sys = ss(A_cl, B_cl, eye(nx+ny), zeros(nx+ny,ny));
                aux = lsim(sys, r(:,i:i+1), t(i:i+1), x_hat(:,i));
                x_hat(:, i+1) = aux(end, :);
            end

        % == LINEAR QUADRATIC GAUSSIAN REGULATOR ==,
        elseif(strcmpi(controller.type, 'lqg'))
            % Calculate the LQR gain matrix
            K = lqr_(A, B, controller.Q, controller.R, 0:deltaT:controller.N);  

            % Calculates the Kalman-Bucy gain matrix
            Ke = kalmanBucy(A, C, controller.Q_k, controller.R_k, 0:deltaT:controller.N); 

            % == Simulation for each timestep (i) on time signal (t) ==
            for i = 1:numel(t)-1
                % Calculates the input signal
                u(:,i) = r(:,i) - K{i}*x_hat(:, i);

                % Simulates the actual plant (represented by the nonlinear model)
                if exist('W', 'var')
                    [~, y_aux] = odeSolver(model.model_W, t(i:i+1), [u(:,i) + ue; W(:,i)], x(:,i), 100);
                else
                    [~, y_aux] = odeSolver(model.model, t(i:i+1), u(:,i) + ue, x(:,i), 100);
                end
                
                x(:, i+1) = y_aux(end, :);

                % Defines the closed-loop matrix for this instant
                A_cl = (A - B*K{i} - Ke{i}*C);
                B_cl = [B Ke{i}];

                % Simulates the linear model using the Lagrange formula
                sys = ss(A_cl, B_cl, eye(size(A_cl)), zeros(size(B_cl)));
                aux = lsim(sys, [r(:,i:i+1); C*(x(:,i:i+1)-xe)+z(:,i:i+1)], t(i:i+1), x_hat(:,i));
                x_hat(:, i+1) = aux(end, :);
            end

        % == LINEAR QUADRATIC GAUSSIAN REGULATOR WITH INTEGRAL ACTION ==,
        elseif(strcmpi(controller.type, 'lqgi'))
            % Augments the state auxiliary vector in time and updates the initial conditions
            x_hat = zeros(nx+ny, numel(t));
            x_hat(:,1) = [x_0' - xe; r(:,1) - C*(x_0' - xe) ];

            % Augment the state and input matrices
            A_i = [A, zeros(nx, ny); 
                  -C, zeros(ny)];
            B_i = [B; 
                   zeros(ny, nu)];
            C_i = [C, zeros(ny, ny)];
            D_i = D;
    
            % Calculate the LQRI gain matrix
            K = lqr_(A_i, B_i, controller.Q, controller.R, 0:deltaT:controller.N);

            % Generates and augments the Kalman observer state-space
            Ke = kalmanBucy(A, C, controller.Q_k, controller.R_k, 0:deltaT:controller.N);

            % == Simulation for each timestep (i) on time signal (t) ==
           for i = 1:numel(t)-1
                % Calculates the input signal
                u(:,i) = - K{i} * x_hat(:,i);

                % Simulates the actual plant (represented by the nonlinear model)
                if exist('W', 'var')
                    [~, y_aux] = odeSolver(model.model_W, t(i:i+1), [u(:,i) + ue; W(:,i)], x(:,i), 100);
                else
                    [~, y_aux] = odeSolver(model.model, t(i:i+1), u(:,i) + ue, x(:,i), 100);
                end
                
                x(:, i+1) = y_aux(end, :);

                % Define the closed-loop matrices
                A_cl = [A-B*K{i}(:,1:nx)-Ke{i}*C  -B*K{i}(:,nx+1:end); 
                           zeros(ny,nx)                zeros(ny)];
                B_cl = [zeros(nx, ny)     Ke{i}; 
                         eye(ny)        -eye(ny)];

                % Simulates the linear model
                sys = ss(A_cl, B_cl, eye(size(A_cl)), zeros(size(B_cl)));
                aux = lsim(sys, [r(:,i:i+1); C*(x(:,i:i+1)-xe)+z(:,i:i+1)], t(i:i+1), x_hat(:,i));
                x_hat(:, i+1) = aux(end, :);
           end

        % == SWITCHING MODE LINEAR QUADRATIC REGULATOR CONTROLLER ==,
        elseif(strcmpi(controller.type, 'switch-lqr'))
            % Define matrices to the operation states and inputs
            xe_t = zeros(nx, numel(t)); xe_t(:, 1) = xe;
            ue_t = zeros(nu, numel(t)); ue_t(:, 1) = ue;

            % Create the time-varying matrices structure
            A = zeros(nx, nx, numel(t));
            B = zeros(nx, nu, numel(t));

            A(:,:,1) = model.ss_model.A_l(xe, ue);
            B(:,:,1) = model.ss_model.B_l(xe, ue);

            % Unormalize the reference signal
            r = r + C*xe_t(:,1);

            % == Simulation for each timestep (i) on time signal (t) ==
            for i = 1:numel(t)-1
                % Defines the discrete matrices associated to state-space matrices
                ss_d = c2d(ss(A(:,:,i), B(:,:,i), C, D), deltaT);
                A_d = ss_d.A; B_d = ss_d.B;

                % Calculate the finite-horizon LQR gain matrix
                K = lqrd_(A_d, B_d, Q, R, N); 

                % Calculates the input signal (if still on control horizon)
                if(i < N)
                    u(:,i) = (r(:,i) - C*xe_t(:,i)) - K{1}*x_hat(:, i) ;
                end

                % Simulates the real physical plant using the non-linear model
                [~, y_aux] = odeSolver(model.model, t(i:i+1), u(:,i) + ue_t(:,i), y_aux(end, :), 100);
                x(:, i+1) = y_aux(end, :);

                % Simulates the linear model using the discrete differences equation
                %       x_{k+1} = A x_k + B u_k 
                x_hat(:, i+1) = (A_d - B_d*K{1}) * x_hat(:, i) + B_d*(r(:,i) - C*xe_t(:,i));

                % Switch the linear model and updates the previous state
                xe_t(:,i+1) = xe_t(:,i) + x_hat(:, i+1);
                ue_t(:,i+1) = (model.param.K1 * xe_t(1,i+1) - model.param.K2 * xe_t(2,i+1)) / xe_t(2,i+1);

                A(:,:,i+1) = model.ss_model.A_l(xe_t(:,i+1), ue_t(:,i+1));
                B(:,:,i+1) = model.ss_model.B_l(xe_t(:,i+1), ue_t(:,i+1));
            end

            % Updates the last state
            x_hat = x_hat + xe;
            u = u + ue;

        % == SWITCHING MODE LINEAR QUADRATIC GAUSSIAN CONTROLLER ==,
        elseif(strcmpi(controller.type, 'switch-lqg'))
            % Define matrices to the operation states and inputs
            xe_t = zeros(nx, numel(t)); xe_t(:, 1) = xe;
            ue_t = zeros(nu, numel(t)); ue_t(:, 1) = ue;

            % Unormalize the reference signal
            r = r + C * xe;

            % == Simulation for each timestep (i) on time signal (t) ==
            for i = 1:numel(t)-1
                % Defines the discrete state-space matrices
                ss_d = c2d(ss(A, B, C, D), deltaT);
                A_d = ss_d.A; B_d = ss_d.B; C_d = ss_d.C; D_d = ss_d.D;         

                K = lqr_(A_d, B_d, Q, R, N);

                % Calculates the input signal (if still on control horizon)
                if(i <= N)
                    u(:,i) = K(i, :) * (C_d * (r(:,i) - xe_t(:,i)) - C_d * x_hat(:, i) );
                end

                % Simulates the real physical plant using the non-linear model
                [~, y_aux] = odeSolver(model.model, t(i:i+1), u(:,i) + ue(:,i) + w(i, :), y_aux(end, :), 100);
                x(:, i+1) = y_aux(end, :) + z(i+1, :);

                % Kalman Filter's Prediction
                x_hat(:, i+1) = A_d*x_hat(:, i) + B_d*u(:,i);
                P_k = A_d * P_k * A_d' + Q_k;

                % Kalman Filter's Update
                S_k = C_d * P_k * C_d' + R_k;
                K_k = P_k * C_d' / S_k;
                err = C_d * x(:,i+1) - C_d * (x_hat(:, i+1) + xe(:,i));

                x_hat(:, i+1) = x_hat(:, i+1) + K_k * err;
                P_k = P_k - K_k * S_k * K_k';

                % Switch the linear model and updates the previous state
                xe(:,i+1) = xe(:,i) + x_hat(:, i+1);
                ue(:,i+1) = (model.param.K1 * xe(1,i+1) - model.param.K2 * xe(2,i+1)) / xe(2, i+1);

                A = model.ss_model.A_l(xe(:,i+1), ue(:,i+1));
                B = model.ss_model.B_l(xe(:,i+1), ue(:,i+1));
            end

            % Updates the last state
            x_hat = x_hat + xe;
            u = u + ue;

        end

        tout = t;
        uout = u;
        xout = x;
        yout = x_hat(1:nx, :);
    end
end

