function [ tout, yout, xout, uout ] = simulate( varargin )
%SIMULATE Simulates a non-linear model using the Euler method.

    switch length(varargin)
        case 4
        %   @param model   The non-linear ODEs representing the system model.
        %   @param t       The initial conditions of the simulation.
        %   @param U       The input signal for the simulation.        
        %   @param X_0     The initial conditions of the simulation.
        model = varargin{1}; t = varargin{2}; u = varargin{3}; x_0 = varargin{4};
        
        tout = t; uout = u;
        yout = zeros(numel(x_0), numel(t));
        yout(:,1) = x_0;
        for i = 1:numel(t)-1
            [~, y_aux] = odeSolver(model, t(i:i+1), u(:,i), yout(:,i), 100); 
            yout(:,i+1) = y_aux(end, :);
        end
        
        xout = yout;

        case 5
        %   @param model   The non-linear ODEs representing the system model.
        %   @param t       The initial conditions of the simulation.
        %   @param U       The reference signal for the simulation.        
        %   @param X_0     The initial conditions of the simulation.
        %   @param K       The state feedback gain vector being used.
        %   @param L       The Luenberger gain vector being used.
        model = varargin{1}; t = varargin{2}; U = varargin{3}; X_0 = varargin{4};
        K = varargin{5}; L = varargin{6};

        a = 'TODO: nonlinear model control';
      
        case num2cell(9:1:11)
        %   @param model     The non-linear ODEs representing the system model.
        %   @param idx       The index position of the linear state-space model to use.
        %   @param t         The initial conditions of the simulation.
        %   @param r         The input signal for the simulation.        
        %   @param X_0       The initial conditions of the simulation.
        %   @param type      The type of controller being used to simulate.
        %   @param Q         The weighting matrix to the states signals.
        %   @param R         The weighting matrix to the input signals.
        %   @param N         The size of the control time horizon.
        %   @param L         The Luenberger gain vector being used.
        %   @param w         The matrix of additive disturbance in the nonlinear model.
        
        % == Saving function parameter as variables ==
        model = varargin{1}; oper = varargin{2}; t = varargin{3}; r = varargin{4}; X_0 = varargin{5}; 
        type = varargin{6}; Q = varargin{7}; R = varargin{8}; N = varargin{9};
      
        if(length(varargin) >= 10)
            w = varargin{10};
            z = varargin{11};
        else
            w = zeros(numel(t), model.sizeX);
            z = zeros(numel(t), model.sizeY);
        end
        
        % == Retrives the correspondent State-Space matrices for the model selected index ==
        A = model.ss_model.A_l(oper.xe, oper.ue);  B = model.ss_model.B_l(oper.xe, oper.ue); 
        C = model.ss_model.C;                      D = model.ss_model.D;
        
        % == Auxiliary variables for the simulation ==
        deltaT = t(2) - t(1); % Size of the time step
        nx = model.sizeX;     % Dimension of the states
        nu = model.sizeU;     % Dimension of the inputs
        ny = model.sizeY;     % Dimension of the outputs
        
        ue = oper.ue;       % Input operation values for the linear model
        xe = oper.xe;       % State operation values for the linear model
        
        % Simulation variables
        u        = zeros(nu, numel(t)); 
        x        = zeros(nx, numel(t)); 
        x_hat    = zeros(nx, numel(t));
        
        % Re-scale the reference signal to the linearized system output
        r = r - C * xe;
        
        % Initial States
        x(:,1)     = X_0' + w(1, :)'; 
        x_hat(:,1) = X_0' - xe;
        y_aux      = X_0;
        
        % == LINEAR QUADRATIC REGULATOR ==
        if(strcmp(type, 'lqr'))
            % Calculate the LQR gain matrices
            K = lqr_(A, B, Q, R, 0:deltaT:N); 

            % == Simulation for each timestep (i) on time signal (t) ==
            for i = 1:numel(t)-1
                % Calculates the input signal
                u(:,i) = r(:,i) - K{i}*x_hat(:, i);

                % Defines the closed-loop matrix for this instant
                A_cl = (A - B * K{i});
            
                % Simulates the linear model using the Lagrange formula
                %       x_k = e^{A (t - t_0)} x_0 + \int{ e^{A (t - \tau)} B u_\tau + L( Y_\tau - C x_\tau ) }
                sys = ss(A_cl, zeros(nx,nu), eye(nx), zeros(nx,nu));
                aux = lsim(sys, r(:,i:i+1), t(i:i+1), x_hat(:,i));
                x_hat(:, i+1) = aux(end, :);

            end
        
        % == LINEAR QUADRATIC REGULATOR WITH INTEGRAL ACTION ==
        elseif(strcmp(type, 'lqri'))
            % Augments the state auxiliary vector in time and updates the initial conditions
            x_hat = zeros(nx+ny, numel(t));
            x_hat(:,1) = [X_0' - xe; r(:,1) - C*(X_0' - xe)];

            % Augment the state and input matrices
            A_i = [A, zeros(nx, ny); 
                   -C, zeros(ny)];
            B_i = [B; 
                   zeros(ny, nu)];
            C_i = [C, zeros(ny, ny)];
            D_i = D;
            
            % Calculate the LQRI gain matrix
            K = lqr_(A_i, B_i, Q, R, 0:deltaT:N);

            % Calculate the closed-loop augmented matrix
            B_cl = [zeros(nx, ny); eye(ny)];

            % == Simulation for each timestep (i) on time signal (t) ==
            for i = 1:numel(t)-1
                % Calculates the input signal
                u(:,i) = - K{i} * x_hat(:, i);

                % Defines the closed-loop matrix for this instant
                A_cl = [A-B*K{i}(:,1:nx) -B*K{i}(:,nx+1:end); -C zeros(ny)];
            
                % Simulates the linear model
                sys = ss(A_cl, B_cl, eye(nx+ny), zeros(nx+ny,nu));
                aux = lsim(sys, r(:,i:i+1), t(i:i+1), x_hat(:,i));
                x_hat(:, i+1) = aux(end, :);

            end

        % == LINEAR QUADRATIC GAUSSIAN REGULATOR ==,
        elseif(strcmp(type, 'lqg'))
            % Generates the actual model for the plant
            realSys = ss(A, B, eye(nx), D_i);
            
            % Estimate the covariances for the process and measurement noises
            Q_k = cov(w);
            R_k = cov(z);
            
            % Calculate the LQR gain matrix
            K = lqr_(A, B, Q, R, 0:deltaT:N);  

            % Calculates the Kalman-Bucy gain matrix
            Ke = kalmanBucy(A, C, Q_k, R_k, 0:deltaT:N); 

            % == Simulation for each timestep (i) on time signal (t) ==
            for i = 1:numel(t)-1
                % Calculates the input signal
                u(:,i) = r(:,i) - K{i}*x_hat(:, i);

                % Defines the closed-loop matrix for this instant
                A_cl = (A - B*K{i} - Ke{i}*C);
                B_cl = [zeros(nx,ny) Ke{i}];
                
                % Simulates the real disturbed plant
                aux = lsim(realSys, u(:,i:i+1)+w(i,:), t(i:i+1), x(:,i));
                x(:, i+1) = aux(end, :) + z(i+1,:);

                % Simulates the linear model using the Lagrange formula
                %       x_k = e^{A (t - t_0)} x_0 + \int{ e^{A (t - \tau)} B u_\tau + L( Y_\tau - C x_\tau ) }
                sys = ss(A_cl, B_cl, eye(nx), D);
                aux = lsim(sys, [r(:,i:i+1); C*(x(:,i:i+1)-xe)], t(i:i+1), x_hat(:,i));
                x_hat(:, i+1) = aux(end, :);

            end
            
        % == LINEAR QUADRATIC GAUSSIAN REGULATOR WITH INTEGRAL ACTION ==,
        elseif(strcmp(type, 'lqgi'))
            % Generates the actual model for the plant
            realSys = ss(A, B, eye(nx), D);
            
            % Augments the state auxiliary vector in time and updates the initial conditions
            x_hat = zeros(nx+ny, numel(t));
            x_hat(:,1) = [X_0' - xe; r(:,1) - C*(X_0' - xe) ];

            Q_k = cov(w);
            R_k = cov(z);
            
            % Augment the state and input matrices
            A_i = [A, zeros(nx, ny); 
                  -C, zeros(ny)];
            B_i = [B; 
                   zeros(ny, nu)];
            C_i = [C, zeros(ny, ny)];
            D_i = D;

            % Calculate the LQRI gain matrix
            K = lqr_(A_i, B_i, Q, R, 0:deltaT:N);

            % Generates and augments the Kalman observer state-space
            Ke = kalmanBucy(A, C, Q_k, R_k, 0:deltaT:N);
            
            % Define the closed-loop matrices
            B_cl = [zeros(nx, ny); eye(ny)];

            % == Simulation for each timestep (i) on time signal (t) ==
           for i = 1:numel(t)-1
                % Calculates the input signal
                u(:,i) = - K{i} * x_hat(:, i);
                
                % Define the closed-loop matrices
                A_cl = [A-B*K{i}(1:nx)-Ke{i}*C  -B*K{i}(nx+1:end); 
                                -C                zeros(ny)];
                B_cl = [Ke{i} zeros(nx, ny); zeros(ny), eye(ny)];
                
                % Simulates the real disturbed plant
                aux = lsim(realSys, u(:,i:i+1), t(i:i+1), x(:,i));
                x(:, i+1) = aux(end, :) + w(i+1,:) + z(i+1,:);
                            
                % Simulates the linear model
                n_resp = expm( A_cl * (t(i+1) - t(1)) ) * x_hat(:,1);
                f_resp = 0;
                for j =1:i+1
                    f_resp = f_resp + deltaT * expm( A_cl * (t(i+1) - t(j)) ) * ( B_cl * [C*(x(:,j)-xe); r(:,j)] ) ;
                end

                x_hat(:, i+1) = n_resp + f_resp;
           end
            
        % == SWITCHING MODE LINEAR QUADRATIC REGULATOR CONTROLLER ==,
        elseif(strcmp(type, 'switch-lqr'))
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
        elseif(strcmp(type, 'switch-lqg'))
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

