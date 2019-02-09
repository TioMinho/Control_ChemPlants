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
      
        case num2cell(10:1:12)
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
        model = varargin{1}; idx = varargin{2}; t = varargin{3}; r = varargin{4}; X_0 = varargin{5}; 
        type = varargin{6}; Q = varargin{7}; R = varargin{8}; N = varargin{9}; L = varargin{10};
      
        if(length(varargin) >= 11)
            w = varargin{11};
            z = varargin{12};
        else
            w = zeros(numel(t), size(model.sizeU, 1));
            z = zeros(numel(t), size(model.sizeY, 1));
        end
        
        % == Retrives the correspondent State-Space matrices for the model selected index ==
        A = model.ss_model.A(idx); B = model.ss_model.B(idx); C = model.ss_model.C; D = model.ss_model.D;
        
        % == Auxiliary variables for the simulation ==
        deltaX = t(2) - t(1);   % Size of the time step
        nx = model.sizeX;     % Dimension of the states
        nu = model.sizeU;    % Dimension of the inputs
        ny = model.sizeY;     % Dimension of the outputs
        
        ue = model.oper.U(idx,:)';       % Input operation values for the linear model
        xe = model.oper.X(idx,:)';        % State operation values for the linear model
        
        % Simulation variables
        u        = zeros(nu, numel(t)); 
        x        = zeros(nx, numel(t)); 
        x_hat = zeros(nx, numel(t));
        
        % Re-scale the reference signal to the linearized system output
        r = r - C * xe;
        
        % Initial States
        x(:,1)        = X_0' + z(1, :)'; 
        x_hat(:,1) = X_0' - xe;
        y_aux = X_0;
        
        % == LINEAR QUADRATIC REGULATOR ==
        if(strcmp(type, 'lqr'))
           
            % Case for the Finite-Horizon Discrete-Time
            if(N ~= 'inf') 
                % Defines the discrete matrices associated to state-space matrices
                ss_d = c2d(ss(A, B, C, D), deltaX);
                A_d = ss_d.A; B_d = ss_d.B; C_d = ss_d.C; D_d = ss_d.D; 
                
                % Calculate the finite-horizon LQR gain matrix
                K = lqr_(A_d, B_d, Q, R, N);  

                % == Simulation for each timestep (i) on time signal (t) ==
                for i = 1:numel(t)-1
                    % Calculates the input signal (if still on control horizon)
                    if(i < N)
                        u(:,i) = K(i,:) * ( r(:,i) - C_d * x_hat(:, i) );
                    end

                    % Simulates the real physical plant using the non-linear model
                    [~, y_aux] = odeSolver(model.model, t(i:i+1), u(:,i) + ue + w(i, :), y_aux(end, :), 100);
                    x(:, i+1) = y_aux(end, :) + z(i+1, :);

                    % Simulates the linear model using the discrete differences equation
                    %       x_{k+1} = A x_k + B u_k + L ( Y_k - C x_k )
                    x_hat(:, i+1) = A_d  * x_hat(:, i) + B_d * u(:, i) + L * ( x(:, i) - C_d*(x_hat(:, i) + xe) );
                    
                end

            % Case for the Infinite-Horizon Continuous-Time
            else
                % Calculate the LQR gain matrix
                K = lqr_(A, B, Q, R, N);  
                
                % Define the closed-loop matrices
                A_cl = (A - B * K);
                B_cl = (B * K);
                
                % == Simulation for each timestep (i) on time signal (t) ==
                for i = 1:numel(t)-1
                    % Calculates the input signal
                    u(:,i) = K * ( r(:,i) - C * x_hat(:, i) );

                    % Simulates the real physical plant using the non-linear model
                    [~, y_aux] = odeSolver(model.model, t(i:i+1), u(:,i) + ue + w(i, :), y_aux(end, :), 100);
                    x(:, i+1) = y_aux(end, :) + z(i+1, :);

                    % Simulates the linear model using the Lagrange
                    % formulas
                    %       x_k = e^{A (t - t_0)} x_0 + \int{ e^{A (t - \tau)} B u_\tau + L( Y_\tau - C x_\tau ) }
                    n_resp = expm( A * (t(i+1) - t(1)) ) * x_hat(:, 1);
                    f_resp = 0;
                    for j =1:i+1
                        f_resp = f_resp + deltaX * expm( A * (t(i+1) - t(j)) ) * ( B * u(:, j) + L *(x(:, j) - C * (x_hat(:, j) + xe)) );
                    end

                    x_hat(:, i+1) = n_resp + f_resp;
                    
                end
            end
        
        % == LINEAR QUADRATIC REGULATOR WITH INTEGRAL ACTION ==
        elseif(strcmp(type, 'lqri'))
            % Augments the Luenberger matrix
            L_i = [L; zeros(ny, ny)];
            
            % Case for a Finite-Horizon Discrete-Time
            if(N ~= 'inf') 
                % Defines the discrete state-space matrices
                ss_d = c2d(ss(A, B, C, D), deltaX);
                A_d = ss_d.A; B_d = ss_d.B; C_d = ss_d.C; D_d = ss_d.D;
                
                 % Augment the state and input matrices
                A_i = [A, zeros(nx, ny); 
                          -C, eye(ny)];
                B_i = [B; 
                           zeros(ny, nu)];
                C_i = [C, zeros(ny, ny)];
                D_i = D;

                % Calculate the finite-horizon LQR gain matrix
                K = lqr_(A_i, B_i, Q, R, N);
                
                % == Simulation for each timestep (i) on time signal (t) ==
                for i = 1:numel(t)-1
                    % Calculates the input signal (if still on control horizon)
                    if(2 <= i && i < N)
                        u(:,i) = u(:,i-1) - K(i,:) * [x_hat(:, i-1) - x_hat(:, i); r(:,i) - C_d * x_hat(:,i)];
                    end

                    % Simulates the real physical plant using the non-linear model
                    [~, y_aux] = odeSolver(model.model, t(i:i+1), u(:,i) + ue + w(i, :), y_aux(end, :), 100);
                    x(:, i+1) = y_aux(end, :) + z(i+1, :);

                    % Simulates the linear model using the discrete differences equation
                    %       x_{k+1} = A x_k + B u_k + L ( Y_k - C x_k )
                    x_hat(:, i+1) = A_d*x_hat(:, i) + B_d*u(:,i) + L * (C_d * x(:, i) - C_d * ( x_hat(1:nx, i) + xe ));
                    
                end
                
            % Case for a Infinite-Horizon Continuous-Time
            else
                % Augments the state auxiliary vector in time and updates the initial conditions
                x_hat = zeros(nx+ny, numel(t));
                x_hat(:,1) = [X_0' - xe; r(:,1) - C*(X_0' - xe) ];
            
                % Augment the state and input matrices
                A_i = [A, zeros(nx, ny); 
                          -C, zeros(ny)];
                B_i = [B; 
                           zeros(ny, nu)];
                C_i = [C, zeros(ny, ny)];
                D_i = D;
                
                % Calculate the LQRI gain matrix
                K = lqr_(A_i, B_i, Q, R, N);
                
                % Calculate the Closed-Loop State-Space Representation
                A_cl = A_i - B_i * K - L_i * C_i;
                B_cl = [zeros(nx, ny); eye(ny)];
                
                % == Simulation for each timestep (i) on time signal (t) ==
                for i = 1:numel(t)-1
                    % Calculates the input signal
                    u(:,i) = - K * x_hat(:, i);

                    % Simulates the real physical plant using the non-linear model
                    [~, y_aux] = odeSolver(model.model, t(i:i+1), u(:,i) + ue + w(i, :), y_aux(end, :), 100);
                    x(:, i+1) = y_aux(end, :) + z(i+1, :);

                    % Simulates the linear model using the Lagrange formula
                    %       x_k = e^{A_i (t - t_0)} x_0 + \int{ e^{A_i (t - \tau)} B_i u_\tau + F_i r_\tau + L_i ( Y_\tau - C_i x_\tau ) }
                    n_resp = expm( A_cl * (t(i+1) - t(1)) ) * x_hat(:,1);
                    f_resp = 0;
                    for j =1:i+1
                        f_resp = f_resp + deltaX * expm( A_cl * (t(i+1) - t(j)) ) * (B_cl * r(:,j) +  L_i *(C *  (x(:,j) - xe))  ) ;
                    end
                    
                    x_hat(:, i+1) = n_resp + f_resp;
                    
                end
            end
            
        % == LINEAR QUADRATIC GAUSSIAN REGULATOR ==,
        elseif(strcmp(type, 'lqg'))
            % Estimate the covariances for the process and measurement noises
            Q_k = cov(w);
            R_k = cov(z);
            P_k = zeros(2);
            
            % Case for the Finite-Horizon Discrete-Time
            if(N ~= 'inf') 
                % Defines the discrete matrices associated to state-space matrices
                ss_d = c2d(ss(A, B, C, D), deltaX);
                A_d = ss_d.A; B_d = ss_d.B; C_d = ss_d.C; D_d = ss_d.D;       
                
                K = lqr_(A_d, B_d, Q, R, N);

                % == Simulation for each timestep (i) on time signal (t) ==
                for i = 1:numel(t)-1
                    % Calculates the input signal (if still on control horizon)
                    if(i <= N)
                        u(:,i) = K(i, :) * ( r(:,i) - C_d *  x_hat(:, i) );
                    end

                    % Simulates the real physical plant using the non-linear model
                    [~, y_aux] = odeSolver(model.model, t(i:i+1), u(:,i) + ue + w(i, :), y_aux(end, :), 100);
                    x(:, i+1) = y_aux(end, :) + z(i+1, :);
                    
                    % Kalman Filter's Prediction
                    x_hat(:, i+1) = A_d*x_hat(:, i) + B_d*u(:,i);
                    P_k = A_d * P_k * A_d' + Q_k;
                    
                    % Kalman Filter's Update
                    S_k = C_d * P_k * C_d' + R_k;
                    K_k = P_k * C_d' / S_k;
                    err = C_d * x(:,i+1) - C_d *  (x_hat(:, i+1)+ xe);
                    
                    x_hat(:, i+1) = x_hat(:, i+1) + K_k * err;
                    P_k = P_k - K_k * S_k * K_k';
                    
                end
                
            % Case for the Infinite-Horizon Continuous-Time
            else
                % Calculate the LQR gain matrix
                K = lqr_(A, B, Q, R, N);  
            
                % Generates the Kalman observer state-space
                [~, L] = kalman(ss(A, B, C, D), Q_k, R_k); 
                
                % Define the closed-loop matrices
                A_cl = A - B*K - L*C;
                B_cl = B*K;
                
                % == Simulation for each timestep (i) on time signal (t) ==
                for i = 1:numel(t)-1
                    % Calculates the input signal
                    u(:,i) = K * ( r(:,i) - C * x_hat(:, i) );

                    % Simulates the real physical plant using the non-linear model
                    [~, y_aux] = odeSolver(model.model, t(i:i+1), u(:,i) + ue + w(i, :), y_aux(end, :), 100);
                    x(:, i+1) = y_aux(end, :) + z(i+1, :);

                    % Simulates the linear model using the Lagrange formula
                    %       x_k = e^{A (t - t_0)} x_0 + \int{ e^{A (t - \tau)} B u_\tau + L( Y_\tau - C x_\tau ) }
                    n_resp = expm( A_cl * (t(i+1) - t(1)) ) * x_hat(:,1);
                    f_resp = 0;
                    for j =1:i+1
                        f_resp = f_resp + deltaX * expm( A_cl * (t(i+1) - t(j)) ) * (B_cl * r(:,j) + L * ( x(:,j) - xe) );
                    end

                    x_hat(:, i+1) = n_resp + f_resp;
                    
                end
            end
           
        end

        tout = t;
        uout = u;
        xout = x;
        yout = x_hat(1:nx, :);
    end
   
end

