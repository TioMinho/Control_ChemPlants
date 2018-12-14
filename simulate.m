function [ tout, yout, xout, uout ] = simulate( varargin )
%SIMULATE Simulates a non-linear model using the Euler method.

    switch length(varargin)
        case 4
        %   @param model   The non-linear ODEs representing the system model.
        %   @param t       The initial conditions of the simulation.
        %   @param U       The input signal for the simulation.        
        %   @param X_0     The initial conditions of the simulation.
        model = varargin{1}.model; t = varargin{2}; U = varargin{3}; X_0 = varargin{4};
        
        tout = t; uout = U;
        yout = zeros(numel(X_0), numel(t));
        yout(:,1) = X_0;
        for i = 1:numel(t)
            [~, y_aux] = odeSolver(model, t(i:i+1), U(:,i), yout(:,i), 100); 
            yout(:,i+1) = y_aux(:,end);
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
        
        case num2cell(8:1:9)
        %   @param model     The non-linear ODEs representing the system model.
        %   @param idx       The index position of the linear state-space model to use.
        %   @param t         The initial conditions of the simulation.
        %   @param r         The input signal for the simulation.        
        %   @param X_0       The initial conditions of the simulation.
        %   @param type      The type of controller being used to simulate.
        %   @param K         The state feedback gain vector being used.
        %   @param L         The Luenberger gain vector being used.
        %   @param w         The matrix of additive disturbance in the nonlinear model.
        model = varargin{1}; idx = varargin{2}; t = varargin{3}; r = varargin{4}; X_0 = varargin{5}; 
        type = varargin{6}; K = varargin{7}; L = varargin{8};

        A = model.ss_model.A(idx); B = model.ss_model.B(idx); C = model.ss_model.C; D = model.ss_model.D;
        
        if(length(varargin) == 9)
            w = varargin{9};
        else
            w = zeros(1, size(A, 1));
        end
        
        deltaX = t(2) - t(1);
        
        ss_d = c2d(ss(A,B,C,D), deltaX);
        A_d = ss_d.A; B_d = ss_d.B; C_d = ss_d.C; D_d = ss_d.D;
        
        u     = zeros(size(B,2), numel(t)); 
        y     = zeros(numel(t), size(A,1));
        x     = zeros(size(A,1), numel(t)); 
        x_hat = zeros(size(A,1), numel(t));
        
        y(1,:) = X_0; x(:,1) = X_0; x_hat(:,1) = X_0 - model.oper.X(idx,:);
        r = r - model.oper.X(idx,:)';
        
        if(type == 'lqr')
            % Case for a Finite-Horizon Discrete-Time Linear Quadratic Regulator
            if(size(K, 1) > 1) 
                N = size(K,1);

                for i = 1:numel(t)-1
                    if(i <= N)
                        u(:,i) = K(i,:) * ( r(:,i) - x_hat(:, i) );
                    end

                    [~, y_aux] = odeSolver(model.model, t(i:i+1), u(:,i)+model.oper.U(idx,:), y(i,:), 100);
                    y(i+1,:) = y_aux(end, :) + w(i, :);

                    x(:,i+1) = y(i+1, :);

                    x_hat(:, i+1) = A_d*x_hat(:, i) + B_d*u(:,i) + L*( x(:, i) - C_d*(x_hat(:, i) + model.oper.X(idx,:)') );
                end
                
            % Case for a Infinite-Horizon Continuous-Time Linear Quadratic Regulator
            else                
                for i = 1:numel(t)-1
                    u(:,i) = K * ( r(:,i) - (x_hat(:, i)) );

                    [~, y_aux] = odeSolver(model.model, t(i:i+1), u(:,i)+model.oper.U(idx,:), y(i,:), 100);
                    y(i+1,:) = y_aux(end, :) + w(i, :);

                    x(:,i+1) = y(i+1, :);

                    n_resp = expm( A * (t(i) - t(1)) ) * x_hat(:,1);
                    f_resp = 0;
                    for j =1:i+1
                        f_resp = f_resp + deltaX * expm( A * (t(i+1) - t(j)) ) * (B * u(:,j) + L *(x(:,j) - C * (x_hat(:,j) + model.oper.X(idx,:)')) );
                    end
                    
                    x_hat(:, i+1) = n_resp + f_resp;
                end
            end
        end

        tout = t;
        uout = u;
        xout = x;
        yout = C * x_hat + D * u;

        case num2cell(10:1:11)
        %   @param model     The non-linear ODEs representing the system model.
        %   @param idx       The index position of the linear state-space model to use.
        %   @param t         The initial conditions of the simulation.
        %   @param r         The input signal for the simulation.        
        %   @param X_0       The initial conditions of the simulation.
        %   @param type      The type of controller being used to simulate.
        %   @param L         The Luenberger gain vector being used.
        %   @param w         The matrix of additive disturbance in the nonlinear model.
        %   @param Q         The weighting matrix to the states signals.
        %   @param R         The weighting matrix to the input signals.
        %   @param N         The size of the control time horizon.
        
        model = varargin{1}; idx = varargin{2}; t = varargin{3}; r = varargin{4}; X_0 = varargin{5}; 
        type = varargin{6}; L = varargin{7}; 

        A = model.ss_model.A(idx); B = model.ss_model.B(idx); C = model.ss_model.C; D = model.ss_model.D;
        
        if(length(varargin) == 8)
            w = varargin{8};
        else
            w = zeros(numel(t), size(A, 1));
        end
        
        deltaX = t(2) - t(1);
              
        if(type == 'lqri')
            % Augment the state and input matrices
            A_i = [A, zeros(size(C')); -C, zeros(size(C,1), size(C,1))];
            B_i = [B; zeros(size(C,1), size(B,2))];
            L_i = [L; zeros(size(C,1), size(L,2))];
            F_i = [zeros(size(A,1), size(r,1)); eye(size(r,1))];

            % Calculate the new LQR gain matrix
            Q = varargin{9}; R = varargin{10}; N = varargin{11};
            K = lqr_(A_i, B_i, Q, R, N);

            % Defines the discrete matrices associated to augmented matrices
            % ss_d = c2d(ss(A_i, B_i, C, D), deltaX);
            % A_d = ss_d.A; B_d = ss_d.B; C_d = ss_d.C; D_d = ss_d.D;
            
            u     = zeros(size(B,2), numel(t)); 
            y     = zeros(numel(t), size(A,1));
            x     = zeros(size(A,1), numel(t)); 
            x_hat = zeros(size(A_i,1), numel(t));
                    
            y(1,:) = X_0; x(:,1) = X_0; x_hat(:,1) = [X_0' - model.oper.X(idx,:)'; X_0' - r(:,1)];
            r = r - model.oper.X(idx,:)';

            % Case for a Finite-Horizon Discrete-Time Linear Quadratic Regulator
            if(N ~= 'inf') 
                for i = 1:numel(t)-1
                    if(i <= N)
                        u(:,i) = - K(i,:) * x_hat(:, i);
                    end

                    [~, y_aux] = odeSolver(model.model, t(i:i+1), u(:,i)+model.oper.U(idx,:), y(i,:), 100);
                    y(i+1,:) = y_aux(end, :) + w(i, :);

                    x(:,i+1) = y(i+1, 1:size(A,1));

                    x_hat(:, i+1) = (A_d*x_hat(:, i) + B_d*u(:,i) + F_i*r(:,i)) + L_i*( x(:, i) - C_d*(1:x_hat(size(x,1), i) + model.oper.X(idx,:)') );
                end
                
            % Case for a Infinite-Horizon Continuous-Time Linear Quadratic Regulator
            else                
                for i = 1:numel(t)-1
                    u(:,i) = - K * x_hat(:, i);

                    [~, y_aux] = odeSolver(model.model, t(i:i+1), u(:,i)+model.oper.U(idx,:), y(i,:), 100);
                    y(i+1,:) = y_aux(end, :) + w(i, :);

                    x(:,i+1) = y(i+1, :);

                    n_resp = expm( A_i * (t(i) - t(1)) ) * x_hat(:,1);
                    f_resp = 0;
                    for j =1:i+1
                        f_resp = f_resp + deltaX * expm( A_i * (t(i) - t(j)) ) * (B_i * u(:,j) + F_i*r(:,i) + L_i *(x(:,j) - C * (x_hat(1:size(x,1), j) + model.oper.X(idx,:)')) );
                    end
                    
                    x_hat(:, i+1) = n_resp + f_resp;
                end
            end
        end

        tout = t;
        uout = u;
        xout = x;
        yout = C * x_hat(1:size(C,2), :) + D * u;
    end
   
end

