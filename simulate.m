function [ tout, yout, xout, uout ] = simulate( varargin )
%SIMULATE Simulates a non-linear model using the Euler method.

    switch varargin
        case 4
        %   @param model   The non-linear ODEs representing the system model.
        %   @param t       The initial conditions of the simulation.
        %   @param U       The input signal for the simulation.        
        %   @param X_0     The initial conditions of the simulation.
        model = varargin{1}.model; t = varargin{2}; U = varargin{3}; X_0 = varargin{4};
        
        tout = t; uout = U;
        yout = zeros(size(t, 1), size(X_0, 2));
        yout(1,:) = X_0;
        for i = 1:numel(t)
            [~, y_aux] = odeSolver(model, t(i:i+1), U(i,:), yout(i,:), 100); 
            yout(i+1,:) = y_aux(end,:);
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

        a = 'TODO: nonlinear model control'
            
        
        case num2cell(7:1:8)
        %   @param model     The non-linear ODEs representing the system model.
        %   @param idx       The index position of the linear state-space model to use.
        %   @param t         The initial conditions of the simulation.
        %   @param r         The input signal for the simulation.        
        %   @param X_0       The initial conditions of the simulation.
        %   @param K         The state feedback gain vector being used.
        %   @param L         The Luenberger gain vector being used.
        model = varargin{1}; idx = varargin{2}; t = varargin{3}; r = varargin{4}; X_0 = varargin{5}; 
        K = varargin{6}; L = varargin{7}; w = 0;

        if(varargin == 8)
            w = varargin{8};
        end
        
        A = model.ss_model.A(idx); B = model.ss_model.B(idx); C = model.ss_model.C; D = model.ss_model.D;
        
        u = zeros(size(B, 1), numel(t)); 
        x = zeros(size(A,1), numel(t)); x_hat = zeros(size(A,1), numel(t));

        x(:,1) = X_0; x_hat(:,1) = model.oper.X(idx,:);
        if(size(K, 1) > 1) % Case for a Finite-Horizon Discrete-Time Linear Quadratic Regulator
            N = size(K,1);

            for i = 1:numel(t)-1
                if(i <= N)
                    u(:,i) = K(i,:)*( r(:,i) - (x_hat(:, i) + model.oper.X(idx,:)) );
                end

                [~, y_aux] = odeSolver(model, t(i:i+1), u(i,:)+model.oper.U(idx,:), x(i,:), 100);

                x(:,i+1) = y_aux(:, 1:size(A,1)) + w(:, i+1);
                x_hat(:, i+1) = A*x_hat(:, i) + B*u(:, i) + L*( x(:, i) - C*(x_hat(:, i) + model.oper.X(idx,:)) ); 
            end
        else                % Case for a Infinite-Horizon Continuous-Time Linear Quadratic Regulator
            for i = 1:numel(t)-1
                u(:,i) = K * ( r(:,i) - (x_hat(:, i) + model.oper.X(idx,:)) );
                
                [~, y_aux] = odeSolver(model, t(i:i+1), u(i,:)+model.oper.U(idx,:), x(i,:), 100);

                x(:,i+1) = y_aux(:, 1:size(A,1));
                x_hat(:, i+1) = A*x_hat(:, i) + B*u(:, i) + L*( x(:, i) - C*(x_hat(:, i) + model.oper.X(idx,:)) ); 
            end
        end

        tout = t;
        uout = u;
        xout = x;
        yout = C * x_hat + D * u;
    end
   
end

