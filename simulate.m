function [ tout, yout, xout, uout ] = simulate( varargin )
%SIMULATE Simulates a non-linear model using the Euler method.

    switch varargin
        case 4
        %   @param model   The non-linear ODEs representing the system model.
        %   @param t       The initial conditions of the simulation.
        %   @param U       The input signal for the simulation.        
        %   @param X_0     The initial conditions of the simulation.
        model = varargin{1}; t = varargin{2}; U = varargin{3}; X_0 = varargin{4};
        
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
            
        
        case 10
        %   @param model     The non-linear ODEs representing the system model.
        %   @param A,B,C,D   The matrices representing the linear state-space model.
        %   @param t         The initial conditions of the simulation.
        %   @param r         The input signal for the simulation.        
        %   @param X_0       The initial conditions of the simulation.
        %   @param K         The state feedback gain vector being used.
        %   @param L         The Luenberger gain vector being used.
        model = varargin{1}; A = varargin{2}; B = varargin{3}; C = varargin{4}; D = varargin{3};
        t = varargin{6}; r = varargin{7}; X_0 = varargin{8}; K = varargin{9}; L = varargin{10};
            
        u = zeros(size(B, 1), numel(t)); 
        x = zeros(size(A,1), numel(t)); x_hat = zeros(size(A,1), numel(t));

        x(:,1) = X_0 + (rand()-0.5)*4; x_hat(:,0) = X_0;
        if(size(K, 1) > 1) % Case for a Finite-time Discrete Linear Quadratic Regulator
            N = size(K,1);

            for i = 1:numel(t)-1
                if(i <= N)
                    u(:,i) = (r(:,i) - K(i,:)*x_hat(:,i));
                end

                [~, y_aux] = odeSolver(model, t(i:i+1), u(i,:), x(i,:), 100);

                x(:,i+1) = y_aux;
                x_hat(:, i+1) = (A - L * C)  * x_hat(:, i) + B * u(:, i) + L * C * x(:, i); 
            end
        else
            for i = 1:numel(t)-1
                u(:,i) = (r(:,i) - K*x_hat(:,i));
                
                [~, y_aux] = odeSolver(model, t(i:i+1), u(i,:), x(i,:), 100);

                x(:,i+1) = y_aux;
                x_hat(:, i+1) = (A - L * C)  * x_hat(:, i) + B * u(:, i) + L * C * x(:, i); 
            end
        end

        tout = t;
        uout = u;
        xout = x;
        yout = C * x_hat + D * u;
    end
   
end

