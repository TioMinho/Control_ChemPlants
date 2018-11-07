function [ tout, yout ] = simulate( varargin )
%SIMULATE Simulates a non-linear model using the Euler method.

    switch varargin
        case 4
        %   @param model   The non-linear ODEs representing the system model.
        %   @param t       The initial conditions of the simulation.
        %   @param U          The input signal for the simulation.        
        %   @param X_0       The initial conditions of the simulation.
        model = varargin{1}; t = varargin{2}; U = varargin{3}; X_0 = varargin{4};
        
        tout = t;
        yout = zeros(size(t, 1), size(X_0, 2));
        yout(1,:) = X_0;
        for i = 2:size(t,1)
            [~, y_aux] = odeSolver(model, t(i-1:i), U(i-1,:), yout(i-1,:), 100); 
            yout(i,:) = y_aux(end,:);
        end
        
        case 5
        %   @param model   The non-linear ODEs representing the system model.
        %   @param t       The initial conditions of the simulation.
        %   @param U       The input signal for the simulation.        
        %   @param X_0     The initial conditions of the simulation.
        %   @param ctrl    The type of feedback controller that is being used.
        model = varargin{1}; t = varargin{2}; U = varargin{3}; X_0 = varargin{4};
        ctrl = varargin{5};
            
        case 7
        %   @param A,B,C,D   The matrices representing the linear state-space model.
        %   @param t         The initial conditions of the simulation.
        %   @param U         The input signal for the simulation.        
        %   @param X_0       The initial conditions of the simulation.
        model = varargin{1}; t = varargin{2}; U = varargin{3}; X_0 = varargin{4};
        ctrl = varargin{5};
            
        case 8
        %   @param A,B,C,D   The matrices representing the linear state-space model.
        %   @param t         The initial conditions of the simulation.
        %   @param U         The input signal for the simulation.        
        %   @param X_0       The initial conditions of the simulation.
        %   @param ctrl      The type of feedback controller that is being used.
        model = varargin{1}; t = varargin{2}; U = varargin{3}; X_0 = varargin{4};
        ctrl = varargin{5};
    end
   
end

