function [ tout, yout ] = simulate( model, t, U, X_0 )
%SIMULATE Simulates a non-linear model using the Euler method.
%   @param model   The non-linear ODEs representing the system model.
%   @param X_0       The initial conditions of the simulation.
%   @param U          The input signal for the simulation.

    tout = t;
    yout = zeros(size(t, 1), size(X_0, 2));
    yout(1,:) = X_0;
    for i = 2:size(t,1)
        [~, y_aux] = odeSolver(model, t(i-1:i), U(i-1,:), yout(i-1,:), 100); 
        yout(i,:) = y_aux(end,:);
    end

end

