function [ t, y ] = odeSolver( odefun, tspan, u, y0, N )
% ODESOLVER Solving ordinary differential equations via the Euler method.
    h = (tspan(2) - tspan(1)) / (N - 1);
    t = linspace(tspan(1), tspan(2), N);
    y = zeros(N, numel(y0));
    
    y(1,:) = y0;
    for i = 1:N-1
        y(i+1, :) = y(i,:) + h * odefun(t(i), u, y(i,:));
    end
end

