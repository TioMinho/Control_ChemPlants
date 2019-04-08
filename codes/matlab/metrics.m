function metrics = metrics( y_hat, y_real, t)
%METRICS Calculate the performance indices for a controlled simulation.
%   Detailed explanation goes here

    ISE   = sum((y_hat - y_real) .^ 2, 2);
    IAE   = sum(abs(y_hat - y_real), 2);
    ITSE = sum((y_hat - y_real) .^ 2 .* t', 2);
    ITAE = sum(abs(y_hat - y_real) .* t', 2);

    metrics = [ISE, IAE, ITSE, ITAE];
    
end

