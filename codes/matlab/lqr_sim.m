function [ y_out, x_out, u_out, t_out ] = lqr_sim( A, B, C, D, Q, R, X0, t )
%LQR_SIM Simulates the system with a Continuous-Time Infinite-Horizon LQR Controller 
%   Detailed explanation goes here

    t_out = t'; x_out = zeros(size(A,2), numel(t)); x_hat = zeros(size(A,2), numel(t)); u_out = zeros(size(B,2), numel(t)); X0 = X0';
    
    [K, ~, ~] = lqr(A, B, Q, R);
    
    deltaX = t_out(2) - t_out(1);
    x_out(:, 1) = X0;
    for i = 2:1:numel(t_out)
        u_out(:, i-1) = -K * x_out(:,i -1);
        
        n_resp = expm( A * (t_out(i) - t_out(1)) ) * X0;
        f_resp = 0;
        for j = 1:i
            f_resp = f_resp + deltaX * expm( A * (t_out(i) - t_out(j)) ) * B * u_out(:,j);
        end
        
        x_out(:, i) =  n_resp + f_resp;
    end  

    y_out = C * x_out + D * u_out;
    
    y_out = y_out'; x_out = x_out'; u_out = u_out'; 
end

