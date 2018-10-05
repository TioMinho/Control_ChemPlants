function [ y_out, x_out, u_out, t_out ] = lqr_sim( A, B, C, D, Q, R, T, X0, t )
%LQR_FH_SIM Simulates the system with a Discrete-Time Finite-Horizon LQR Controller 
%   Detailed explanation goes here

    t_out = t'; x_out = zeros(size(A,2), numel(t)); u_out = zeros(size(B,2), numel(t)); X0 = X0';
    deltaX = t_out(2) - t_out(1); 
    
    ss_d = c2d(ss(A,B,C,D), deltaX); 
    A = ss_d.A; B = ss_d.B; C = ss_d.C; D = ss_d.D;
    
    P = zeros(size(Q, 1), size(Q, 2), numel(t)); P(:,:,numel(t)) = Q;
    for i = T:-1:2
        P(:,:,i-1) = Q + A' * P(:,:,i) * A - A' * P(:,:,i) * B * pinv(R + B' * P(:,:,i) * B) * B' * P(:,:,i) * A;
    end
    
    x_out(:, 1) = X0;
    for i = 1:1:numel(t_out)-1
        if(i < T)
        u_out(:,i) = - pinv(R + B' * P(:,:,i+1) * B) * B' * P(:,:,i+1) * A * x_out(:,i);
        end
        
        x_out(:, i+1) =  A * x_out(:,i) + B * u_out(:, i);
    end  

    y_out = C * x_out + D * u_out;
    
    y_out = y_out'; x_out = x_out'; u_out = u_out'; 
end

% function [ y_out, x_out, u_out, t_out ] = lqr_sim( A, B, C, D, Q, R, T, X0, t )
% %LQR_FH_SIM Simulates the system with a Discrete-Time Finite-Horizon LQR Controller 
% %   Detailed explanation goes here
% 
%     t_out = t'; x_out = zeros(size(A,2), numel(t)); u_out = zeros(size(B,2), numel(t)); X0 = X0';
%     deltaX = t_out(2) - t_out(1); x_out(:, 1) = X0;
%     
%     ss_d = c2d(ss(A,B,C,D), deltaX); 
%     A = ss_d.A; B = ss_d.B; C = ss_d.C; D = ss_d.D;
%     
%     for k = 1:numel(t)-1
%         P = zeros(size(Q, 1), size(Q, 2), numel(t)); P(:,:,numel(t)) = Q;
%         for i = T-1:-1:k+1
%             P(:,:,i-1) = Q + A' * P(:,:,i) * A - A' * P(:,:,i) * B * pinv(R + B' * P(:,:,i) * B) * B' * P(:,:,i) * A;
%         end
%         
%         u_out(:,k) = - pinv(R + B' * P(:,:,2) * B) * B' * P(:,:,2) * A * x_out(:,k);
%         x_out(:, k+1) =  A * x_out(:,k) + B * u_out(:, k);
%     end
%     y_out = C * x_out + D * u_out;
%     
%     y_out = y_out'; x_out = x_out'; u_out = u_out'; 
% end


