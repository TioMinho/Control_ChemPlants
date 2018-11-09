function [ K, P ] = lqr_( A, B, Q, R, N )
%LQR_FH_SIM Simulates the system with a Discrete-Time Finite-Horizon LQR Controller 
%   Detailed explanation goes here    
    if(N ~= 'inf')
        P = zeros(size(Q, 1), size(Q, 2), N); P(:,:,N) = Q;
        K = zeros(N-1, size(A,1));

        for i = N:-1:2
            P(:,:,i-1) = Q + A' * P(:,:,i) * A - A' * P(:,:,i) * B * pinv(R + B' * P(:,:,i) * B) * B' * P(:,:,i) * A;
            K(i-1, :) = - pinv(R + B' * P(:,:,i) * B) * B' * P(:,:,i) * A;
        end

    elseif(N == 'inf')
        [K, P, ~] = lqr(A, B, Q, R);
    end
    
end