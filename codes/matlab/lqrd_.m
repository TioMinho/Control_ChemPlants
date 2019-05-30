function [ K, P ] = lqrd_( A, B, Q, R, N )
%LQRD_ Linear Quadratic Regulator for Discrete-Time and Finite-Horizon
%optimal controllers
%
%   [K,P] = LQRD_(A,B,Q,R,N) Solves the LQR finite-horizon problem defined by
%
%         \sum_{n=0}^{N} (x[n]'Qx[n] + u[n]'Ru[n]) + x[N]'Qx[N] ,
%
%   given the resulting optimal control law K(t) = - (R + B'P[k+1]B)^{-1} B' P[k+1] A, where
%   P[k+1] is the solution of the matrix Riccati differences equation:
%   
%         P[k] = A' P[k+1] A - A' P[k+1] B (R + B' P[k+1] B)^{-1} B' P[k+1] A
%
    P = zeros(size(Q, 1), size(Q, 2), N); P(:,:,N) = Q;
    K = cell(N-1);
    
    for i = N:-1:2
        P(:, :, i-1) = Q + A' * P(:, :, i) * A - A' * P(:, :, i) * B * pinv(R + B' * P(:, :, i) * B) * B' * P(:, :, i) * A;
        K{i-1} = pinv(R + B' * P(:,:,i) * B) * B' * P(:,:,i) * A;
    end
end

