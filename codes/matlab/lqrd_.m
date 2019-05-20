function [ K, P ] = lqrd_( A, B, Q, R, N )
%LQRD_ Summary of this function goes here
%   Detailed explanation goes here
    P = zeros(size(Q, 1), size(Q, 2), N); P(:,:,N) = Q;
    K = cell(N-1);
    
    for i = N:-1:2
        P(:, :, i-1) = Q + A' * P(:, :, i) * A - A' * P(:, :, i) * B * pinv(R + B' * P(:, :, i) * B) * B' * P(:, :, i) * A;
        K{i-1} = pinv(R + B' * P(:,:,i) * B) * B' * P(:,:,i) * A;
    end
end

