function [ K, P ] = lqr_( A, B, Q, R, N )
%LQR_FH_SIM Simulates the system with a Discrete-Time Finite-Horizon LQR Controller 
%   Detailed explanation goes here    
    Q_f = reshape(Q, size(A));
    [~, P] = ode45(@(t,P) mRiccati(t, P, A, B, Q, R), fliplr(N), Q_f);

    [m, n] = size(P);
    P = mat2cell(P, ones(m,1), n);
    fh_reshape = @(p) reshape(p, size(A));
    fh_getK    = @(p) pinv(R) * B' * p;

    P = cellfun(fh_reshape, P, 'UniformOutput', false);
    K = cellfun(fh_getK, P, 'UniformOutput', false);

end