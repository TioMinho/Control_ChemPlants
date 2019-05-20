function [ Ke, Pe ] = kalmanBucy( A, C, Q, R, N )
%LQR_FH_SIM Simulates the system with a Discrete-Time Finite-Horizon LQR Controller 
%   Detailed explanation goes here    
    Q_i = reshape(Q, size(A));
    [~, Pe] = ode45(@(t,P) -mRiccati(t, P, A', C', Q, R), N, Q_i);
    
    [m, n] = size(Pe);
    Pe = mat2cell(Pe, ones(m,1), n);
    fh_reshape = @(p) reshape(p, size(A));
    fh_getK    = @(p) p * C' * pinv(R);

    Pe = cellfun(fh_reshape, Pe, 'UniformOutput', false);
    Ke = cellfun(fh_getK, Pe, 'UniformOutput', false);

end