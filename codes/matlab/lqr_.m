function [ K, P ] = lqr_( A, B, Q, R, N )
%LQR_ Linear Quadratic Regulator for Continuous-Time and Finite-Horizon
%optimal controllers
%
%   [K,P] = LQR_(A,B,Q,R,N) Solves the LQR finite-horizon problem defined by
%
%         \int_{0}^{N} (x(t)'Qx(t) + u(t)'Ru(t) dt) + x(N)'Qx(N) ,
%
%   given the resulting optimal control law K(t) = - R^{-1} B' P(t), where
%   P(t) is the solution of the matrix Riccati differential equation:
%   
%         (dP(t)/dt) = A' P(t) + P(t) A - P(t) B R^{-1} B' P(t) + Q
%
    Q_f = reshape(Q, size(A));
    [~, P] = ode45(@(t,P) mRiccati(t, P, A, B, Q, R), fliplr(N), Q_f);

    [m, n] = size(P);
    P = mat2cell(P, ones(m,1), n);
    fh_reshape = @(p) reshape(p, size(A));
    fh_getK    = @(p) pinv(R) * B' * p;

    P = cellfun(fh_reshape, P, 'UniformOutput', false);
    K = cellfun(fh_getK, P, 'UniformOutput', false);

end