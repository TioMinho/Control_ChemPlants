function dXdt = mRiccati(t, P, A, B, Q, R)
%MRICCATI The matrix Riccati differential equation.
%
%   dXdt = MRICCATI(P,A,B,Q,R) Computes the matrix Riccati differential
%   equation formulated as:
%
%       (dP(t)/dt) = A' P(t) + P(t) A - P(t) B R^{-1} B' P(t) + Q
%

    P = reshape(P, size(A)); 
    dXdt = -(A.'*P + P*A - P*B*pinv(R)*B.'*P + Q); 
    dXdt = dXdt(:); 

end

