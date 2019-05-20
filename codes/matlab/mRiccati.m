function dXdt = mRiccati(t, P, A, B, Q, R)
%MRICCATI The matrix Riccati differential equation
    P = reshape(P, size(A)); 
    dXdt = -(A.'*P + P*A - P*B*pinv(R)*B.'*P + Q); 
    dXdt = dXdt(:); 
end

