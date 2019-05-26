function exp_param = generateRandomExperiment(xe)
% generateRandomExperiment
%   Detailed explanation goes here    
    types = {'lqg', 'lqgi'};
    Q_values = [1 5:5:500];
    Qi_values = 12:1:15;
    R_values = [1 5:5:500];
    w_values = [linspace(0.1, 2); linspace(0.1, 2); linspace(0.01, 1); linspace(0.01, 1)]';
    z_values = [linspace(0.1, 2); linspace(0.01, 1)]';
    r_values_track = [xe(2)*linspace(0.5, 1.5); xe(3)*linspace(0.9, 1.1)];
    x0_values = [xe(1)*linspace(0.1, 1.1); xe(2)*linspace(0.1, 1.1); xe(3)*linspace(0.99, 1.01); xe(4)*linspace(0.99, 1.01)];
    
    R = diag(R_values(randperm(numel(R_values), 2)));
    
    type = types(randperm(numel(types),1));
    if(strcmp(type{1}(end), 'i'))
        Qi = Qi_values(randperm(numel(Qi_values), 2));
        Q = diag([Q_values(randperm(numel(Q_values), 4)) 10^Qi(1) 10^Qi(2)]);
    
        t = 0:0.1:23.9; T = numel(t);
        
        idx = randperm(size(r_values_track, 2), 4);
        rr = r_values_track;
        r = [ones(2,T/4).*[xe(2); xe(3)] ones(2,T/4).*[rr(1,idx(1)); rr(2,idx(2))] ones(2,T/4).*[rr(1,idx(3)); rr(2,idx(4))] ones(2,T/4).*[xe(2); xe(3)]];
        x0 = xe;
    else
        Q = diag(Q_values(randperm(numel(Q_values), 4)));
        
        t = 0:0.0025:0.5975; T = numel(t);
        r = ones(2,T).*[xe(2); xe(3)];
        
        idx = randperm(size(x0_values,2), 4);
        x0 = [x0_values(1,idx(1)) x0_values(2,idx(2)) x0_values(3,idx(3)) x0_values(4,idx(4))];
    end
    
    idx = randperm(size(w_values,1), 4);
    w = randn(T, 4) .* [w_values(idx(1), 1) w_values(idx(2), 2) w_values(idx(3), 3) w_values(idx(4), 4)];
    idx = randperm(size(z_values,1), 2);
    z = randn(T, 2) .* [z_values(idx(1), 1) z_values(idx(2), 2)];
    
    exp_param.w = w; exp_param.z = z;
    exp_param.t = t'; exp_param.r = r; exp_param.x_0 = x0;
    exp_param.type = type{1}; exp_param.Q = Q; exp_param.R = R; 
    exp_param.xout = []; exp_param.yout = []; exp_param.uout = [];
    
end