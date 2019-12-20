function [c,f,s] = pdefun(x,~,u,dudx,params)
    J_B = params(1); %(s^-1)
    J_R = params(2); %(s^-1)
    mu = params(3); %(s^-1)
    mu_hb = params(4); %(s^-1)
    mu_Hb = params(5); %(s^-1)
    alpha = params(6); %(s^-1)
    beta = params(7); %(um^-1 s^-1)
    K = params(8); %(um^-1)
    x_B = params(9); %(um)
    x_R = params(10); %(um)
    nu = params(11); %(um s^-1)
    eta = params(12);
    D = params(13); %(um^2 s^-1)
    D_hb = params(14); %(um^2 s^-1)
    D_Hb = params(15); %(um^2 s^-1)
    gamma = params(16);
    frac = params(17);
    delta_radius = params(18); %(um)
    
    c = ones(5,1);
    f = [D; D; D; D_hb; D_Hb] .* dudx;

    s = zeros(5,1);
    s(1) = -mu*u(1)-nu*u(1)*u(2)+J_B*(abs(x-x_B) < delta_radius)/(2*delta_radius);
    s(2) = -mu*u(2)-nu*u(1)*u(2)+frac*J_R*(abs(x-x_R) < delta_radius)/(2*delta_radius)...
            +(1-frac)*J_R/max(x);
    s(3) = -mu*u(3)+nu*u(1)*u(2);
    s(4) = -mu_hb*u(4)+(beta*u(1)^3*(eta*u(5)+gamma*K))/(K^4+u(1)^3*(eta*u(5)+K));
    s(5) = -mu_Hb*u(5)+alpha*u(4);
end