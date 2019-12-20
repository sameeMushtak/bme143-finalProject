function [pl,ql,pr,qr] = pdebc(~,~,~,~,~)
% Zero flux boundary conditions on both ends of the embryo
    pl = zeros(5,1);
    ql = ones(5,1);
    pr = zeros(5,1);
    qr = ones(5,1);
end