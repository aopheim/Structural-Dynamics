function [p, num_t_entries, t] = createForceAndTime(h, omega, F_max, deltaT)

T = 1 / 2 * (1 / omega);       % The period of the force
num_t_entries = int64((T + deltaT) / h);

t = zeros(num_t_entries, 1);
p = zeros(num_t_entries, 1);

for i = 1 : 1 : size(t,1)
    if i ~= num_t_entries
        t(i + 1,1) = t(i,1) + h;
    end
    
    if (t(i,1) <= (T / 2))
        p(i,1) = F_max / (T/2) * t(i, 1);
    elseif ((t(i,1) > T/2) && (t(i,1) <= T ))
        p(i,1) = - F_max / (T/2) * t(i, 1) + 2 * F_max;
    else
        p(i,1) = 0;
    end
end






