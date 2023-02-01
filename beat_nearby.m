function [minf, point] = beat_nearby(delta, point, prebest, nvars, objective_function)
minf = prebest;
z = point;
for i = 1 : nvars
    z(i) = point(i) + delta(i);
    ftmp = objective_function(z);
    if ftmp < minf
        minf = ftmp;
    else
        delta(i) = 0 - delta(i);
        z(i) = point(i) + delta(i);
        ftmp = objective_function(z);
        if ftmp < minf
            minf = ftmp;
        else
            z(i) = point(i);
        end
    end
end
point = z;
end