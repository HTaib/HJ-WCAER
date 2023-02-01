function [iters, endpt,fnew] =hooke(nvars, startpt, rho, epsilon, itermax, objective_function)
    [newx, xbefore] = deal(startpt)
    delta = abs(startpt.*rho);
    delta(delta == 0) = rho;
    iadj = 0;
    steplength = rho;
    iters = 0;
    fbefore = objective_function(newx);

    while iters < itermax && steplength > epsilon
        iters = iters + 1;
        iadj = iadj + 1;
        fprintf('\n f(x) =  %12.8f at\n', fbefore);
        newx = xbefore;
        [newf, newx] = beat_nearby(delta, newx, fbefore, nvars, objective_function);
        keep = 1;
        while newf < fbefore && keep == 1
            iadj = 0;
            for i = 1 : nvars
                if newx(i) <= xbefore(i)
                    delta(i) = 0 - abs(delta(i));
                else
                    delta(i) = abs(delta(i));
                end
                tmp = xbefore(i);     
                xbefore(i) = newx(i);
                newx(i) = newx(i) * 2 - tmp;
            end
            fbefore = newf;
            [newf, newx] = beat_nearby(delta, newx, fbefore, nvars, objective_function);
            if newf >= fbefore
                break
            end
            keep = 0;
            for i = 1 : nvars
                keep = 1;
                if abs(newx(i) - xbefore(i)) > 0.5 * abs(delta(i))
                    break
                else
                    keep = 0;
                end
            end            
        end
        if steplength >= epsilon && newf >= fbefore
            steplength = steplength * rho;
            delta = delta * rho;
        end
    end
    endpt = newx;
    fnew=newf;
end
