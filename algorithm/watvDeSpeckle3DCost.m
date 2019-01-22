function [ cost ] = watvDeSpeckle3DCost( y, x, w, s, J, wave, thresh, rho, lam, beta, gamma, dim)
%WATV3DOCTCOST Summary of this function goes here
%   Detailed explanation goes here

beta1 = beta(1);
beta2 = beta(2);
beta3 = beta(3);
a = 1 ./ lam;

penalty0 = rho .* norm(s(:), 1);

switch thresh
    case 'L1'
        fpen = @(x, a) abs(x);
    case 'log'
        fpen = @(x, a) 1/a * log(1 + a*abs(x));
    case 'MC'
        fpen = @(x, a) (abs(x) - a/2 * x.^2).*(abs(x) <= 1/a) + 1/(2*a)*(abs(x) > 1/a);
end

penalty1 = 0;
if strcmp(wave, 'cdt') == 1
    for j = 1:J
        for m = 1:2
            for n = 1:2
                for p = 1:2
                    for d = 1:7
                        penalty1 = penalty1 + ...
                            lam(j).*sum(sum(sum(fpen(w{j}{m}{n}{p}{d}, a(j)))));
                    end
                end
            end
        end
    end
    
else
    for j = 1:J
        for k = 1:7
            penalty1 = penalty1 + lam(j).*sum(sum(sum(fpen(gather(w{j}{k}), a(j)))));
        end
    end

end

tmp = diff(x,1,1);
tvx1 = norm(tmp(:), 1);
tmp = diff(x,1,2);
tvx2 = norm(tmp(:), 1);
tmp = diff(x,1,3);
tvx3 = norm(tmp(:), 1);
penalty2 = beta1*tvx1 + beta2*tvx2 + beta3*tvx3;

sy = sum(y, dim) ./ size(y, dim);
sx = sum(x, dim) ./ size(y, dim);
penalty3 = gamma/2 .* norm(sy(:) - sx(:))^2;

cost = 0.5 .* norm(y(:) - x(:) - s(:))^2 + penalty0 ...
    + penalty1 + penalty2 + penalty3;

end

