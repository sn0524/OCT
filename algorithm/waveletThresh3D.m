function [ v0 ] = waveletThresh3D( z0, J, K, waveType, threshType, lam, a )
%WAVELETTHRESH Summary of this function goes here
%   Detailed explanation goes here

switch waveType
    case 'haar'
        % daubf(1);
        W = @(x) haar3D(x, J);
        WT = @(w) ihaar3D(w, J);
        
    case 'dwt'
        [h0, h1, g0, g1] = daubf(K);
        af = [h0(:) h1(:)];
        sf = [g0(:) g1(:)];
        W = @(x) dwt3D(x, J, af); 
        WT = @(w) idwt3D(w, J, sf);
        
    case 'dwt_l'
        Dim = size(z0);
        [h0, h1, g0, g1] = daubf(K);
        af = [h0(:) h1(:)];
        sf = [g0(:) g1(:)];
        W = @(x) dwt3D_l(x,J,af); 
        WT = @(w) idwt3D_l(w,J,sf,Dim);

    case 'cdt'
        [Faf, Fsf] = FSfarras;
        [af, sf] = dualfilt1;
        W = @(x) cplxdual3D(x, J, Faf, af);
        WT = @(w) icplxdual3D(w, J, Fsf, sf);

end

switch threshType
    case 'L1'
        thresh = @(x, T, a) soft(x, T);
    case 'log'
        thresh = @(x, T, a) (abs(x)/2 - 0.5/a + sqrt( ( (abs(x)/2 + 0.5/a).^2 )- T/a ) ) ...
            .* sign(x) .* (abs(x)-T >= 0);
    case 'MC'
        thresh = @(x, T, a) (min(abs(x), max( 1/(1-a*T) * (abs(x) - T), 0))) .* sign(x);
end

%%
Wz = W(z0);    
tmp = Wz;

if strcmp(waveType,'cdt')
    % cdt
    for j = 1:J
        for m = 1:2
            for n = 1:2
                for p = 1:2
                    for d = 1:7        
                        tmp{j}{m}{n}{p}{d} = thresh(Wz{j}{m}{n}{p}{d}, lam(j), a(j)) ...
                            - Wz{j}{m}{n}{p}{d};
                    end
                end
            end
        end
    end
    j = j + 1;
    nz = size(Wz{j}{1}{1}{1});
    for m = 1:2
        for n = 1:2
            for p = 1:2
                tmp{j}{m}{n}{p} = zeros(nz);
            end
        end
    end
else
    % dwt (including haar)
    for j = 1:J
        for k = 1:7            
            tmp{j}{k} = thresh(Wz{j}{k}, lam(j), a(j)) - Wz{j}{k};
        end
    end
    j = j + 1;
    tmp{j} = Wz{j} - Wz{j};
end
v0 = z0 + WT(tmp);

end

