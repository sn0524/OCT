function [ y_g ] = icdt3DGPU( w_g, Fsf_g, sf_g, J, volSize, rowsKernel, columnsKernel, beamsKernel )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

for j = 1:J
    for m = 1:7
        [w_g{j}{1}{1}{1}{m}, w_g{j}{2}{2}{1}{m}, w_g{j}{2}{1}{2}{m}, w_g{j}{1}{2}{2}{m}] = ...
            pm4inv(w_g{j}{1}{1}{1}{m}, w_g{j}{2}{2}{1}{m}, w_g{j}{2}{1}{2}{m}, w_g{j}{1}{2}{2}{m});
         [w_g{j}{2}{2}{2}{m}, w_g{j}{1}{1}{2}{m}, w_g{j}{1}{2}{1}{m}, w_g{j}{2}{1}{1}{m}] = ...
            pm4inv(w_g{j}{2}{2}{2}{m}, w_g{j}{1}{1}{2}{m}, w_g{j}{1}{2}{1}{m}, w_g{j}{2}{1}{1}{m});
    end
end

y_g = gpuArray.zeros(volSize);

for m = 1:2
    for n = 1:2
        for p = 1:2
            lo_g = w_g{J+1}{m}{n}{p};
            lpf_g{1} = sf_g{m}(:,1);
            lpf_g{2} = sf_g{n}(:,1);
            lpf_g{3} = sf_g{p}(:,1);

            hpf_g{1} = sf_g{m}(:,2);
            hpf_g{2} = sf_g{n}(:,2);
            hpf_g{3} = sf_g{p}(:,2);

            for j = J:-1:2
                lo_g = sfb3DGPU_l2(lo_g, w_g{j}{m}{n}{p}, lpf_g, hpf_g, rowsKernel, columnsKernel, beamsKernel);
                [hiRowSize, hiColSize, hiBeaSize] = size(w_g{j-1}{1}{1}{1}{1});
                lo_g = lo_g(1:hiRowSize,1:hiColSize,1:hiBeaSize);
            end
            
            lpf_g{1} = Fsf_g{m}(:,1);
            lpf_g{2} = Fsf_g{n}(:,1);
            lpf_g{3} = Fsf_g{p}(:,1);

            hpf_g{1} = Fsf_g{m}(:,2);
            hpf_g{2} = Fsf_g{n}(:,2);
            hpf_g{3} = Fsf_g{p}(:,2);
            
            lo_g = sfb3DGPU_l2(lo_g, w_g{1}{m}{n}{p}, lpf_g, hpf_g, rowsKernel, columnsKernel, beamsKernel);
            y_g = y_g + lo_g(1:volSize(1), 1:volSize(2), 1:volSize(3));
        end
    end
end
y_g = y_g/sqrt(8);
end

