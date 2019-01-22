function [ w_g ] = cdt3DGPU( x_g, Faf_g, af_g, J, rowsKernel, columnsKernel, beamsKernel )
%CDT3DGPU Summary of this function goes here
%   Detailed explanation goes here

% normalization
x_g = x_g/sqrt(8);
w_g = cell(J+1,1);
for m = 1:2
    for n = 1:2
        for p = 1:2
            Flpf_g{1} = Faf_g{m}(:,1);
            Flpf_g{2} = Faf_g{n}(:,1);
            Flpf_g{3} = Faf_g{p}(:,1);

            Fhpf_g{1} = Faf_g{m}(:,2);
            Fhpf_g{2} = Faf_g{n}(:,2);
            Fhpf_g{3} = Faf_g{p}(:,2);
            
            [lo_g, w_g{1}{m}{n}{p}] = afb3DGPU_l2(x_g, Flpf_g, Fhpf_g, rowsKernel, columnsKernel, beamsKernel);
            for j = 2:J
                lpf_g{1} = af_g{m}(:,1);
                lpf_g{2} = af_g{n}(:,1);
                lpf_g{3} = af_g{p}(:,1);

                hpf_g{1} = af_g{m}(:,2);
                hpf_g{2} = af_g{n}(:,2);
                hpf_g{3} = af_g{p}(:,2);
                [lo_g, w_g{j}{m}{n}{p}] = afb3DGPU_l2(lo_g, lpf_g, hpf_g, rowsKernel, columnsKernel, beamsKernel);
            end
            w_g{J+1}{m}{n}{p} = lo_g;
        end
    end
end

for j = 1:J
    for m = 1:7
        [w_g{j}{1}{1}{1}{m}, w_g{j}{2}{2}{1}{m}, w_g{j}{2}{1}{2}{m}, w_g{j}{1}{2}{2}{m}] = ...
            pm4(w_g{j}{1}{1}{1}{m}, w_g{j}{2}{2}{1}{m}, w_g{j}{2}{1}{2}{m}, w_g{j}{1}{2}{2}{m});
         [w_g{j}{2}{2}{2}{m}, w_g{j}{1}{1}{2}{m}, w_g{j}{1}{2}{1}{m}, w_g{j}{2}{1}{1}{m}] = ...
            pm4(w_g{j}{2}{2}{2}{m}, w_g{j}{1}{1}{2}{m}, w_g{j}{1}{2}{1}{m}, w_g{j}{2}{1}{1}{m});
    end
end

end

