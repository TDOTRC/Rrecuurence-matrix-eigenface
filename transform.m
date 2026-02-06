function [out] = transform(in)
% this is for the transformation between distribution and plots
out=[];
for i=1:size(in)
    out=[out,i*ones(1,in(i))];
end

end