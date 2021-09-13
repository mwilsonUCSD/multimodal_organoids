function [output_curve]=f_norm1(input_curve)
if (max(input_curve)-min(input_curve))==0
    warning('Standard normalization not possible!')
    output_curve=input_curve/mean(input_curve)-1;
else
    output_curve=(input_curve-min(input_curve))/(max(input_curve)-min(input_curve));
end