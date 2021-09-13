function [output_curve]=f_clean_pulses(input_curve)
output_curve=input_curve;
for i=length(output_curve):-1:2
    if output_curve(i)==output_curve(i-1)+1
        output_curve(i)=[];
    end
end
end