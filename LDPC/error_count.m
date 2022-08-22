function [md,fa]=error_count(B,B_hat)
    fa = sum(~ismember(B_hat,B,'rows'))/size(B,1);  %false alarm
    md = sum(~ismember(B,B_hat,'rows'))/size(B,1); %misdetection   
end