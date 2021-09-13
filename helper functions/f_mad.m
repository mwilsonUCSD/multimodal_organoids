function [absDevNorm,cutoff,ind]=f_mad(in, threshold)

tmpDummy=1:length(in);
absDev=abs(in-median(in));
absDevMedian=1.4826*median(absDev);
absDevNorm=absDev/absDevMedian;
cutoff=threshold/absDevMedian;
ind=f_clean_pulses(tmpDummy(absDevNorm>cutoff));
end