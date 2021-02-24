function [thrs,medcurvemodel] = estimateFCthreshols(XX,perc)
medcurvemodel = median(XX,3);
dims = size(XX);
thrs = nan(dims(1:2));
for ii=1:dims(1)
    for jj=1:dims(2)
        xx = squeeze(XX(ii,jj,:));
        [vals,bins] = hist(double(xx),1e3);
        [x12,locmax] = max(vals);
        x1 = [bins(locmax) x12]';
        if perc
            % find the percentile that has %95 of right hand data
            idx2 = find((cumsum(vals)/sum(vals(:)))>=perc,1,'first');
        else
            idx2 = length(bins);
        end
        x2 = [bins(idx2) 0]';
        x0 = [bins(locmax+1:end);vals(locmax+1:end)] ;
        % make sure solution is in the convex region (this is necessary for perc
        % calculations)
        x0 = x0(:,x0(1,:)<x2(1) & x0(1,:)>x1(1));
        [maxval,maxarg,d]=dist2line(x1,x2,x0);
        if numel(maxarg)
            threshold = maxarg(1);
        else
            maxIn = max(xx(:));
            threshold = max(1,graythresh(xx/maxIn)*maxIn); % for heavy peaked distributions, OTSU returns 0
        end
        thrs(ii,jj)=abs(threshold-medcurvemodel(ii,jj));
    end
end
