function [X_matched, Y_matched, varargout] = fun_matched_pair_rejection_by_fcestimate_vessel(X_matched, Y_matched, fc_params)



if isempty(X_matched)
    if nargout > 2
        varargout{1} = [];
    end
    return;
end
% Absolute error tolorance
if isfield(fc_params, 'max_disp_error')
    abs_error_th = fc_params.max_disp_error;
else
    abs_error_th = 2;
end
iadj = fc_params.y_sub_idx;
dims = [1024 1536 251];

if ~isfield(fc_params, 'viz')
    viz = false;
else
    viz = fc_params.viz;    
end


dispvec = X_matched - Y_matched;
x = X_matched(:,fc_params.x_sub_idx);
y = dispvec(:, fc_params.y_sub_idx);

y_est = feval(fc_params.model, fc_params.fitting_params, x);
fit_err = abs(y - y_est);

inliers_Q = fit_err < abs_error_th;

% histogram(fit_err)

x_range = min(x) : max(x);
y_range = feval(fc_params.model, fc_params.fitting_params, x_range);

X_matched = X_matched(inliers_Q, :);
Y_matched = Y_matched(inliers_Q, :);

if nargout == 3
    varargout{1} = inliers_Q;
end
if viz
    figure;
    if iadj==1
        ran = [min(dispvec);max(dispvec)];
        xl = ran(:,iadj);
        xl = xl(:)';
        yl = [0 dims(fc_params.x_sub_idx)];
        subplot(2,2,[1 3])
        cla
        xlim(xl)
        ylim(yl)
        plot(y,x,'+')
        hold on
        % Outliers
        plot(y(~inliers_Q),x(~inliers_Q),'ro')
        plot(y_range,x_range,'g-')
        plot(feval(fc_params.model,fc_params.fitting_params,x),x,'go')
        daspect([1 50 1])   
    elseif iadj==2
        subplot(2,1,iadj)
        cla
        plot(x,y,'+')
        hold on
        % Outliers
        plot(x(~inliers_Q),y(~inliers_Q),'ro')
        plot(x_range,y_range,'g-')
        plot(x,y_est,'go')
    end
end
end