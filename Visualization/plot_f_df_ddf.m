function plot_f_df_ddf(X,Y, y_scale, x_scale)
if nargin < 3
    x_scale = 'linear';
    y_scale = 'linear';
elseif nargin < 4
    x_scale = 'linear';
end

figure;
subplot(3,1,1)
plot(X,Y)
set(gca, 'XScale', x_scale);
set(gca, 'YScale', y_scale);
grid on
title('Y(x)', 'Interpreter', 'latex');
subplot(3,1,2)
dx = movmean(X, 2, 'Endpoints', 'discard');
dy_dx = diff(Y) ./ dx;
plot(dx, dy_dx);
set(gca, 'XScale', x_scale);
set(gca, 'YScale', y_scale);
grid on
title('$$\frac{dY}{dx}$$', 'Interpreter', 'latex');
subplot(3,1,3)
ddx = movmean(dx,2, 'Endpoints', 'discard');
ddy = diff(dy);
plot(ddx, ddy);
set(gca, 'XScale', x_scale);
set(gca, 'YScale', y_scale);
grid on
title('$$\frac{d^2Y}{dx^2}$$', 'Interpreter', 'latex');

end
