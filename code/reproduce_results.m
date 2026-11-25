function reproduce_results(T_max)
    l = 10;
    n = 200;
    dx = 2 * l / n;
    tau = 0.03;
    g = 1;
    bc = 0.2;
    u0 = ones(n+1, 1);
    h0 = ones(n+1, 1);

    x = linspace(-l, l, n+1).';
    v0 = u0 .* x;
    q0 = q(x, bc);
    H = Hfunc(x, bc);
    h0 = 1 - H;

    [u, h, q, v] = shallowWater(dx, l, tau, T_max, u0, h0, q0, v0, H, g);

    subplot(1, 2, 1);
    plot(x, u(:, end));
    xlabel('x');
    ylabel('u');
    title(sprintf('u at t = %.2f', T_max));

    subplot(1, 2, 2);
    plot(x, h(:, end) + H);
    hold on;
    plot(x, H, 'r--');
    xlabel('x');
    ylabel('h + H');
    title(sprintf('h + H at t = %.2f', T_max));
    hold off;
end

function results = q(x, bc)
    results = x;
    left = x <= -2;
    middle = x > -2 & x < 2;
    right = x >= 2;
    results(left) = x(left) + 0.266;
    results(middle) = x(middle) - bc*(x(middle) - x(middle).^3/12);
    results(right) = x(right) - 0.266;
end

function results = Hfunc(x, bc)
    results = zeros(size(x));
    inside = x >= -2 & x <= 2;
    results(inside) = bc*(1 - x(inside).^2/4);
end