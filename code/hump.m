function hump(T_max)
    l = 10;
    n = 200;
    dx = 2 * l / n;
    tau = 0.03;
    g = 1;
    bc = 0.2;
    lc = 2;
    u0 = ones(n+1, 1);
    h0 = ones(n+1, 1);

    x = linspace(-l, l, n+1).';
    v0 = u0 .* x;
    q0 = q0_hump(x, bc, lc);
    H = H_hump(x, bc, lc);
    h0 = 1 - H;

    [u, h, q, v] = shallowWater(dx, l, tau, T_max, u0, h0, q0, v0, H, g);

    subplot(2, 1, 1);
    plot(x, u(:, end));
    xlabel('x');
    ylabel('u');
    title(sprintf('u at t = %.2f', T_max));

    subplot(2, 1, 2);
    plot(x, h(:, end) + H);
    hold on;
    plot(x, H, 'r--');
    xlabel('x');
    ylabel('h + H');
    title(sprintf('h + H at t = %.2f', T_max));
    hold off;
end

function results = H_hump(x, bc, lc)
    results = zeros(size(x));
    inside = abs(x) <= lc/2;
    results(inside) = bc*(1 - 4*x(inside).^2/lc^2);
end


function results = q0_hump(x, bc, lc)

    cst = bc * lc/3;

    results = x;
    left = x <= -lc/2;
    inside = abs(x) < lc/2;
    right = x >= lc/2;

    results(left) = x(left) + cst;
    results(inside) = (1 - bc) * x(inside) + (4 * bc / (3 * lc^2)) * x(inside).^3;
    results(right) = x(right) - cst;

end