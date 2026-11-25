function sin(T_max)
    l = 10;
    n = 200;
    dx = 2 * l / n;
    tau = 0.03;
    g = 1;
    bc = 0.2;
    freq = 2;
    u0 = ones(n+1, 1);
    h0 = ones(n+1, 1);

    x = linspace(-l, l, n+1).';
    v0 = u0 .* x;
    q0 = q0_sin(x, bc, freq);
    H = H_sin(x, bc, freq);
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

function results = H_sin(x, bc, freq)
    results = bc * sin(2 * pi * freq * x);
end


function results = q0_sin(x, bc, freq)
    results = x + (bc / (2 * pi * freq)) * (cos(2 * pi * freq * x) - 1);
end