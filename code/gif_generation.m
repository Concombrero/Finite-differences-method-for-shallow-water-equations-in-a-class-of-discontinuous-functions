function gif_generation()
    l = 10;
    n = 200;
    dx = 2 * l / n;
    tau = 0.03;
    T_max = 0.9;
    frameDelay = 0.15;
    animationFile = 'shallow_water_animation.gif';
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

    t = linspace(0, T_max, size(u, 2));
    fig = figure;

    axU = subplot(2, 1, 1, 'Parent', fig);
    uPlot = plot(axU, x, u(:, 1), 'b', 'LineWidth', 1.5);
    xlabel(axU, 'x');
    ylabel(axU, 'u');
    title(axU, sprintf('u, t = %.2f', t(1)));
    xlim(axU, [x(1), x(end)]);
    grid(axU, 'on');

    axH = subplot(2, 1, 2, 'Parent', fig);
    surfacePlot = plot(axH, x, h(:, 1) + H, 'b', 'LineWidth', 1.5);
    xlabel(axH, 'x');
    ylabel(axH, 'h + H');
    title(axH, sprintf('h + H, t = %.2f', t(1)));
    xlim(axH, [x(1), x(end)]);
    grid(axH, 'on');

    for k = 1:length(t)
        set(uPlot, 'YData', u(:, k));
        set(surfacePlot, 'YData', h(:, k) + H);
        title(axU, sprintf('u, t = %.2f', t(k)));
        title(axH, sprintf('h + H, t = %.2f', t(k)));
        drawnow;
        frame = getframe(fig);
        gifOptions = {'gif', 'DelayTime', frameDelay};
        if k == 1
            gifOptions = [gifOptions, {'LoopCount', Inf}];
        else
            gifOptions = [gifOptions, {'WriteMode', 'append'}];
        end
        if exist('OCTAVE_VERSION', 'builtin')
            imwrite(frame.cdata, animationFile, gifOptions{:});
        else
            [frameData, colorMap] = rgb2ind(frame.cdata, 256);
            imwrite(frameData, colorMap, animationFile, gifOptions{:});
        end
        pause(frameDelay);
    end
    fprintf('Animation enregistrée dans %s\n', animationFile);
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