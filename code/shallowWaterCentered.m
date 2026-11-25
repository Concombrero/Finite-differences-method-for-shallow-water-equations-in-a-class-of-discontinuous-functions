function [u, h, q, v] = shallowWaterCentered(dx, l, tau, T_max, u0, h0, q0, v0, H, g)
    x = linspace(-l, l, 2*l/dx + 1);
    t = linspace(0, T_max, T_max/tau + 1);

    u = zeros(length(x), length(t));
    h = zeros(length(x), length(t));
    q = zeros(length(x), length(t));
    v = zeros(length(x), length(t));

    u(:,1) = u0;
    h(:,1) = h0;
    q(:,1) = q0;
    v(:,1) = v0;

    for k = 1:length(t)-1
        v(:,k+1) = v(:,k) - tau/2 * u(:,k).^2 - tau*g*h(:,k) - tau*g*H(:);
        q(:,k+1) = q(:,k) - tau * u(:,k) .* h(:,k);

        u(2:end-1,k+1) = (v(3:end,k+1) - v(1:end-2,k+1)) / (2*dx);
        u([1,end],k+1) = (v(2,k+1) - v(end-1,k+1)) / (2*dx);

        h(2:end-1,k+1) = (q(3:end,k+1) - q(1:end-2,k+1)) / (2*dx);
        h([1,end],k+1) = (q(2,k+1) - q(end-1,k+1)) / (2*dx);
    end
end
