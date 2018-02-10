function [ wavelength, Ip ] = rayleigh( L, theta, wavelength_range)
%RAYLEIGH Calculate Rayleigh-Gans approximation of scattered lights
%intensity
%   Detailed explanation goes here]
L = L.* 10^-9;
m = 1.05;

% calculate wavelengths from range
wavelength = linspace(wavelength_range(1), wavelength_range(2), 10000);
wavelength = wavelength.* 10^-9;

k = 2 * pi ./ wavelength;

for x = 1:length(L)
    u(x,:) = 2 .* k * L(x) * sin(theta/2);
end

% form factor
for x = 1:length(L)
    P(x,:) = ((3 ./ u(x,:).^3).* (sin(u(x,:)) - u(x,:) .* cos(u(x,:)))).^2;
end

% volume of scatterer
r = L./ 2;
V = (4*pi/3) .* (r).^3;

% scattered light intensity
for x = 1:length(L)
    Ip(x,:) = ( (pi^2 * V(x)^2) ./ (r(x)^2 * wavelength.^4)) .* (m -1).^2 * (cos(theta)^2) .* P(x,:);
    
    %normalize
    Ip(x,:) = Ip(x,:)./ Ip(x,1);
end

% plot
for x = 1:length(L)
    loglog(wavelength, Ip(x,:));
    hold on
end

ylabel('Scattered light intensity');
xlabel('Wavelength (meters)');
end
