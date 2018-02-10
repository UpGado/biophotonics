function [ wavelength, Os ] = vandehulst( d, wavelength_range )
%VAN-DE-HULST Summary of this function goes here
%   Detailed explanation goes here]
d = d .* 10^-6;
m = 1.42 / 1.36;

% calculate wavelengths from range
wavelength = linspace(wavelength_range(1), wavelength_range(2), 10000);
wavelength = wavelength.* 10^-9;

% calculate sigma
for x = 1:length(d)
    sigma(x,:) = (2 * pi ./ wavelength) * d(x) * abs(m -1);
end

% calculate Os

for x = 1:length(d)
   t1 = 2- 4 * sin(sigma(x,:)) ./ ( sigma(x,:));
   t2 = 4 * (1 - cos(sigma(x,:))) ./ (sigma(x,:).^2);
   Os(x,:) = (pi * d(x)^2 / 4) * (t1 + t2); 
   
   % normalize
   Os(x,:) = Os(x,:)./ Os(x,1);
end

% plot
for x = 1:length(d)
    loglog(wavelength, Os(x,:));
    hold on
end

ylabel('Scattered light intensity');
xlabel('Wavelength (meters)');
end

