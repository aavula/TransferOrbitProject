% FILE: Orbit Model of Earth
% NAME: Akhil Avula, Aaron Abeyta

clear all; close all; clc;

a_E = 1.496e8; %km
e_E = 0.0167; % unitless

nu = 0:(2*pi)/365:(2*pi);

r = (a_E.*(1-e_E.^2))./(1+e_E.*cos(nu));
r_pos = [r.*cos(nu); r.*sin(nu)];

figure;
while true
for k = 1:length(nu)
    f = r_pos(:,k);
    plotv(f,'o-')
    grid
    axis square
    xlim([-2e8 2e8])
    ylim([-2e8 2e8])
    drawnow 
end 
end 







