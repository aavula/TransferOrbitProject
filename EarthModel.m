% FILE: Orbit Model of Earth
% NAME: Akhil Avula, Aaron Abeyta

clear all; close all; clc;
T = 10; %Number of Earth years to be modelled. 
t = 0:365*T; %Increment of days
a_E = 1.496e8; %km
e_E = 0.0167; % unitless
b_E = a_E*sqrt(1-e_E^2); %semiminor axis
nu = 0:(2*pi)/365*T:(2*pi);
a_M = 2.2792e8;
e_M = 0.0935;
b_M = a_M*sqrt(1-e_M^2);
omegaE = 2*pi/365.2;
omegaM = 2*pi/687;
nuE = omegaE*t;
nuM = omegaM*t;
r_E = (a_E.*(1-e_E.^2))./(1+e_E.*cos(nuE));
r_M= (a_M.*(1-e_M.^2))./(1+e_M.*cos(nuM));

r_posE = [r_E.*cos(nuE); r_E.*sin(nuE)];

r_posM = [r_M.*cos(nuM); r_M.*sin(nuM)];

[xE,yE] = pol2cart(nuE,r_E); %Plots the actual ellipses
[xM,yM] = pol2cart(nuM,r_M);

figure(1)
for k = 1:length(t)
    f = r_posE(:,k);
    plotv(f,'o-')
    hold on 
    g = r_posM(:,k);
    plotv(g,'o-')
    grid
    axis square
    xlim([-3e8 3e8])
    ylim([-3e8 3e8])
    plot(xE,yE, 'b')
    plot(xM,yM, 'r')
    drawnow
    hold off
    
end 





