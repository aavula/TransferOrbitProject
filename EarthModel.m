% FILE: Orbit Model of Earth
% NAME: Akhil Avula, Aaron Abeyta

clear all; close all; clc;
T = 3; %Number of Earth years to be modelled. 
t = 0:365*T; %Increment of days
a_E = 1.496e8; %km
e_E = 0.0167; % unitless
b_E = a_E*sqrt(1-e_E^2); %semiminor axis
a_M = 2.2792e8;
e_M = 0.0935;
b_M = a_M*sqrt(1-e_M^2);
mu = 6.67408e-11*1.989e30*(1e-9/1.33959e-10);
lambda = 4*pi^2/mu;
omegaE = 2*pi/365.2;
omegaM = 2*pi/687;
nuE = omegaE*t;
nuM = omegaM*t;
nuM2 = omegaM*t+pi./4; %pi/4 is an arbitary phase difference between starting earth and end mars.
r_E = (a_E.*(1-e_E.^2))./(1+e_E.*cos(nuE));
r_M = (a_M.*(1-e_M.^2))./(1+e_M.*cos(nuM));
r_M2= (a_M.*(1-e_M.^2))./(1+e_M.*cos(nuM2));

r_posE = [r_E.*cos(nuE); r_E.*sin(nuE)];

r_posM = [r_M.*cos(nuM); r_M.*sin(nuM)];

r_posM2 = [r_M2.*cos(nuM2); r_M2.*sin(nuM2)];

[xE,yE] = pol2cart(nuE,r_E); %Plots the actual ellipses
[xM,yM] = pol2cart(nuM,r_M);
%deltaT = zeros(1:length(t));

figure(1)
for k = 1:length(t)
    E_M = acos((e_M+cos(nuM(k)))/(1+e_M*cos(nuM(k))));
    E_M2 = acos((e_M+cos(nuM2(k)))/(1+e_M*cos(nuM2(k))));
    f = r_posE(:,k);
    g = r_posM(:,k);
    h = r_posM2(:,k);
    
    fnuE = nuE(k)/(2.*pi) - fix(nuE(k)./(2.*pi)).*2.*pi;
    fnuM2 = nuM2(k)/(2.*pi) - fix(nuM2(k)./(2.*pi)).*2.*pi;
    deltanu = ((nuM(k)/(2.*pi) - fix(nuM(k)./(2.*pi)).*2.*pi))-(nuE(k)/(2.*pi) - fix(nuM(k)./(2.*pi)).*2.*pi);
    deltanu2 = ((nuM2(k)/(2.*pi) - fix(nuM2(k)./(2.*pi)).*2.*pi))-(nuE(k)/(2.*pi) - fix(nuM(k)./(2.*pi)).*2.*pi);
    deltaT = abs((deltanu-deltanu2)/(-omegaM));
    
    
    %if nuM2(k) < nuM(k)
    %deltaT = abs(sqrt(a_M.^3./mu).*(2*pi + E_M2-e_M.*sin(E_M2) - (E_M - e_M.*sin(E_M)))); 
    %else 
    %deltaT = abs(sqrt(a_M.^3./mu).*(E_M2-e_M.*sin(E_M2) - (E_M - e_M.*sin(E_M)))); 
    %end 
    a_s = ((2*pi*deltaT./deltanu2).^2./lambda).^(1/3);
    e_s = (r_E(k).*cos(fnuE) -sqrt(r_E(k).^2.*cos(fnuE).^2+4.*a_s.*(abs(a_s-r_E(k))))./(-2.*a_s));
    
    if fnuE < fnuM2
    nus = fnuE:.0001:fnuM2;
    elseif fnuE > fnuM2
        nus = fnuM2:.0001:fnuE;
    else 
        continue 
    end 
        
    r_s = (a_s.*(1-e_s.^2))./(1+e_s.*cos(nus));
    [xs,ys] = pol2cart(nus,r_s);
    plotv(f,'bo-')
    hold on 
    plotv(g,'ro-')
    plotv(h,'mo-')
    grid
    axis square
    xlim([-3e8 3e8])
    ylim([-3e8 3e8])
    plot(xE,yE, 'b')
    plot(xM,yM, 'r')
    plot(xs,ys)

    drawnow
    hold off
    
end
    
