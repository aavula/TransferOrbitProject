% FILE: Orbit Model of Earth
% NAME: Akhil Avula, Aaron Abeyta

clear all; close all; clc;
T = 1; %Number of Earth years to be modelled. 
a_E = 1.496e8; %km
e_E = 0.0167; % unitless
b_E = a_E*sqrt(1-e_E^2);
nu = 0:(2*pi)/365*T:(2*pi)*T;
a_M = 2.2792e8;
e_M = 0.0935;
b_M = a_M*sqrt(1-e_M^2);

r_E = (a_E.*(1-e_E.^2))./(1+e_E.*cos(nu));
r_M= (a_M.*(1-e_M.^2))./(1+e_M.*cos(nu));

r_posE = [r_E.*cos(nu); r_E.*sin(nu)];

r_posM = [r_M.*cos(nu); r_M.*sin(nu)];

xE = -a_E:1e5:a_E;
y_posE = b_E./a_E.*sqrt(a_E.^2-xE.^2);
y_negE = -b_E./a_E.*sqrt(a_E.^2-xE.^2);

xM = -a_M:1e5:a_M;
y_posM = b_M./a_M.*sqrt(a_M.^2-xM.^2);
y_negM = -b_M./a_M.*sqrt(a_M.^2-xM.^2);
figure(1)
while true
for k = 1:length(nu)
    f = r_posE(:,k);
    plotv(f,'o-')
    hold on 
    g = r_posM(:,k);
    plotv(g,'o-')
    grid
    axis square
    xlim([-3e8 3e8])
    ylim([-3e8 3e8])
    plot(xE,y_posE,'b')
    plot(xE,y_negE,'b')
    plot(xM,y_posM,'r')
    plot(xM,y_negM,'r')
    drawnow
    hold off
    
end 

end 







