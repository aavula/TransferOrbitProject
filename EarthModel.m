% FILE: Orbit Model of Earth
% NAME: Akhil Avula, Aaron Abeyta

clear all; close all; clc;

a_E = 1.496e8; %km
e_E = 0.0167; % unitless

nu = 0:(2*pi)/365:(2*pi);

r = (a_E.*(1-e_E.^2))./(1+e_E.*cos(nu));
r_pos = [r.*cos(nu); r.*sin(nu)];

figure;
plotv(r_pos,'-');

% M = moviein(length(nu));
% for j=1:length(nu)
%     plot(x,u(j,:)), M(:,j)=getframe;
% end
% movie(M);







