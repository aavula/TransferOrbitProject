% An Introduction to the Mathematical Theory of Waves
% Roger Knobel

clear all; close all; clc;

x = -10:0.1:10; 
t = 0:0.3:6;
[X,T] = meshgrid(x,t);
u = exp(-(X-T).^2);
M = moviein(length(t));
for j=1:length(t)
    plot(x,u(j,:)), M(:,j)=getframe;
end
movie(M);

