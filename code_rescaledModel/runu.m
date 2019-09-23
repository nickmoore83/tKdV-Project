clear all; close all;

% Model parameters from experiments
epsi0 = 0.017;  % amplitude-to-depth ratio
del0 = 0.22;    % depth-to-wavelength ratio
Drat = 0.24;    % depth ratio.
% Simulation parameters
Lambda = 16;
Nw = 8;
MC = 1e1;
% time step parameters
dt = 5E-4;
nout = 10; % Default 100?
tfin = 1;

% Choose upstream or downstream
down = true;
if down == true
    theta = 13;
    gibd = 0.24;
    fi = 1;
else
    theta = 15;
    gibd = 1.;
    fi = 0;
end

% Run the simulation.
tic;
C2 = (2/3)*pi^2*del0/Nw^2
C3 = 3/2*pi^(1/2)*epsi0/del0
SolverKdV_SymplecticM4a_MC(C2,C3,Drat,Lambda,MC,theta,gibd,fi,dt,nout,tfin);
cputime = round(toc/60, 2, 'significant')

% Load the output file.
load('output.mat')
JJ = size(uarray,1);
MC = size(uarray,2);
Ntsteps = size(uarray,3);
% Flatten the array of u into a single list.
ulist = reshape(uarray,1,[]);
dulist = reshape(duarray,1,[]);
% Write output to a text file.
fileID = fopen('ulist.txt','w');
fprintf(fileID,'# Input parameters: ');

%fprintf(fileID,'value of Nw, theta, tfin.\n');

fprintf(fileID,'%d\n',[]);
%fprintf(fileID,'%d\n',[JJ,MC,Ntsteps,Nw,theta,tfin]);

fprintf(fileID,'# Values of u\n'); fprintf(fileID,'%9.5f\n', ulist);
fprintf(fileID,'# Values of du\n'); fprintf(fileID,'%9.5f\n', dulist);
fclose(fileID);

% Make a plot to see the waves.
%x=(-pi:2*pi/J:pi-2*pi/J)';
%contourf(T,x,uarray(:,1,:));

% Make a histogram of u.
figure(1); histogram(ulist);
figure(2); histogram(ulist); set(gca,'yscale','log');
skewu = skewness(ulist)