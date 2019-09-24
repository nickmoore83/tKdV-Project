clear all; close all;

% Model parameters from experiments
epsi0 = 0.017;  % amplitude-to-depth ratio
del0 = 0.22;    % depth-to-wavelength ratio
Drat = 0.24;    % depth ratio.
% Simulation parameters
Lambda = 16;    % Reference 16
Nw = Lambda/2;  % Reference Lambda/2
MC = 1E3;       % Reference 1E4
% time step parameters
dt = 5E-4;  % Reference 5E-4
nout = 100; % Reference 100
tfin = 10.;  % Reference 10

% Choose upstream or downstream
down = false;
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
C2 = (2/3)*pi^2*del0/Nw^2;
C3 = 3/2*pi^(1/2)*epsi0/del0;
[uarray, duarray] = SolverKdV_SymplecticM4a_MC(...
    C2,C3,Drat,Lambda,MC,theta,gibd,fi,dt,nout,tfin);

% Save output
runtime = round(toc/60, 3, 'significant')
save('output.mat', 'C2','C3','Drat','Lambda',...
    'MC','theta','gibd','fi','dt','nout','tfin',...
    'Nw','runtime','uarray', 'duarray');

% Make a histogram of u.
ulist = reshape(uarray,1,[]);
figure(1); histogram(ulist);
figure(2); histogram(ulist); set(gca,'yscale','log');
skewu = skewness(ulist)
% Make a plot to see the waves.
%x=(-pi:2*pi/J:pi-2*pi/J)';
%contourf(T,x,uarray(:,1,:));