function plot_travis_profiles(filename)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    temp = importdata(filename,' ',2);
    temp = temp.data;
    profiles.reff = temp(:,1);
    profiles.ne = temp(:,2);
    profiles.te = temp(:,3);
    profiles.Zeff = temp(:,4);
f=figure;
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
nexttile
ax{1}=gca;
hold on
nexttile
ax{2}=gca;
hold on
    plot(ax{1},profiles.reff,profiles.te)
    plot(ax{2},profiles.reff,profiles.ne/1e19)
xlabel(ax{1},'r_{eff}/a [-]')
xlabel(ax{2},'r_{eff}/a [-]')
ylabel(ax{2},'n_e 10^{19} [m^{-3}]')
ylabel(ax{1},'T_e [keV]')  
end