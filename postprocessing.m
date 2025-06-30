clear; clc
% This script plots figures 7 and 8 in the accompanying article and exports
% the corresponding data to csv files.
% 
% (c) Jannik Theyssen, LVA INSA Lyon, 2025 (GNU GPLv3)


load('final_study.mat')
l = calc(1).p(1).lambda;

d = zeros(1, length(calc));
d_l = zeros(1, length(calc));
eta = zeros(length(calc),1);
W1_ref = zeros(length(calc),1);
W1_rad = zeros(length(calc),1);
W1_sea = zeros(length(calc),1);

E2_ref = zeros(length(calc),1);
E2_rad = zeros(length(calc),1);
E2_sea = zeros(length(calc),1);

for i = 1:length(calc)
    d(i) = norm(calc(i).spring_plate1 - calc(i).source);
    d_l(i) = d(i)/l;
    
    eta(i,:) = calc(i).p(1).mat.eta;

    g(i) = find_grid_point_idx(calc(i).field_grid, calc(i).spring_plate1(1), calc(i).spring_plate1(2));
    W1_ref(i) = calc(i).ref.energy_field(1, g(i));
    W1_rad(i) = calc(i).rad.energy_field(g(i));
    W1_sea(i) = calc(i).sea.energy(1)/calc(i).p(1).geo.A;

    E2_ref(i,:) = calc(i).ref.energy(2);
    E2_rad(i,:) = calc(i).rad.energy(2);
    E2_sea(i,:) = calc(i).sea.energy(2);
end


% Filter mask based on loss factor
idx01 = eta' == 0.001;  % loss factor 0.1 %
idx10 = eta' == 0.1;  % loss factor 10 % 


% Plot and export data for results with 10% loss factor

figure(7); clf
tiledlayout(3,1,'TileSpacing','tight')
nexttile
h10 = histogram(d_l(idx10),'BinWidth', 0.02/l);
xlim([0, 5.3])
ylabel('Count (-)')

dist_cat10 = h10.BinEdges(1:end-1);
res10 = zeros(13, length(dist_cat10));

for i = 1:length(dist_cat10)
    idx = d_l < h10.BinEdges(i+1) & d_l > h10.BinEdges(i);

    res10(:,i) = [dist_cat10(i);
                prctile(W1_ref(idx & idx10), 15.9);
                prctile(W1_ref(idx & idx10), 84.1);
                mean(W1_ref(idx & idx10));
                mean(W1_rad(idx & idx10)); 
                mean(W1_sea(idx & idx10));    
                prctile(E2_ref(idx & idx10), 15.9);
                prctile(E2_ref(idx & idx10), 84.1);
                mean(E2_ref(idx & idx10));
                mean(E2_rad(idx & idx10)); 
                mean(E2_sea(idx & idx10));
                h10.BinCounts(i);
                h10.BinEdges(i)];
end


nexttile
plot(dist_cat10, 10*log10(res10(2,:)),'-.r')
hold on
plot(dist_cat10, 10*log10(res10(3,:)),'-.r')
plot(dist_cat10, 10*log10(res10(4,:)),'r','LineWidth',2)
plot(dist_cat10, 10*log10(res10(5,:)),'k','LineWidth',2)
plot(dist_cat10, 10*log10(res10(6,:)),'-o','LineWidth',2)
xlim([0, 5.3])
xlabel('Wavelengths (-)')
ylabel('W1(c1) (dB re 1 J/m2)')

nexttile
plot(dist_cat10, 10*log10(res10(7,:)),'-.r')
hold on
plot(dist_cat10, 10*log10(res10(8,:)),'-.r')
plot(dist_cat10, 10*log10(res10(9,:)),'r','LineWidth',2)
plot(dist_cat10, 10*log10(res10(10,:)),'k','LineWidth',2)
plot(dist_cat10, 10*log10(res10(11,:)),'-o','LineWidth',2)
xlim([0, 5.3])
xlabel('Wavelengths (-)')
ylabel('E2 (dB re 1 J)')


labels = ["dist", "W1_ref_prctile16", "W1_ref_prctile84", "W1_ref_mean", "W1_rad_mean", "W1_sea_mean", "E2_ref_prctile16", "E2_ref_prctile84", "E2_ref_mean", "E2_rad_mean", "E2_sea_mean", "binCounts", "binStartEdge"];
T10 = array2table(res10', 'VariableNames', labels);
writetable(T10,'result_eta_10%.txt')


% Plot and export data for results with 0.1% loss factor

figure(8); clf
tiledlayout(3,1,'TileSpacing','tight')
nexttile
h01 = histogram(d_l(idx01),'BinWidth', 0.02/l);
xlim([0, 5.3])
ylabel('Count (-)')

dist_cat01 = h01.BinEdges(1:end-1);
res01 = zeros(13, length(dist_cat01));

for i = 1:length(dist_cat01)
    idx = d_l < h01.BinEdges(i+1) & d_l > h01.BinEdges(i);

    res01(:,i) = [dist_cat01(i);
                 prctile(W1_ref(idx & idx01), 15.9);
                 prctile(W1_ref(idx & idx01), 84.1);
                 mean(W1_ref(idx & idx01));
                 mean(W1_rad(idx & idx01)); 
                 mean(W1_sea(idx & idx01));
                 prctile(E2_ref(idx & idx01), 15.9);
                 prctile(E2_ref(idx & idx01), 84.1);
                 mean(E2_ref(idx & idx01));
                 mean(E2_rad(idx & idx01)); 
                 mean(E2_sea(idx & idx01));
                 h01.BinCounts(i);
                 h01.BinEdges(i)];
end

nexttile
plot(dist_cat01, 10*log10(res01(2,:)),'-.r')
hold on
plot(dist_cat01, 10*log10(res01(3,:)),'-.r')
plot(dist_cat01, 10*log10(res01(4,:)),'r','LineWidth',2)
plot(dist_cat01, 10*log10(res01(5,:)),'k','LineWidth',2)
plot(dist_cat01, 10*log10(res01(6,:)),'-o','LineWidth',2)
xlabel('Wavelengths (-)')
ylabel('W1(c1) (dB re 1 J/m2)')
xlim([0, 5.3])


nexttile
plot(dist_cat01, 10*log10(res01(7,:)),'-.r')
hold on
plot(dist_cat01, 10*log10(res01(8,:)),'-.r')
plot(dist_cat01, 10*log10(res01(9,:)),'r','LineWidth',2)
plot(dist_cat01, 10*log10(res01(10,:)),'k','LineWidth',2)
plot(dist_cat01, 10*log10(res01(11,:)),'-o','LineWidth',2)
xlim([0, 5.3])
xlabel('Wavelengths (-)')
ylabel('E2 (dB re 1 J)')


T01 = array2table(res10', 'VariableNames', labels);
writetable(T01,'result_eta_01%.txt')
