%% Spectrin estimates vs. data

cell2_mean = [0.6999 0.3001];
cell2_stdev = [0.0236 0.0236];
cell2_est = [0.700 0.300];

cell4_mean = [0.4921 0.2830 0.1284 0.0965];
cell4_stdev = [0.0577 0.0349 0.0275 0.0350];
cell4_est = [0.475 0.275 0.125 0.125];

cell8_mean = [0.3272 0.2549 0.1417 0.0950 0.0567 0.0402 0.0396 0.0446];
cell8_stdev = [0.0471 0.0320 0.0252 0.0155 0.0303 0.0232 0.0259 0.0180];
cell8_est = [0.3583 0.2250 0.125 0.125 0.0417 0.0417 0.0417 0.0417];

cell16_mean = [0.2052 0.1659 0.1057 0.1088 0.0629 0.0534 0.0516 0.0486 0.0307 0.0274 0.0245 0.0225 0.0217 0.0264 0.0230 0.0216];
cell16_stdev = [0.0233 0.0181 0.0184 0.0185 0.0120 0.0099 0.0142 0.0122 0.0132 0.0065 0.0088 0.0093 0.0084 0.0070 0.0101 0.0088];
cell16_est = [0.24 0.16 0.100 0.100 0.050 0.050 0.050 0.050 0.025 0.025 0.025 0.025 0.025 0.025 0.025 0.025];


figure; box on;
errorbar(1:2,cell2_mean, cell2_stdev, 'ko', 'MarkerSize',14, 'LineWidth',1.5, 'MarkerFaceColor','k'); hold on;
plot(1:2,cell2_est,'o','Color',[0.9 0.4 0.2],'MarkerSize',14, 'LineWidth',1.5,'MarkerFaceColor',[0.9 0.4 0.2]);
axis([0.8 2.2 0.25 0.75])
xlabel('$\#$ Cells','interpreter','latex')
xticks([1 2])
ylabel('Fusome Vol. Ratios','interpreter','latex')
h = gca;
h.FontSize = 36;

figure; box on;
errorbar(1:4,cell4_mean, cell4_stdev, 'ko', 'MarkerSize',14, 'LineWidth',1.5, 'MarkerFaceColor','k'); hold on;
plot(1:4, cell4_est,'o','Color',[0.9 0.4 0.2],'MarkerSize',14, 'LineWidth',1.5,'MarkerFaceColor',[0.9 0.4 0.2]);
axis([0.8 4.2 0.05 0.6])
xlabel('$\#$ Cells','interpreter','latex')
xticks(1:4)
ylabel('Fusome Vol. Ratios','interpreter','latex')
h = gca;
h.FontSize = 36;

figure; box on;
errorbar(1:8, cell8_mean, cell8_stdev, 'ko', 'MarkerSize',14, 'LineWidth',1.5, 'MarkerFaceColor','k'); hold on;
plot(1:8, cell8_est,'o','Color',[0.9 0.4 0.2],'MarkerSize',14, 'LineWidth',1.5,'MarkerFaceColor',[0.9 0.4 0.2]);
axis([0.8 8.2 0 0.4])
xlabel('$\#$ Cells','interpreter','latex')
xticks(1:8)
ylabel('Fusome Vol. Ratios','interpreter','latex')
h = gca;
h.FontSize = 36;

figure; box on;
errorbar(1:16, cell16_mean, cell16_stdev, 'ko', 'MarkerSize',14, 'LineWidth',1.5, 'MarkerFaceColor','k'); hold on;
plot(1:16, cell16_est,'o','Color',[0.9 0.4 0.2],'MarkerSize',14, 'LineWidth',1.5,'MarkerFaceColor',[0.9 0.4 0.2]);
axis([0.8 16.2 0 0.25])
xlabel('$\#$ Cells','interpreter','latex')
xticks(1:16)
ylabel('Fusome Vol. Ratios','interpreter','latex')
h = gca;
h.FontSize = 36;
legend('Data','Model','location','north')
