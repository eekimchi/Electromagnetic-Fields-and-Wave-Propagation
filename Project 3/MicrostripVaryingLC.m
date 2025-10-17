clc
clear
close all


fgraph = linspace(2,6,701);

firstcase = readmatrix('0pF 0nH.csv');
firstcase = firstcase';
secondcase = readmatrix('0pF 20nH.csv');
secondcase = secondcase';
thirdcase = readmatrix('20pF 0nH.csv');
thirdcase = thirdcase';
fourthcase = readmatrix('20pF 20nH.csv');
fourthcase = fourthcase';

figure(1)
title('S11 Parameters','FontSize', 20)
xlabel('Frequency, GHz','FontSize', 14)
ylabel('Phase, ang','FontSize', 14)
xlim([5.69 5.76])
ylim([150 180])
hold on
plot(fgraph, firstcase(4,:), 'LineWidth', 3, 'Color', 'cyan')
hold on
plot(fgraph, secondcase(4,:), 'LineWidth', 3, 'Color', 'magenta')
hold on
plot(fgraph, thirdcase(4,:), 'LineWidth', 3, 'Color', 'green')
hold on
plot(fgraph, secondcase(4,:), 'LineWidth', 3, 'Color', 'blue')
legend('0pF 0nH', '0pF 20 nH', '20pF 0nH', '20pF 20nH','Location', 'best')

figure(2)
title('S11 Parameters','FontSize', 20)
xlabel('Frequency, GHz','FontSize', 14)
ylabel('Phase, ang','FontSize', 14)
xlim([5.73294 5.73296])
ylim([165.2 165.21])
hold on
plot(fgraph, firstcase(4,:), 'LineWidth', 3, 'Color', 'cyan')
hold on
plot(fgraph, secondcase(4,:), 'LineWidth', 3, 'Color', 'magenta')
hold on
plot(fgraph, thirdcase(4,:), 'LineWidth', 3, 'Color', 'green')
hold on
plot(fgraph, secondcase(4,:), 'LineWidth', 3, 'Color', 'blue')
legend('0pF 0nH', '0pF 20 nH', '20pF 0nH', '20pF 20nH','Location', 'best')



