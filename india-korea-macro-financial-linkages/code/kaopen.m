T = readtable('kaopen_2023.xls');

disp('File loaded. Columns:')
disp(T.Properties.VariableNames)

disp('Country names with India:')
disp(unique(T.country_name(contains(T.country_name,'India','IgnoreCase',true))))
disp('Country names with Korea:')
disp(unique(T.country_name(contains(T.country_name,'Korea','IgnoreCase',true))))

india = T(contains(T.country_name,'India','IgnoreCase',true), :);
korea = T(contains(T.country_name,'Korea','IgnoreCase',true), :);
if numel(unique(korea.country_name)) > 1
    korea = korea(contains(korea.country_name,'Rep','IgnoreCase',true) | ...
        contains(korea.country_name,'South','IgnoreCase',true), :);
end
india = india(india.year>=2000 & india.year<=2023, :);
korea = korea(korea.year>=2000 & korea.year<=2023, :);

fprintf('India rows: %d, Korea rows: %d\n', height(india), height(korea));
fprintf('India ka_open: mean=%.2f, 2000=%.2f, 2023=%.2f\n', mean(india.ka_open,'omitnan'), india.ka_open(1), india.ka_open(end));
fprintf('Korea ka_open: mean=%.2f, 2000=%.2f, 2023=%.2f\n', mean(korea.ka_open,'omitnan'), korea.ka_open(1), korea.ka_open(end));

figure('Position',[100 100 800 450]);
plot(india.year, india.ka_open, '-o','LineWidth',1.5); hold on;
plot(korea.year, korea.ka_open, '-s','LineWidth',1.5);
ylim([0 1.05]); grid on;
xlabel('Year'); ylabel('Chinn-Ito index (0 = closed, 1 = fully open)');
legend('India','South Korea','Location','east');
title('Chinn-Ito capital-account openness, 2000-2023');
saveas(gcf,'Figure3_openness.png');
fprintf('Saved Figure3_openness.png\n');