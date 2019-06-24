

tempColors = num2cell([255,105,97;97,168,255;178,223,138;117,112,179;255,179,71;255,231,152]/255,2);
col_vec=[1 4 2 3 5 6];
theColors = tempColors(col_vec);
fh = figure('color','white');

ax = subplot(2,1,1);
h1 = raincloud_plot(data2{1}, 'box_on', 1, 'color', theColors{1}, 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
     'box_col_match', 0);

h2 = raincloud_plot(data2{2}, 'box_on', 1, 'color', theColors{2}, 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0);

h3 = raincloud_plot(data2{3}, 'box_on', 1, 'color', theColors{3}, 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .55, 'box_col_match', 0);

for pl=2:6
	set(h1{pl},'visible','off');
	set(h2{pl},'visible','off');
	set(h3{pl},'visible','off');
end

ax.FontSize = 18;
ax.YLim = [1e-8 1];
ax.YScale = 'log';
ylabel('Density');

axes(ha(2));
h1rd = raincloud_plot(data2{1}, 'box_on', 1, 'color', theColors{1}, 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
     'box_col_match', 0);

h2rd = raincloud_plot(data2{2}, 'box_on', 1, 'color', theColors{2}, 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0);

h3rd = raincloud_plot(data2{3}, 'box_on', 1, 'color', theColors{3}, 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .55, 'box_col_match', 0);

ylabel('Box plots');
set(gca,'YLim', [-0.35 0],'XLim',[-5 10],'fontSize',18,'YTick',[]);

set(h1rd{1},'visible','off');
set(h2rd{1},'visible','off');
set(h3rd{1},'visible','off');

legend([h1{1} h2{1} h3{1}], noiseOptions);
