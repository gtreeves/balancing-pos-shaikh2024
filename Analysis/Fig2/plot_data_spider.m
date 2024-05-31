% Initialize data points
% D1 = [5 3 9 1 2];
% D2 = [5 8 7 2 9];
% D3 = [8 2 1 4 6];
% P = [D1; D2; D3];
figure
PO = [trise;NAR;phi]';
PO([2,4,6,7,8,9],:) = [];
% Delete variable in workspace if exists
if exist('s', 'var')
    delete(s);
end

% Spider plot
s = spider_plot_class(PO);

% Legend properties
s.LegendLabels = {'D1', 'D2', 'D3'};
s.LegendHandle.Location = 'northeastoutside';