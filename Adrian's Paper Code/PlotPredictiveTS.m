function PlotPredictiveTS(YQ, QQ, Time, yh, varargin)

%% Variable arguments
for iArg = 1:2:length(varargin)
    switch varargin{iArg}
        case 'filename'
            filename = varargin{iArg + 1};
        case 'xlims'
            xlims = varargin{iArg + 1};
        case 'ylims'
            ylims = varargin{iArg + 1};
        case 'legendLocation'
            legendLocation = varargin{iArg + 1};
        otherwise
            error(['Unrecognized variable argument name:', varargin{iArg}])
    end
end

if ~exist('xlims', 'var')
    xlims = [Time(1), Time(end)];
end

if ~exist('legendLocation', 'var')
    legendLocation = 'NorthEast';
end

%% select a few important quantiles: 5%, 10%, 25%, 50%, 75%, 90%, 95%
Qsel = [.05 .1 .25 .50 .75 .90 .95];
for j=1:length(Qsel)
    [~, J] = min(abs(QQ-Qsel(j)));
    PosSel(j) = J;
end
mat_quant=YQ(:, PosSel);

%%
matm = mat_quant;
for i = 2:size(matm, 2)
    matm(:, i) = matm(:, i) - mat_quant(:, i - 1);
    if(isnan(matm(end, i)))
        matm(end, i) = matm(end - 1, i);
    end
end
if(isnan(matm(end, 1)))
    matm(end, 1) = matm(end - 1, 1);
end

f = figure;
% Generate plot
h = area(Time, matm);
r = 0.8;
g = 0.8;
b = 0.8;
set(h, 'LineStyle','none')
set(h(1), 'FaceColor', [1 1 1]) % white
set(h(2), 'FaceColor', 0.95 * [r g b])
set(h(3), 'FaceColor', 0.90 * [r g b])
set(h(4), 'FaceColor', 0.75 * [r g b])
set(h(5), 'FaceColor', 0.75 * [r g b])
set(h(6), 'FaceColor', 0.90 * [r g b])
set(h(7), 'FaceColor', 0.95 * [r g b])
hold on
p1=plot(Time, yh, '--b', 'LineWidth', 1); %realized time series
hold on
p2=plot(Time, mat_quant(:, 4), '-k', 'LineWidth', 1);  % median forecast in black
hold off;
set(gcf, 'Color', 'w')

%set(gca,'XTick', Time(1:16:end),'XMinorTick','on')
XTickTimes = Time(rem(year(Time),5) == 0 & month(Time) == 1);
set(gca,'XTick', XTickTimes)

datetick('x', 'YYYY','keeplimits', 'keepticks')
axis tight; box on;
xlim(xlims)
if exist('ylims', 'var')
    ylim(ylims)
end

%%% NOTE: will need to modify this?
%legend([p1 p2], 'Realized','Median', 'Location','Best')
legend([p1 p2], 'Realized','Median', 'Location', legendLocation)

%%% NOTE: will need to modify this
if exist('filename', 'var')
    printpdf(f,filename);
end
