function PlotQRbands(BQ, bBQ, B2, bB2, QQ, varargin)

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


qq = [.025 0.05 .16 .5 .84 0.95 .975]; %quantiles for error bands

sBQ = sort(bBQ,2);
sB2 = sort(bB2);

Nboot = length(bB2);

qBQ = sBQ(:,ceil(qq*Nboot));
qB2 = sB2(ceil(qq*Nboot));



% PLOT DISTRIBUTION OF COEFFICIENTS
mat_quant=squeeze(qBQ);
matm=mat_quant;
for i=2:size(matm,2)
    matm(:,i)=matm(:,i)-mat_quant(:,i-1);
    if(isnan(matm(end,i)))
        matm(end,i)=matm(end-1,i);
    end
end
if(isnan(matm(end,1)))
    matm(end,1)=matm(end-1,1);
end
f=figure;
%f = figure
h=area(QQ,matm);
r=0.8;
g=0.8;
b=0.8;
set(h,'LineStyle','none')
set(h(1),'FaceColor',[1 1 1]) % white
set(h(2), 'FaceColor', .95*[r g b])
set(h(3),'FaceColor',.9*[r g b])
set(h(4),'FaceColor',.75*[r g b])
set(h(5),'FaceColor',.75*[r g b])
set(h(6),'FaceColor',.9*[r g b])
set(h(7),'FaceColor',.95*[r g b])
hold on
p1=plot(QQ,BQ,'*-r','LineWidth',2); %fit coefficients
hold on
p2=plot(QQ,mat_quant(:,4),'--k','LineWidth',2);  % median coefficients in black
hold on
p3=plot(QQ,B2*ones(size(QQ)),'--b','LineWidth',2);  % OLS coefficients in blue
hold on
p4=plot(QQ,zeros(size(QQ)),'-k','LineWidth',1);  % zero in black
set(gcf,'Color','w')
axis tight; box on;
hold off;
xlim([QQ(1) QQ(end)])
if exist('ylims', 'var')
    ylim(ylims)
end
xlabel('\tau')
ylabel('\beta(\tau)')
legend([p1 p2 p3], 'In-sample fit','Median', 'OLS', 'Location','Best')

if exist('filename', 'var')
    printpdf(f,filename);
end
