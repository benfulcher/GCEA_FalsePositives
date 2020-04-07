% propLitCFPR

% Parameters:
equiprobableBins = true;
plotInOne = true;

% Plot:
f = figure('color','w');
if plotInOne
    hold('on')
    PropLitFPSR('mouse','all',equiprobableBins,false)
    PropLitFPSR('human','cortex',equiprobableBins,false)
    ax = gca();
    ax.XScale = 'log';
    f.Position = [744   871   300   216];
else
    subplot(1,2,1);
    PropLitFPSR('mouse','all',equiprobableBins,false)
    subplot(1,2,2);
    PropLitFPSR('human','cortex',equiprobableBins,false)
    f.Position = [744   871   533   179];
end

% Save:
if plotInOne
    fileName = fullfile('OutputPlots','CFPR_Lit_Together_InOne.svg');
else
    fileName = fullfile('OutputPlots','CFPR_Lit_Together.svg');
end
saveas(f,fileName,'svg')
fprintf(1,'Saved to %s\n',fileName);
