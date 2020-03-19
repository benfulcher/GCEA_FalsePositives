
% Parameters:
equiprobableBins = true;

% Plot:
f = figure('color','w');
subplot(1,2,1);
PropLitFPSR('mouse','all',equiprobableBins,false)
subplot(1,2,2);
PropLitFPSR('human','cortex',equiprobableBins,false)
f.Position = [744   871   533   179];

% Save:
fileName = fullfile('OutputPlots','CFPR_Lit_Together.svg');
saveas(f,fileName,'svg')
fprintf(1,'Saved to %s\n',fileName);
