clear all; close all; clc;

%%%%%%%%%%%
% Run this after obtaining some simulations using main.jl!
%%%%%%%%%%%

%%%%%%%%%% COLORS! %%%%%%%%%%%%%%%%%
clrs(1,:)  = [31 120 180]/256; %[0 0.4470 0.7410]; dark blue
clrs(2,:)  = [115, 156, 112]/256; %[0.8500 0.3250 0.0980]; light green
clrs(3,:)  = [11, 163, 0]/256; %[0.9290 0.6940 0.1250]; dark green
clrs(4,:)  = [166 206 227]/256; %[0.4940 0.1840 0.5560]; light blue
clrs(5,:)  = [173, 101, 101]/256; %[240,128,128]/256; %new light red to not conflict with pink %[251 154 153]/256; %[0.4660 0.6740 0.1880]; light red
clrs(6,:)  = [179 50 50]/256; %[227 26 28]/256; %[0.3010 0.7450 0.9330]; red
clrs(7,:)  = [256 128 0]/256; %[253 141 1.11]/256; %[0.6350 0.0780 0.1840]; organge

clrs(8,:)  = [75,0,130]/256; %indigo
clrs(9,:)  = [243, 204, 255]/256; %[147,112,219]/256; %light purp
clrs(10,:) = [255,20,147]/256; %deep pink
clrs(11,:) = [255,182,193]/256; %light pink
clrs(12,:) = [87, 179, 176]/256; %dark cyan
clrs(13,:) = [72,209,204]/256; %light cyan

clrs(14,:) = [18, 70, 105]/256; %v dark blue
clrs(15,:) = [36, 102, 100]/256; %v dark cyan
clrs(16,:) = [171, 107, 44]/256; %v dark orange
clrs(17,:) = [102, 14, 14]/256; %v dark red
clrs(18,:) = [24, 97, 20]/256; %v dark green

clrs(19,: ) = [250, 178, 105]/256; %light orange

%%%%%%%%%% COLORS! %%%%%%%%%%%%%%%%%

fs = 16;

clrOrder = [1,3,13,8,6,10,7];
clrOrderShifts = [1,4,3,2,13,12,8,6,10,7,19];

% Create figure
fig = figure;
set(fig, 'Position', [1700 50 1000 1000]);

% Define the directory path
directoryPath = './simulations';
files = dir(directoryPath);
filenames = {files.name};

% Filter filenames to include those containing 'tgt' and exclude 'RESTART'
filteredFilenames = filenames(contains(filenames, 'tgt') & ~contains(filenames, 'RESTART'));

for k = 1:length(filteredFilenames)
	dmy2 = dir(fullfile(directoryPath, filteredFilenames{k}));
	dmy2 = dmy2(contains({dmy2.name},filteredFilenames{k}(end-5:end)));

	ax = subplot(5, 3, [4:6]);
	hold(ax, 'on');

	for ii = length(dmy2):-1:1
	    ld = load(fullfile(dmy2(ii).folder, dmy2(ii).name));
	    plot(ax, ld.tme/(60*1000), ld.V, 'Color', [0 0.4470 0.7410 0.2], 'LineStyle', '-');
	    plot(ax, ld.tme/(60*1000), (ld.alfa'*60)-50, 'Color', clrs(6,:), 'LineStyle', '-', 'LineWidth', 2.5);
	    plot(ax, ld.tme/(60*1000), ld.EKvals+30, '-k', 'LineWidth', 0.9);
		if ii == length(dmy2)
			ax_tr = subplot(5,3,3);
			xlineLoc = zoomPlot(ax_tr,ld,fs,'C');
			xline(ax, xlineLoc,'linewidth',1.2,'linestyle','--','LineWidth', 1.5);
			statStruct = measureTrace(ld.V,ld.tme,5910000,6000000);
			% Get Activity Measurements of Steady State
			data = [
			    statStruct.burst.vals(1) statStruct.burst.cv(1);
			    statStruct.burst.vals(2) statStruct.burst.cv(2);
			    statStruct.burst.vals(3) statStruct.burst.cv(3);
			    statStruct.burst.vals(4) statStruct.burst.cv(4);
			    statStruct.burst.vals(5) statStruct.burst.cv(5);
			    statStruct.burst.vals(6) statStruct.burst.cv(6)
			];
			data = data';

			% Create the table
			T = table(data(:, 1), data(:, 2), data(:, 3), data(:, 4), data(:, 5), data(:, 6), ...
			          'VariableNames', {'Period (ms)', 'Burst Duration (ms)', 'Interburst Interval (ms)', 'Spike Height (mV)', 'Max Hyperpolarization (mV)', 'Slow Wave Amplitude (mV)'});

			% Add a new column at the beginning of the table for the row labels
			T.RowLabels = {'Mean'; 'CV'};

			% Rearrange the columns to make 'RowLabels' the first column
			T = [T(:, end) T(:, 1:end-1)];

			% Display the table
			disp(T);
		end
		if ii == 5
			ax_tr = subplot(5,3,1);
			xlineLoc = zoomPlot(ax_tr,ld,fs,'A');
			xline(ax, xlineLoc,'linewidth',1.2,'linestyle','--','LineWidth', 1.5);
		end
		if ii == 8
			ax_tr = subplot(5,3,2);
			xlineLoc = zoomPlot(ax_tr,ld,fs,'B');
			xline(ax, xlineLoc,'linewidth',1.2,'linestyle','--','LineWidth', 1.5);
			text(ax_tr,0.3,45,filteredFilenames{k}(end-5:end), 'Fontsize', 24);
		end
	end

	% Apply settings to ax
	ylim(ax, [-70 15]);
	yticks(ax, [-50 -20 10]);
	yticklabels(ax, {'-50', '', '10'});
	set(ax, 'XTickLabel', []);
	if k == 1
		ylab1 = ylabel(ax, 'Trace (mV)', 'FontSize', 16);
	end
	set(ylab1,'position',[-48 -20 -1]);
	set(ax,'fontsize',fs);
	ax2 = axes('Position', ax.Position,'YAxisLocation', 'right', 'fontsize',fs,'Color', 'none','YColor', clrs(6,:), 'XTick', [], 'Box', 'off');
	yticks(ax2, [.2353 .5882 .9412]);
	yticklabels(ax2, {'0', '', '1'});
	ylabAlfa = ylabel('Alpha (au)','fontsize',fs,'Color',clrs(6,:));
	set(ylabAlfa,'position',[1.016 .55 0]);

	text(ax, 58,22,'A','Fontsize',15);
	text(ax,158,22,'B','Fontsize',15);
	text(ax,956,22,'C','Fontsize',15);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Conductance Plots
	set(gca, 'XTickLabel', {});

	for ii = length(dmy2):-1:1
		ld = load(strcat(dmy2(ii).folder,'/',dmy2(ii).name));

		pltSpcC1 = [13:15];
		pltSpcC2 = [16:18];

		axC1 = subplot(10, 3, pltSpcC1);
		for jj = 1:7
			h1 = plot(axC1,ld.tme(1:end)/(60*1000),ld.conds(jj,1:end),'linewidth',1.7,'color',clrs(clrOrder(jj),:));
			hold on;
		end
		box off;
		axC2 = subplot(10, 3, pltSpcC2);
		for jj = 1:7
			h2 = plot(axC2,ld.tme(1:end)/(60*1000),ld.conds(jj,1:end),'linewidth',1.7,'color',clrs(clrOrder(jj),:));
			hold on;
		end
		box off;
		axC1.XAxis.Visible = 'off'; set(axC2, 'XTickLabel', []);
		if ii == length(dmy2)
			finConds = mean(ld.conds(:,end-2000,end),2);
		end
	end
	ylim(axC1, [1 15]);
	ylim(axC2, [0 1]);
	xlim(axC1,[0 1000]); 
	xlim(axC2,[0 1000]);

	posMat = get(axC2,'Position');
	posMat(4) = posMat(4) + 0.0117;
	posMat = get(axC1,'Position');
	posMat(2) = posMat(2) - .023; posMat(4) = posMat(4) + .0317;
	set(axC1,'Position',posMat);


	% Apply settings to ax
	yticks(axC1, [5 10 15]);
	yticklabels(axC1, {'5', '10', '15'});
	if k == 1
		ylab2 = ylabel(axC1, 'Conductance (nS)', 'FontSize', fs);
	end
	set(ylab2,'Position',[-48 3 -1]);
	yticks(axC2, [0 .5 1]);
	yticklabels(axC2, {'0','','1'});
	set(axC1,'fontsize',fs)
	set(axC2,'fontsize',fs)

	% Add break marks using annotation
	breakPosition = [posMat(1) posMat(2)];  % Adjust this as needed
	annotation('line', [breakPosition(1)-0.002, breakPosition(1)+0.006], [breakPosition(2)-0.003, breakPosition(2)-0.001], 'Color', 'k','linewidth',1.1);
	annotation('line', [breakPosition(1)-0.002, breakPosition(1)+0.006], [breakPosition(2)+0.001, breakPosition(2)+0.003], 'Color', 'k','linewidth',1.1);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Half Activation Plot
	lnType{1} = '-'; lnType{2} = '--';
	ix = [1 1;1 2; 3 1;3 2;13 1;13 2; 8 1; 6 1; 5 1; 7 1; 7 2];

	axH = subplot(5, 3,[10:12]);
	hold(axH, 'on');

	for ii = length(dmy2):-1:1
		ld = load(strcat(dmy2(ii).folder,'/',dmy2(ii).name));
		for jj = 1:11
			plot(axH,ld.tme(1:end)/(60*1000),ld.halfAcs(jj,1:end),'linewidth',1.7,'linestyle',lnType{ix(jj,2)},'color',clrs(ix(jj,1),:));
		end
		if ii == length(dmy2)
			finShifts = mean(ld.halfAcs(:,end-2000,end),2);
		end
	end

	ylim([-20 20])
	if k == 1
		yLabH = ylabel(axH, 'Voltage Shift (mV)', 'FontSize', 16);
	end
	xlabel(axH, 'Time (min)', 'FontSize', 16);
	set(gca,'fontsize',fs);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	relShift = .03;

	axFinCond = subplot(5, 3, [13]);
	hold(axFinCond, 'on');
	barHandle = bar(axFinCond, finConds);
	for ii = 1:7
	    barHandle.FaceColor = 'flat'; % This line is crucial for using custom colors
	    barHandle.CData(ii,:) = clrs(clrOrder(ii),:);
	end
	xlim([.5 7.5])
	ylabel(axFinCond, 'Conductances (nS)')
	xticklabels(axFinCond, {'Na','CaT','CaS','H','Kd','KCa','A'});
	yticks(axFinCond, [0 5 10]);
	yticklabels(axFinCond, {'0', '5', '10'});
	set(axFinCond,'fontsize',13);
	posMatCond = get(axFinCond,'Position');
	posMatCond(3) = posMatCond(3) + relShift;
	posMatCond(2) = posMatCond(2) - .04;
	set(axFinCond,'Position',posMatCond)

	axFinShifts = subplot(5, 3, [14:15]);
	hold(axFinShifts, 'on');
	mVVals = [25.5, 48.9, 27.1, 32.1, 33, 66, 70,12.3, 28.3, 27.2, 56.9]' * -1;
	barz = mVVals-finShifts;
	barHandle2 = bar(axFinShifts,barz);
	yticks([-80 -60 -40 -20 0]);
	yticklabels({'-80','','-40','','0'});
	xlim([.5 11.5])
	xticks([1:11])
	xticklabels({'Na_m','Na_h','CaT_m','CaT_h','CaS_m','CaS_h','H_m','Kd_m','KCa_m','A_m','A_h'});
	for ii = 1:11
	    barHandle2.FaceColor = 'flat'; % This line is crucial for using custom colors
	    barHandle2.CData(ii,:) = clrs(clrOrderShifts(ii),:);
	end
	yLabHndl = ylabel(axFinShifts,'Half-(in)activations (mV)');
	set(axFinShifts,'fontsize',12);
	posMatShifts = get(axFinShifts,'Position');
	posMatShifts(1) = posMatShifts(1) + relShift;
	posMatShifts(3) = posMatShifts(3) - relShift;
	posMatShifts(2) = posMatShifts(2) - .04;
	set(axFinShifts,'Position',posMatShifts)
end

saveas(fig, strcat('./simulations/',filteredFilenames{k}(end-5:end),'.eps'), 'epsc');

% Define local functions at the end of the script
function xDashLoc = zoomPlot(ax_tr,ldStore,fs,ttl)

	ld = ldStore;

	xDashLoc = ld.tme(1)/(60*1000);

	plot(ax_tr, ld.tme(1:2000)/1000-ld.tme(1)/1000,ld.V(1:2000),'linewidth',1.2); hold on;
	plot(ax_tr, ld.tme(1:2000)/1000-ld.tme(1)/1000,ld.EKvals(1:2000)+30,'-k','linewidth',.9,'linestyle','-'); 
	axis off;
	ylim(ax_tr, [-70 15]); xlim([0 1.05]);

	xlims = xlim();
	ylims = ylim();

	% Calculate positions for scale bars
	xLblStart = xlims(2) - 0.4 * (xlims(2) - xlims(1)); % 60% from the right
	xLblHeight = ylims(1);
	yLblStart = ylims(1) + 0.4 * (ylims(2) - ylims(1)); % 60% from the top
	yLblLeftAdjustment = xlims(2) - .02 * (xlims(2) - xlims(1));

	xScaleLength = .25; % Length of the x-axis scale bar
	yScaleLength = 40; % Length of the y-axis scale bar

	% Draw scale bars
	plot(ax_tr, [xLblStart, xLblStart + xScaleLength], [xLblHeight, xLblHeight], 'k', 'LineWidth', 1.2); % X scale bar
	plot(ax_tr, [yLblLeftAdjustment, yLblLeftAdjustment], [yLblStart, yLblStart + yScaleLength], 'k', 'LineWidth', 1.2); % Y scale bar
	% Labels
	text(ax_tr, xLblStart + xScaleLength*.4, xLblHeight + 0.0*yScaleLength, [num2str(xScaleLength) ' sec'], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top','fontweight','bold','FontSize',fs-2);
	text(ax_tr, yLblLeftAdjustment + 1*xScaleLength, yLblStart + .5*yScaleLength, [num2str(yScaleLength) ' mV'], 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle','fontweight','bold','FontSize',fs-2);

	%ps = get(ax_tr,'Position'); ps(2) = ps(2) - .02;
	%set(ax_tr,'Position',ps);
	ttl = title(ax_tr,ttl,'FontSize',fs,'fontweight','bold');

end