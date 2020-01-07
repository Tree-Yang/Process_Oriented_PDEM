function PDEM_post(tm, rm, prob, type, time, rlim, Respon_Name)
%PDEM_post - Description
%solve the GDEE by finite difference method with TVD scheme

%by J. S. Yang
%date: 2019-11-05
%====================================================================
    %INUPT:
    %----------------------------
    %tm         : mesh of time
    %rm         : mesh of displacement
    %prob       : matrix of time variant probability density
    %type       : type of post-process
    %                'surf'    : surface of probability density function
    %                'contour' : contour of probability density function
    %                'curve'   : probability density function
    %time       : time instance for 'curve', and time interval for 'surf' and 'contour'
    %rlim       : limit of response (xlim for figures)
    %Respon_Name: xlabel for figures
%=====================================================================

    if length(rlim) ~= 2
        error('The <rlim> must be a 2-component positive array!');
    end
    rlim(2) = max(rlim);
    rlim(1) = min(rlim);
    if strcmp(type, 'surf')   
        if length(time) ~= 2 || min(time) < 0
            error('The <time> must be a 2-component positive array, indicating the starting and ending time!');
        end
        time(2) = max(time);
        time(1) = min(time);
        t0      = tm(1,:);
        [~,n_time1] = min(abs(t0-time(1)));
        [~,n_time2] = min(abs(t0-time(2)));
        figure;
        surf(tm(:, n_time1:n_time2), rm(:, n_time1:n_time2), prob(:, n_time1:n_time2));
        shading interp; 
        colormap jet;
        colorbar;
        xlim([t0(n_time1),t0(n_time2)]);
        ylim([rlim(1), rlim(2)]);
        zlim([0, 1.25*max(max(prob(:, n_time1:n_time2)))]);
        xlabel('Time[s]');
        ylabel(Respon_Name);
        zlabel('PDF');
%         box on; grid off;
        set(gca, 'FontName', 'Arial', 'FontSize', 14);  
    elseif strcmp(type, 'contour')
        if length(time) ~= 2 || min(time) < 0
            error('The <time> must be a 2-component positive array, indicating the starting and ending time!');
        end
        time(2) = max(time);
        time(1) = min(time);
        t0      = tm(1,:);
        [~,n_time1] = min(abs(t0-time(1)));
        [~,n_time2] = min(abs(t0-time(2)));
        figure;
        isoline = linspace(0,max(max(prob(:, n_time1:n_time2))),15);
        contour(tm(:, n_time1:n_time2), rm(:, n_time1:n_time2), prob(:, n_time1:n_time2), isoline);
        colormap jet;
        colorbar;
        xlim([t0(n_time1),t0(n_time2)]);
        ylim([rlim(1), rlim(2)]);
        xlabel('Time[s]');
        ylabel(Respon_Name);
        box on; grid off;
        set(gca, 'FontName', 'Arial', 'FontSize', 14); 
    elseif strcmp(type, 'curve')
        color = [80,81,79; 242,95,92; 255,223,38; 36,123,159; 112,193,179]/255;
        linsty = {'-','--','-.',':'};
        figure; hold on;
        for ii = 1:1:length(time)
            [~,n_time] = min(abs(tm(1,:)-time(ii)));
            plot(rm(:,n_time), prob(:,n_time), 'LineStyle', linsty{mod(ii, 4)+1}, 'Color',color(mod(ii,5)+1,:),'LineWidth', 1.5);
            lgd{ii} = ['\itt = \rm', num2str(time(ii),'%.1f'),'s'];
        end
        xlim([rlim(1), rlim(2)]);
        xlabel(Respon_Name);
        ylabel('PDF');
        legend(lgd,'Location','best');
        %----------------------------------------------------
%         %method II
%         lgd_str = [];
%         for ii = 1:1:length(time)
%             lgd_str = [lgd_str, '"\itt = \rm', num2str(time(ii)),'s",'];
%         end
%         lgd_str = ['legend(',lgd_str,'"Location","best")'];
%         eval(lgd_str);
        %----------------------------------------------------
        set(gca, 'FontName', 'Arial', 'FontSize', 14);
        box on; grid off;
        hold off;
    else
        error('The <type> must be "surf", "contour" or "curve".')
    end
end