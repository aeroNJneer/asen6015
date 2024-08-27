function plotNeesNis(t, NEES, NIS, alpha, numTrials, title, xmax,neesmax, nismax)
    if iscell(NEES)
        NEES_bar = mean(cell2mat(NEES),1);
        NIS_bar = mean(cell2mat(NIS),1);
    else 
        NEES_bar = mean(NEES, 1);
        NIS_bar = mean(NIS, 1);
    end

    
    r1nees = chi2inv(alpha/2,4*numTrials)/numTrials * ones(size(t(2:end)));
    r2nees = chi2inv(1-alpha/2,4*numTrials)/numTrials * ones(size(t(2:end)));
    r1nis = chi2inv(alpha/2,3*numTrials)/numTrials * ones(size(t(2:end)));
    r2nis = chi2inv(1-alpha/2,3*numTrials)/numTrials * ones(size(t(2:end)));

    figure
    subplot(2,1,1)
    hold on
    sgtitle(title)
    scatter(t(2:end),NEES_bar)
    plot(t(2:end), r1nees, t(2:end),r2nees)
    xlabel('Time (s)'); xlim([0 xmax])
    ylabel('NEES statistic'); ylim([0 neesmax])
    
    subplot(2,1,2)
    hold on
    scatter(t(2:end),NIS_bar)
    plot(t(2:end), r1nis,t(2:end), r2nis);
    xlabel('Time (s)'); xlim([0 xmax])
    ylabel('NIS statistic'); ylim([0 nismax])

end