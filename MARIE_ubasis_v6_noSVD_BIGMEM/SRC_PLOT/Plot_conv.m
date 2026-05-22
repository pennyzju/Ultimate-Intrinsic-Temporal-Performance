function Plot_conv(rv1,rv1_p,solver)

figure()

    semilogy(0:length(rv1)-1,rv1/rv1(1),'--k','LineWidth',2.0);
    hold on
    semilogy(0:length(rv1_p)-1,rv1_p/rv1_p(1),'-k','LineWidth',2.0);
    hold off
    xlabel('Iteration number');
    ylabel('Relative residual');
    grid on;
    legend('w/o prec.','w/  prec.')
    title(solver)