%% Sensitivity Analyisis

tiledlayout(9,3) % Change dimentions for more graphs. NB :As a large tile of graphs the produced table can only be seen clearly on a large laptop screen or a monitor

for i=1:0.25:1.5%Inital flow rate kmol/s NB: As previously determined outside of  range 0.9-1.7kmol/s is unacheivable for specified conditions of previous assesment. However avaliable to be changed for the analyist.
    
    for j=0.3:0.1:0.5 %Bed Voidage
        
        
        for k=6*10^(-3):1*10^(-3):8*10^(-3) %Dp Equivalent Diameter for assuming catalyist is a sphere m
            yo=[600 ; 450; i*0.53;i*0.43;0;i*0.02;i*0.02 ]; % Matrix of inital Conditions T(K) P(Bar) Component Flow Rates (Kmol/s) CO H M Me W  
            lspan=[0 120]; %Investigsting FBR length l for a resonable range of l 0-120m to allow desired flow rate to be reached
            [l,y]= ode45(@(l,y)FBR_Group43_28022020_V0(l,y,j,k) ,lspan,yo); % Solving the simultaneous equation relative to l and inserting into a matrix length and variables
            nexttile
            hold on
            plot(l,y(:,5))
            xlim([0 120])
            ylim([0 0.1])
            xlabel('Length of reactor (m)')
            ylabel('FM (kmol/s)') % Where FM is the outlet flow rate of methanol
            titlelayout=sprintf('F: %.2f E: %.2f Dp: %.5f', [i, j, k]); % Creating automatic table labeling relavant to a particular cycle and graph.
            title(titlelayout)
            plot([0, max(l)],[0.07766 0.07766]) %Plot of desired flow rate 0.07766Kmol/s
            hold off
           
        end
    end
 
end

 legend('FM Production rate at specified conditions','Target Flow Rate of Methanol FM') % Left out of the loop so the legend is only printed once.


            
            