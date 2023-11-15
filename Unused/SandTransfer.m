x=[60:20:260]; %set x axis ticks
y=rand(11);  %get something to plot
h1=subplot(2,1,2); %setup subplot1
plot(x,y,'-.'); %plot subplot1
box on  %  leave only x and y axes
xlim([60 260])  %setup some x axis 
set(h1,'Xtick',x) %set the x axis ticks to show only x
h1_pos = get(h1,'Position'); %get the position data for sublot1.
y2 = 10*y.^2;  %make something up for subplot2
h2=subplot(2,1,1);  %make subplot2
plot(x,10*y,'-.'); %plot subplot2
box on
set(h2,'Xcolor',[1 1 1]) %make the Y axis line white
set(h2,'Xtick',[])
xlim([60 260])  %setup some x axis 
h2_pos=get(h2,'Position'); 
set(h2,'Position',[h1_pos(1) h1_pos(2)+.1+h1_pos(4) h2_pos(3:end)])