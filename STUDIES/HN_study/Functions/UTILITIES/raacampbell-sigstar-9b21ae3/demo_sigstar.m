function demo_sigstar
% demo function for sigstar
%
% function demo_sigstar
%
% No input or output arguments.
%
%
% 


rows=2;
cols=3;

clf



subplot(rows,cols,1)
barData=[2,1.2,1.9,1,2,3];
H=bar(barData);
set(H,'FaceColor',[0.5,1,0.5])
groups={[3,6],[1,2],[5,4]};

%overlay random data
hold on
for ii=1:length(barData)
	r=randn(1,10)*0.1+barData(ii);
	p(ii)=plot(zeros(size(r))+ii, r ,'b.');
end
H=sigstar(groups,[0.001,0.05,0.04]);
%now jitter the points
for ii=1:length(p)
	x=get(p(ii),'XData');
	set(p(ii),'XData',x+randn(size(x))*0.05);
end
xlim([0.5,6.5])



subplot(rows,cols,2)
barData=[10,7,6,5];
H=bar(barData);
set(gca,'XTickLabel',{'X','a','b','c'})
set(H,'FaceColor',[0.5,1,0.5])
groups={{'X','c'},...
 		{1,'b'},... %note you can mix and match notions
		{'X','a'}};


%overlay random data
hold on
for ii=1:length(barData)
	r=randn(1,10)*0.4+barData(ii);
	p(ii)=plot(zeros(size(r))+ii, r ,'k.');
end

%H=sigstar(groups,[0.001,0.05,0.04]);

%now jitter the points
for ii=1:length(p)
	x=get(p(ii),'XData');
	set(p(ii),'XData',x+randn(size(x))*0.1);
end


H=sigstar(groups,[0.0001,0.005,0.04]);
set(H,'Color','b')

%Extend the last arms down:
Y=get(H(1,1),'YData'); Y(4)=7.5; set(H(1,1),'YData',Y);
Y=get(H(2,1),'YData'); Y(4)=6.5; set(H(2,1),'YData',Y);
Y=get(H(3,1),'YData'); Y(4)=5.5; set(H(3,1),'YData',Y);
ylim([0,13])
xlim([0.5,4.5])



subplot(rows,cols,3)
%A box plot
d=randn([20,3]);
d(:,2)=d(:,2)+2;
boxplot(d)
ylim([-3,6])
H=sigstar({[1,2],[2,3]},[]);


subplot(rows,cols,4)
%A line plot
x=1:12;
y=randn(size(x));
y(5)=y(5)-5;
y(6)=y(6)+5;

plot(x,y,'o-r','MarkerFaceColor',[1,0.5,0.5])
sigstar({[5,6]},[0.05]);
xlim([1,12])



subplot(rows,cols,5)
%Grouped bar chart
y = [2 2 3; 2 5 6; 2 8 9; 6 13 15];
bar(y)
sigstar({[1,4],[2,4],[3,4]},[0.05,0.05,0.05]);
xlim([0.5,4.5])



if exist('notBoxPlot')
	subplot(rows,cols,6)
	R=randn(15,3);
	R(:,2)=	R(:,2)+2;
	h=notBoxPlot(R);
	set([h.data],'MarkerSize',3)
	sigstar({[1,2],[3,2]},[0.05,0.05]);	
	box on
end


