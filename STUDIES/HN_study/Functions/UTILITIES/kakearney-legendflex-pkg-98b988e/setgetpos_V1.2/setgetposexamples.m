function setgetposexamples(n)
% SETGETPOSEXAMPLES Launch SETPOS/GETPOS examples.
%   SETGETPOSEXAMPLES(N) launches the Nth example for setpos/getpos 
%   commands. Edit SETGETPOSDEMO.M for more informations about 
%   each example.
%

%   This demo file will be improved in the future (...I hope)

%   Author: Jérôme Briot, Matlab 6.1.0.450 (R12.1)
%   Contact: dutmatlab@yahoo.fr
%   Revision: 1.0 (21-Feb-2007)
%   Comments:
%

% Check the number of input arguments

narginchk(0,1);


if nargin==0
    n=1;   
end

switch n
    
	case 1%       Create a figure with "Position" property sets to 
          %       [0.25{normalized} 100{pixels} 0.5{normalized} 300{pixels}] :

			h=figure;
			setpos(h,'0.25nz 100px 0.5nz 300px');   
            
            disp('To create a figure with "Position" property sets to')
            disp('[0.25{normalized} 100{pixels} 0.5{normalized} 300{pixels}]')
	case 2%       Only set the width of a figure to 0.5{normalized} :

			h=figure('units','pixels','position',[0 200 200 400]);
			pause(.5)
			setpos(h,'# # .5nz #');

    case 3%       Add 200{pixels} to the width of a figure :

			h=figure('units','normalized','position',[.1 .1 .5 .8]);
			pause(.5)
			setpos(h,'# # +200px #');

    case 4%       Use SETPOS as SET(H,'Position',...)   

			h=figure('units','pixels');
			setpos(h,[100 100 300 200]);


    case 5%       TSet the position of a pushbutton according to the current 
          %       axes position (instead of the figure parent)

			figure
			axes
			u(1)=uicontrol('string','(0,0) gcf');
			setpos(u(1),'0 0 60px 40px')
			u(2)=uicontrol('string','(0,0) gca');
			setpos(u(2),'0 0 60px 40px',gca)
   
    case 6%       Create a uniformly-spaced group of buttons

			figure
			u(1)=uicontrol;
			for n=2:5
			    u(n)=uicontrol;
			    setpos(u(n),'# 30px # #',u(n-1));
			end
            
            
    case 7%       Get the Left&Bottom position in {Pixels} and the
          %       Width&Height position in {Points} of a figure object :

			h=figure('units','normalized','position',[.1 .1 .5 .8]);
			pos=getpos(h,'px px pt pt')

    case 8%       Get the default "Position" of the figure object in all units

			h=figure;
			pos=[getpos(h,'px')
                 getpos(h,'nz')
                 getpos(h,'in')
                 getpos(h,'cm')
                 getpos(h,'pt')
                 getpos(h,'ch')]

    case 9%       Only get the width in {Normalized} of a figure

			h=figure('units','pixels','position',[100 100 200 400]);
			pos=getpos(h,'# # nz #')
			%or
			pos=getpos(h,'# # nz #','compact')

    case 10%      Get the position of one button to another one

			figure
			u(1)=uicontrol('units','pixels','position',[50 50 100 30]);
			u(2)=uicontrol('units','pixels','position',[200 150 100 30]);
			pos=getpos(u(1),'px',u(2))
			pos=getpos(u(1),'# px # px',u(2),'compact')

end