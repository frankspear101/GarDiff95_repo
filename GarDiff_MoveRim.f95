	PROGRAM GarDiff_MoveRim95
! 
!     FINITE DIFFERENCE PROGRAM FOR CALCULATING DIFFUSION PROFILE IN
!     INITIALLY HOMOGENEOUS quartz
!     (1) this version assumes spherical geometry
!     (2) the grid spacing is fixed
!     (3)  SURFACE COMPOSITION is variable and modeled as a linear function of core and rim
!     (4) DIFFUSION COEFFICIENT IS A FUNCTION OF TIME(TEMPERATURE).
! 
	USE AWE_INTERFACES
	implicit none
	TYPE(AWE_Canvas) :: XYPlot
	TYPE(AWE_CanvasPen) :: pen
	Type(AWE_CanvasBrush) :: brush

	Include "PlotStuff.inc"
	include "GarDiff_MoveRim.inc"
! c*********************************************************
!     local variables
	REAL*8 RJoules,Tinterval,dTemp,Dta2Stop,Dta2
	real*8 change,geometry,ppmStart
	real*8 kCarlson,aZeroCarlson
	real*8 tempToList
	REAL*4 X,Y
	real*8 xPoint(500),yPoint(500),slope(500),intercept(500)
	real*4 Tintv4,Tt(20,5)
	CHARACTER*32 modelTitle,fileIn,titleline
	integer*4 i,k,iopt,iup,iout,iok,inttintv,Lgeometry,ioutput,stopOnDta2,gridStart
	integer*4 newProblem,numTt,numDone,numPoints,initialType,coreBoundary
	real*4 rimPositionSlope,rimPositionStart,rimCompSlope,rimComp,rimCompStart,rimPosition,rimPositionChange
	integer*4 rimGridStart
	integer*4 status,myColor
	character*32 FLNM,coefName
	character*64 aline
	data RJoules/8.3144D0/
! c*******************************************************************
! c**********Start program here***********************
! c---------Block to initialize windows and menus-----------------
	PSFILE = 0
	FLNM = 'Zoning.out'

	OPEN(12,FILE = 'OutPut',ACCESS = 'window, 400, 600')
	ColorFileName = "GarDiffColorFile.txt"
	call ReadColorFile
	currentColor = 1 		! black
	myColor = 1
	RVAL = 0.1
!	set some plotting defaults
	Pltitle = 'Difusion'
	XLAB  = 'Distance'
	xor = 50
	xmin = 0.
	xmax = 1000.
	xlen = 20.		! units are cm
	nxstep = 10
	nxdec = 0
	YLAB = 'Composition'
	yor = 50
	ymin = 0
	ymax = 1.
	ylen = 15
	nystep = 10
	nydec = 1

1     	continue
	write(*,*)' Main menu options'
	write(*,*)' 1 = set up new problem'
	write(*,*)' 2 = Draw plot '
	write(*,*)' 3 = compute profiles'
	write(*,*)' 4 = plot results'
	write(*,*)' 5 = Save graphic as Illustrator file'
	read(*,*)(iopt)
2     	continue
	if(iopt.eq.4)then
		call tcplot(XYPlot)
		go to 1
		endif
	if(iopt.eq.5)then
		PSFileName = Trim(filein)//'.ps'        	! default ps file name
		call psopcl(5)    				! save old PostScript scratch file
		go to 1
		endif
! c------------Option 1: set up problem----------------------------     
	if(iopt.eq.1)then
	tempToList = 0
	newProblem = 0
	open(2,file='',status='old', IOSTAT=status)
	if(status.ne.0)then
		call FSS_Alert('ALERT','Problem opening input file ... abort')
		write(*,*)'Problem opening input file'
		go to 1
		endif
	INQUIRE (2, NAME=filein)
	read(2,*)TitleLine
	read(2,*)               ! Stars
	read(2,'(a)')Coefname
	read(2,*)              ! names of elements
	write(*,*)'Reading diffusion coefficients'
	read(2,*)(dZero)
	read(2,*)(eActiv)
	read(2,*)(vActiv)
	read(2,*)kCarlson,aZeroCarlson
!     DCONV IS CONVERSION FACTOR FOR Cm2/sec to u2/MILLION YEARS
	SECMA=3600.D0*24.D0*365.D0*1.D06
	DCONV=SECMA*1.D08
!     DZERO IN u2/Ma
	Dzero=Dzero*DCONV
	read(2,*)               ! Stars

!	Initial conditions
	write(*,*)'Reading initial conditions'
	read(2,*)modelTitle			! model title
	read(2,*)qtzRadius			! quartz radius
	read(2,*)dGrid				! grid spacing
	read(2,*)Lgeometry
	read(2,*)coreBoundary			! 0 = mirror, 1 = fixed
! set initial grid
	nGrid = int((qtzRadius+.5)/dGrid)
	write(*,*)'nGrid = ',nGrid
	read(2,*)         ! dashes
	read(2,*)				! dummy title line
	write(*,*)'Reading initial profile data'
	read(2,*)initialType

	select case (initialType)
	case(0)			! homogeneous composition
		read(2,*)ppmStart				! if ppmStart<0 then we are reading the initial profile from the input file - otherwise it is a homogeneous phase
		qtzRad(0)=0
		qtzComp(0) = ppmStart
		do 17 i = 1,maxdim		! fill up the entire array so there are no zero entries
		qtzRad(i) = i*dGrid
		qtzComp(i) = ppmStart
17		continue
	case(1)			! read initial composition (must be even grid spaces)
		read(2,*)ngrid		! note that we are resetting ngrid from what was calculated above
		do 18 i = 0,ngrid
		read(2,*)qtzRad(i),qtzComp(i)
		write(12,*)i,qtzRad(i),qtzComp(i)
		!write(*,*)i,qtzrad(i),qtzComp(i)
18		continue		
	case(2)			! read points and interpolate for starting conditions
		read(2,*)numPoints
		i = 0
20		continue
!		do 20 i = 1,numPoints
		read(2,'(A64)')aline
		write(*,*)aline
		if(aline(1:5).eq.'-----')go to 21
		i = i + 1
		read(aline,*)xPoint(i),yPoint(i)
		go to 20
21		continue
		Write(*,*)'Points read successfully'
		numPoints = i
		do 22 i = 2,numPoints
		slope(i) = (yPoint(i)-yPoint(i-1))/(xPoint(i)-xPoint(i-1))
		intercept(i) = yPoint(i) - slope(i)*xPoint(i)
22		continue
		do 24 i = 1,maxdim		! fill up the entire array so there are no zero entries
		qtzRad(i) = i*dGrid
24		continue
		k = 2
		do 25 i = 0,ngrid
		if(qtzRad(i).ge.xPoint(k))then
			k=k+1
			if(k.gt.numPoints)then
				write(*,*)'k,i,qtzRad(i),xPoint(k) ',k,i,qtzRad(i),xPoint(k)
				write(*,*)' Not enough points defined to interpolate the entire grid'
				write(*,*)' Add some more points or distance to the interpolation points'
				pause 'hit return to continue'
				endif
			endif
		qtzComp(i) = intercept(k) + slope(k)*qtzRad(i)
		write(*,*)'i,qtzRad(i),qtzComp(i) ',i,qtzRad(i),qtzComp(i)
25		continue

	case default
	end select

	read(2,*)         ! dashes
	write(*,*)'Reading T-t information'
	read(2,*)numTt			! number of T-t points
	write(*,*)'number T-t points',numTt
	read(2,*)aline				! label	
	do 15 i = 1,numTt
	read(2,*)Tt(i,1),Tt(i,3),Tt(i,4),Tt(i,5)
	write(*,*)'Tt Point# = ',i,Tt(i,1),Tt(i,3),Tt(i,4),Tt(i,5)
	!Tt(i,1) = T; Tt(i,2) = time; Tt(i,3) = cooling rate
	!Tt(i,4) = Comp; Tt(i,5) = radius
15	continue
	read(2,*)		! line of dashes
	write(*,*)'Reading plotter information'
	READ(2,'(a)')Pltitle
	READ(2,'(a)')XLAB
	read(2,*)xor,xmin,xmax,xlen,nxstep,nxdec
	READ(2,'(a)')YLAB
	read(2,*)yor,ymin,ymax,ylen,nystep,nydec
	write(12,*)xor,xmin,xmax,xlen,nxstep,nxdec
	write(12,*)yor,ymin,ymax,ylen,nystep,nydec
	read(2,*)		! line of dashes
	write(*,*)'Reading measured profile'
	read(2,*)		! header
	read(2,*)		! header
	i = 0
14	continue
	i = i+1	
	read(2,*,end = 13)dataToFit(i,1),dataToFit(i,2)
	go to 14
13	continue
	numDataToFit = i-1
	close(2)
	Tt(1,3)=0		! cooling rate for initial conditions must be zero
	Tt(1,2) = 0		! initial time must be zero
	do 16 i=2,numTt
	Tt(i,2) = dabs((Tt(i-1,1)-Tt(i,1))/Tt(i,3)) ! calculate the time at each temperature point
16	continue
	Rjoules = 8.314		! gas constant in joules
	geometry=Dfloat(Lgeometry)
	if(Lgeometry.lt.0.or.Lgeometry.gt.2)go to 1
	newProblem = 1
	if(Lgeometry.eq.0)then
		gridstart = 1			! 1 = core - linear inner boundary starts at core
		else
		gridStart = 2			! cylinder and sphere inner boundary starts at 1 in from core
		endif

	go to 1
	endif
! c------------Option 2: Plotting ------------------------------
	if(iopt.eq.2)then
!	xmin = 400
!	xmax = 1000.
!	ymin = 0.
!	ymax = 12.
!	xor = 50
!	xlen = 20
!	nxstep = 12
!	nxdec = 0
!	xlab = 'T_C'
!	yor = 50
!	ylen = 20
!	nystep = 24
!	nydec = 0
!	ylab = 'P_kb'
!	pltitle = 'GeoBaRamanTry'
	
	myColor = 1
	CurrentColor = 1	! plot the axes using black
	call SetPlot(iok)
	if(iok.eq.1)go to 1  ! cancel button was selected
	CALL PlotAxes(XYPlot)
	go to 1
	endif
! c------------Option 3: Calculate ------------------------------
	if (iopt.eq.3)then
!     ASSIGN INITAL VALUES TO PARAMETERS
	Ti = Tt(1,1)
	Tfinal = Tt(2,1)
	coolrt = Tt(2,3)
	rimPositionSlope = (Tt(2,5)-Tt(1,5)) /(Tfinal-Ti)
	rimPositionStart = Tt(1,5)
	rimPosition = rimPositionStart
	rimGridStart = ngrid
	rimCompSlope = (Tt(2,4)-Tt(1,4))/(Tfinal-Ti)		! this sets the change rate of the rim comp
	!rimCompStart = ppmStart 
	rimCompStart = qtzComp(ngrid)				! this sets the initial rim comp
	rimComp = rimCompStart
	numDone = 2			! number of intervals done
	PBAR=5000.D0            ! nominal pressure for T calculation
	TC = Ti
	TK=TC+273.15
	PBar=5000
	iup=0

	myColor = 1
	call PickAWEColor(myColor)
	CurrentColor = myColor				! CurrentColor is in PlotStuff.inc common block

	DO 122 I=0,ngrid
	Y=qtzComp(i)  			! note qtzcomp and qtzrad are real*8 and plot requires real*4
	x=qtzRad(i)
	CALL PLOT(XYPlot,x,y,iup)
	CALL pplot(x,y,IUP)
	iup=1
122 	continue
	CALL ppenup
	iup = 0
	DO 72 i=1,numDataToFit
	x=dataToFit(i,1)
	Y=dataToFit(i,2)  			! note qtzcomp and qtzrad are real*8 and plot requires real*4
	CALL PLOT(XYPlot,x,y,iup)
	CALL pplot(x,y,IUP)
	iup=1
72 	continue
	CALL ppenup
	!     CALCULATE DIFFUSION COEFFICIENT
	DC = Dzero*exp(-eActiv/(Rjoules*TK) + kCarlson*AzeroCarlson)
	!     OUTPUT initial data
	iout=12
117   continue
	write(iout,*)'Ti in quartz Diffusion output'
	write(iout,*)
	write(iout,*)modelTitle
	write(iout,*)
	write(iout,1318)Coefname
	write(iout,*)'dZero',Dzero
	write(iout,*)'eActiv',eActiv
1318 	format(' Diffusion coefficients: ',(a80))
	if(LGeometry.eq.0)then
		write(iout,*)' Geometry: Linear'
		elseif(LGeometry.eq.1)then
		write(iout,*)' Geometry: Cylindrical'
		elseif(LGeometry.eq.2)then
		write(iout,*)' Geometry: Spherical'
		endif            
	!      write(iout,*)
	write(iout,*)'Quartz: Initial data'
	write(iout,*)'RADIUS                 =',qtzRADIUS
	write(iout,*)'# grid points          =',ngrid
	write(iout,*)' T     time         coolRate     RimComp    RimPos'
	do 318 i=1,numTt
	write(iout,*)Tt(i,1),Tt(i,2),Tt(i,3),Tt(i,4),Tt(i,5)
318	continue
	write(iout,*)' R value               =',RVAL
	!       pause 'Pausing...'
	iout=12
	tMyrsTotal=0.D0
	tMyrsInterval=0.D0
	tSec=0.D0
	index=0
	ioutput = 2
	!     BEGIN MAIN LOOP
	stopOnDta2 = 0
	Dta2 = 0.0d0
349   continue
	tcount=0.D0
344   continue
	write(*,*)
	write(*,*)
	write(*,*)'Current T-t interval  = ',numdone
	write(*,*)'Tfinal,coolrate       = ',tfinal,coolrt
	write(*,*)'Current T             = ',TC
	write(*,*)'Total time          = ',tMyrsTotal
	write(*,*)'Current time interval = ',tMyrsInterval
	write(*,*)'Time interval to stop = ',Tt(numdone,2)
	if(ioutput.eq.0)write(*,*)' Output option = No Output'
	if(ioutput.eq.1)write(*,*)' Output option = Index T only Output'
	if(ioutput.eq.2)then
	write(*,*)' Output option = Full Output'
	write(*,*)'dTime =    ',dtime
	write(*,*)'DiffCoef = ',DC
	write(*,*)'Dt/a**2Stop =  ',stopOnDta2,Dta2Stop
	write(*,*)'Dt/a**2 =  ',Dta2
	endif
	write(*,*)'-----------------'
	write(*,*)' Starting time step'
	write(*,*)' Step number = ',numDone
	write(*,*)Tt(numDone,1),Tt(numDone,2),Tt(numDone,3),Tt(numDone,4),Tt(numDone,5)
	write(*,*)'Specify time interval (in Ma) to next plot'
	write(*,*)' Single step =  0'
	write(*,*)' List grid   = -1'
	write(*,*)' Change Rval = -2'
	write(*,*)' Abort       = -3'
	write(*,*)' Set output  = -4'
	write(*,*)' Write debug = -5'
	write(*,*)' Stop On Dt/a^2 = -6'
	read(*,*)(tintv4)
	tintv=tintv4
	Tinterval=TC            ! used in routine writeout as a 1 degree counter
345   continue
	inttintv=int(tintv4)
	if(inttintv.eq.-6)then
		stopOnDta2 = 1
		write(*,*)'Input value of Dt/a**2 to stop on (e.g. 0.01 to .5)'
		read(*,*)Dta2Stop
		tintv4 = 100		! set to a large value so we can get to Dta2Stop
		tintv=tintv4
		Tinterval=TC            ! used in routine writeout as a 1 degree counter
		go to 350
		endif
	if(inttintv.eq.-5)then
		call listem
		write(*,*)'Ti,Tfinal,coolRt',Ti,Tfinal,coolRt
		write(*,*)'rimPositionStart,rimPosition,rimPositionSlope',rimPositionStart,rimPosition,rimPositionSlope
		write(*,*)'rimGridStart,nGrid',rimGridStart,nGrid
		write(*,*)'rimCompStart,rimComp,rimCompSlope',rimCompStart,rimComp,rimCompSlope
		go to 349
		endif
	if(inttintv.eq.-4)then
		write(*,*)
		write(*,*)' 0 = no output'
		write(*,*)' 1 = index only'
		write(*,*)' 2 = full output (slowest)'
		read(*,*)(ioutput)
		go to 349
		endif
	if(inttintv.eq.-3)then
		go to 800
		endif
	if(inttintv.eq.-1)then
		call listem
		!            pause 'pausing...'
		go to 349
		endif
	if(inttintv.eq.-2)then
		write(*,*)
		write(*,*)'Rval = ',Rval
		write(*,*)'Input a new value for R:'
		read(*,*)Rval
		go to 349
		endif
	!      write(*,*)' Hit escape to abort'
350   CONTINUE
	if(abs(tempToList).gt.1)then
		write(*,*)TC
		tempToList = 0
		endif
	!	DC = Dzero*exp(-eActiv/(Rjoules*TK))
	DC = Dzero*exp(-eActiv/(Rjoules*TK) + kCarlson*AzeroCarlson)
	!     compute dtime from diffusion coefficient and R value
	!     R = ∆t*D/∆x**2
	!     so ∆t = Rval*∆x**2/D
	!     (DTIME in Ma)
	DTIME=dGrid*dGrid*RVAL/DC 
	!     compute ∆T from DTIME and cooling rate
	DTEMP = COOLRT * DTIME  
	!     set maximum ∆T at 5 C per step
	if(DTEMP.GT.5)then
		DTEMP=5.D0
		DTIME=DTEMP/COOLRT
		ENDIF
	!     Continue with calculations
351   continue
	tempToList = tempTolist + dtemp
	index=index+1  
	tMyrsTotal=tMyrsTotal+DTIME
	tMyrsInterval=tMyrsInterval+DTIME
	tSEC=tMyrsTotal*SECMA
	!     compute new Temperature and diffusion coeffieient
	TC=TC+DTEMP
	!      IF(TC.LE.Tfinal)then
	IF(tMyrsInterval.GE.Tt(numDone,2))then
		! check to see if there are more T-t intervals to work on
		if(numDone.eq.numTt)go to 800		! exit routine - all done
		! if here, then we keep going
		! new values for model
		Ti = Tt(numDone,1)
		Tfinal = Tt(numDone+1,1)
		coolrt = Tt(numDone+1,3)
		rimPositionSlope = (Tt(numDone+1,5)-Tt(numDone,5)) /(Tfinal-Ti)
		rimPositionStart = Tt(numDone,5)
		rimGridStart = ngrid
		rimCompSlope = (Tt(numDone+1,4)-Tt(numDone,4)) /(Tfinal-Ti)
		!		rimCompStart = qtzComp(ngrid)
		rimCompStart = Tt(numDone,4)
		numDone = numDone+1
		tMyrsInterval = 0
		write(*,*)'-----------------'
		write(*,*)' Starting new time step'
		write(*,*)' Step number = ',numDone
		write(*,*)Tt(numDone,1),Tt(numDone,2),Tt(numDone,3),Tt(numDone,4),Tt(numDone,5)
		endif
	TK=TC+273.0D0
	! compute Diffusion coeff
	DC = Dzero*exp(-eActiv/(Rjoules*TK) + kCarlson*AzeroCarlson)
	Dta2 = Dta2 + DC*dtime/qtzRadius**2			! use the initial radius for the length scale
	!	write(12,*)Dta2
	!	DC = Dzero*exp(-eActiv/(Rjoules*TK))


	! numdone is the model number	
	! gridRim is the grid point on the rim
	! determine the rim grid point (always round down = truncate)

	rimComp = rimCompStart + rimCompSlope * (TC - Ti)
	rimPositionChange = rimPositionSlope * (TC - Ti)
	rimPosition = rimPositionStart + rimPositionChange
	!	rimGridChange = int(rimPositionChange / dGrid)
	!	rimGrid = rimGridStart + rimGridChange
	!	nGrid = rimGrid	
	nGrid = rimPosition/dGrid
	!	write(*,*)'nGrid =',rimgridStart,ngrid,rimcompStart,rimcomp
	!     CALCULATE RIM COMPOSITION OF quartz
	!	qtzRimppm = 10**(slopeppm*(1./(TC+273)) + intcptppm)
	!	qtzComp(0) = qtzRimppm
	qtzComp(nGrid) = rimComp


	!     CALCULATE NEW CONCENTRATION PROFILE
	!     using an explicit finite difference formulation
	!	qtzComp(ngrid+1) = qtzComp(ngrid-1)		!symmetrical condition at quartz core
	select case (coreBoundary)
		case(0)				! mirror boundary
			qtzComp(0) = qtzComp(2)		!symmetrical condition at core
		case(1)
			!		qtzComp(0) = qtzComp(1) + (qtzComp(1)-qtzComp(2))		!fixed core composition - make curvature = 0
			qtzComp(0) = 2.0d0*qtzComp(1) - qtzComp(2)		!fixed core composition - make curvature = 0
		case default
		end select
	do 360 i=1,NGRID-1		! change to nGrid-1 because the rim comp is fixed above
	curve(i)=(qtzcomp(i+1)-2.*qtzcomp(i) + qtzcomp(i-1))/dgrid**2
	!      gradient(i)=(qtzcomp(i+1) - qtzcomp(i-1))/(2*dgrid)
	gradient(i)=(qtzcomp(i+1) - qtzcomp(i))/(dgrid)
360   continue
	if(Lgeometry.eq.2)then		! spherical geometry
		qtzcomp(1) = qtzcomp(1) + (6.0d0*(qtzcomp(2)-qtzcomp(1))/dgrid**2)*DC*dtime
		go to 362
		endif
	if(Lgeometry.eq.1)then		! cylindrical geometry
		qtzcomp(1) = qtzcomp(1) + (4.0d0*(qtzcomp(2)-qtzcomp(1))/dgrid**2)*DC*dtime
		go to 362
		endif

362	continue
	do 361 i=gridstart,ngrid-1		! change to nGrid-1 because the rim comp is fixed above
	!     Geometry=0(linear), =1(cylindrical) = 2(spherical)
	change = curve(i) + (geometry/qtzrad(i))*gradient(i)         ! the geometry factor accommodates linear, cylindrical and spherical geometries
	! check to be sure this is the correct formulation      
	qtzcomp(i)=qtzcomp(i) + change*DC*dtime 
361   continue


	!	check time interval and plot, if required
	if(stopOnDta2.eq.1)then
		if(Dta2.gt.Dta2stop)then
			!mcolor = 4
			call PlotGrid(XYPlot,mycolor)
			go to 349
			endif
		endif
	tcount=tcount+dtime
	if(tcount+dtime.lt.tintv)go to 350
	if(tcount.ge.tintv.or.inttintv.eq.0)then
		!icolor = 4
		call plotGrid(XYPlot,mycolor)			! (4) is the color = red
		go to 349
		endif

	if(tcount+dtime.gt.tintv)then
		dtime=tintv-tcount
		!           compute ∆T from new DTIME and cooling rate
		DTEMP = COOLRT * DTIME  
		go to 351
		endif
800   continue
	!icolor = 5
	call PlotGrid(XYPlot,mycolor)
	write(*,*)
	write(*,*)' Problem is done'
	write(*,*)
	write(*,*)'Total time          = ',tMyrsTotal
	write(*,*)
	iout=12
211   continue
	pause 'All done -- press return to continue'
	go to 1
	endif
	go to 1
	END

! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	subroutine PlotGrid(XYPlot,mycolor)
	USE AWE_INTERFACES
	implicit none
	TYPE(AWE_Canvas) :: XYPlot
	TYPE(AWE_CanvasPen) :: pen
	Type(AWE_CanvasBrush) :: brush
	integer*4 i,iup,myColor
	real*4 x,y
	Include 'PlotStuff.inc'
	include "GarDiff_MoveRim.inc"
	! c*******************************************************

	!	call scolor(4)
!	call scolor(icolor)
!	myColor = 1		default is to use the last color (= mycolor)
	call PickAWEColor(myColor)
	CurrentColor = myColor				! CurrentColor is in PlotStuff.inc common block
	iup=0
	DO 10 i=0,ngrid
	Y=qtzcomp(i)
	x=qtzrad(i)
	CALL PLOT(XYPlot,x,Y,iup)
	CALL pplot(x,y,IUP)
	iup=1
	10	continue
	!CurrentColor = 1		! black
	CALL ppenup
	return
	end

! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	subroutine listem
	implicit none
	integer*4 i
	include "GarDiff_MoveRim.inc"
	! c*******************************************************
	write(*,*)
	write(*,*)
	write(*,*)
	write(*,*)' i        delx         radius         conc'
	!      do 10 i=0,ngrid
	do 10 i=0,ngrid
	write(12,100)i,dGRId,Qtzrad(i),QtzCOMP(i)
100   	format (i8,4f12.6)
10    	continue
	return
	end
! c***********************************************************************
! c***********************************************************************
	subroutine tcplot(XYPlot)
	USE AWE_INTERFACES
	implicit none
	TYPE(AWE_Canvas) :: XYPlot
	TYPE(AWE_CanvasPen) :: pen
	Type(AWE_CanvasBrush) :: brush

	Include "PlotStuff.inc"
	include "GarDiff_MoveRim.inc"
	integer*4 iopt,iok,i,j,status,iup,iplot
	real*4 minc,maxc,mint,maxt,radius,x,y
	! c*********************************************************
	character*32 FLNM
	!     initialize plotter defaults
	!     XLEN=0 gives variable plot sizing
	data xor,xmin,xmax,xlen,nxstep,nxdec/70,0.,1.,0,10,1/
	data xlab/'Radius (um)'/
	data yor,ymin,ymax,ylen,nystep,nydec/70,.7,1.,6.6,8,2/
	data ylab/'Fe/Fe+Mg'/
	data PLTITLE/'garnet'/
	!------------------------------------------
	iplot=0
1     continue
	write(*,*)' Plotting menu options'
	write(*,*)' 0 = return to main program'
	write(*,*)' 1 = open and read data file'
	write(*,*)' 2 = specify plotting information'
	write(*,*)' 3 = make a plot'
	read(*,*)(iopt)
2     continue
	if (iopt.eq.0)return
! c**************option 1***********************
	if(iopt.eq.1)then
		!      iokL = stdfil(3,vref,' ',FLNM,1,'TEXT')	! 3 is getfile and set vref
		!      if (.not.iokL) go to 1
		OPEN(2,FILE='',STATUS='OLD')
		if(status.ne.0)go to 1
		INQUIRE (2, NAME=Flnm)
		do 20 i=1,35
		read(2,*)         ! read the useless stuff
20    		continue    
		read(2,32)radius
32    		format(25x,f10.3)
		do 23 i=37,60
		read(2,*)         ! read the useless stuff
23    		continue
		minT=1.e25
		maxT=-1.e25
		minC=1.e25
		maxC=-1.e25
		i=0   
25    		continue
		read(2,*,end=26)j,qtzrad(i),qtzcomp(i)    ! array element 0 is the core
		if(qtzcomp(i).gt.maxC)maxC=qtzcomp(i)
		if(qtzcomp(i).Lt.minC)minC=qtzcomp(i)
		i=i+1
		go to 25
26    		continue
		ngrid=i-1
		write(*,*)' Number of grid points = ',ngrid
		write(*,*)'Minimum Ti ppm = ',minc
		write(*,*)'Maximum Ti ppm = ',maxC
		!      pause 'pausing...'
		go to 1
		endif
! c**********option = 2*************************
	if(iopt.eq.2) then
		write(*,*)' Please specify the plot type'
		write(*,*)' 1 = quartz composition vs radius'
		read(*,*)iplot
		xlen=10 
		if(iplot.eq.1)then
		xmin=0
		xmax=radius
		ymin=.7
		ymax=1.
		nxstep=10
		nystep=15
		nxdec=0
		nydec=2
		ylab=' Fe/(Fe+Mg)'
		endif
	call SetPlot(iok)
	if(iok.eq.1)go to 1  ! cancel button was selected
	go to 1
	endif
! c**********option 3**************************
	if(iopt.eq.3)then
		if(iplot.eq.0)then
			write(*,*)'You must do option 2 before option 3'
			!            pause 'pausing...'
			go to 1
			endif
		Call PlotAxes(XYPlot)
		iup=0
		DO 337 I=0,ngrid
		x=qtzrad(i)
		if(iplot.eq.1)  Y=qtzcomp(i)  
		CALL PLOT(XYPlot,x,y,iup)
		CALL pplot(x,y,IUP)
		iup=1
337   		continue
		CALL ppenup
		go to 1
		endif
	end

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	subroutine FSS_alert(title, text)
    	use AWE_Interfaces
!	interface Alert
!	subroutine AWE_alertBox(title, text)
	character(len=*) :: title, text
	call AWE_alertBox(title,text)
	end
!	end subroutine AWE_alertBox
!	end interface
!	title is used as the title of the alert box text is the text that will be displayed in it.

