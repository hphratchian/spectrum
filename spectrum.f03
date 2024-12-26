INCLUDE "spectrum_mod.f03"
      program spectrum
!
!     This program reads Franck-Condon progression data from a Gaussian
!     formatted checkpoint file and generates {x,y} plotting data with gaussian
!     enveloping.
!
!     Hrant P. Hratchian, 2024.
!     hhratchian@umcerced.edu
!     University of California, Merced
!
!
      USE spectrum_mod
      implicit none
      integer(kind=int64)::i,nElements,nPeaksFC,nPlotPoints
      real(kind=real64)::minPeakPosition,maxPeakPosition,  &
        plotStepSize=0.10,maxFWHM=1500.0,scaleFactor
      real(kind=real64),dimension(:),allocatable::tmpVector,fcDataLinear
      real(kind=real64),dimension(:,:),allocatable::fcDataTable,plotData
      character(len=512)::fchkFileName,tmpChar
      logical::fail=.false.,OK
      type(MQC_Gaussian_FChk_File)::myFChk
      type(spectrumData)::fcProgression
!
!     Format statements.
!
 1000 format(1x,'Program spectrum.')
 4000 format(/,1x,'Stick Spectrum Data',/,  &
        10x,'Energy (cm^-1)',15x,'Intensity',15x,'Dipole Strength')
 4100 format(5x,f20.5,5x,f20.5,5x,f20.5)
 5000 format(/,1x,'Generating the plotting data.',/,  &
        4x,'Step Size = ',f8.5,3x,'Starting Point = ',f20.5,/,  &
        27x,'Ending   Point = ',f20.5)
 5100 format(/,1x,'Spectral Pointing Data')
 5110 format(20x,'X',25x,'Y')
 5120 format(5x,f20.5,5x,f20.5)
!
!
!     Begin the program.
!
      write(iOut,1000)
      call myFChk%openFile('test.fchk',fChkUnit,OK)
      open(Unit=simOut,file='simulated.dat')
      write(iOut,*)
      write(iOut,*)' OK = ',OK
      write(iOut,*)
!
!     Read in the FC data from the FChk file.
!
      call Find_FChk_Entry(myFChk,'FCHT RAssign',OK,tmpChar,nElements,  &
        Vector_Real=fcDataLinear)
      write(iOut,*)' Hrant - After Find_FChk_Entry...'
      write(iOut,*)'           OK           = ',OK
      write(iOut,*)'           tmpChar      = ',tmpChar
      write(iOut,*)'           nElements    = ',nElements
      write(iOut,*)'           SIZE(fcData) = ',SIZE(fcDataLinear)
      write(iOut,*)
      if(Mod(nElements,3).ne.0) call mqc_error('nElements is not divisible by 3!')
      nPeaksFC = nElements/3
      Allocate(fcDataTable(3,nPeaksFC))
      fcDataTable = Reshape(fcDataLinear,[ 3,nPeaksFC ])
      call mqc_print(fcDataTable,iOut=iOut,header='This is the FC Data Table')
      write(iOut,4000)
      do i = 1,nPeaksFC
        write(iOut,4100) fcDataTable(:,i)
      endDo
!
!     Initialize the fcProgression object and fill the object with spectral
!     features read from the fchk file.
!
      call fcProgression%init(nPeaksFC,unitsPeakPositions='cm-1',unitsPeakIntensities='arbitrary')
      do i = 1,nPeaksFC
        call fcProgression%addPeak(peakPosition=fcDataTable(1,i),  &
          peakIntensity=fcDataTable(3,i),peakFWHM=maxFWHM)
      endDo
!
!     Now, evaluate the plot values for the simulated spectrum.
!

!hph+
      write(iOut,*)
      write(iOut,*)
      write(iOut,*)' Hrant - Here is the test code...'
      do i = 1,nPeaksFC
        write(iOut,*) fcDataTable(1,i),fcDataTable(3,i),spectrumData_plot_value(fcProgression,fcDataTable(1,i))
      endDo
      write(iOut,*)
      write(iOut,*)
      Allocate(tmpVector(1))
      tmpVector = MaxLoc(fcDataTable(3,:))
      i = tmpVector(1)
      DeAllocate(tmpVector)
      write(iOut,*)' The position of the largest intensity is ',i
      scaleFactor = fcDataTable(3,i)/spectrumData_plot_value(fcProgression,fcDataTable(1,i))
      write(iOut,*)' scaleFactor = ',scaleFactor
      write(iOut,*)
      write(iOut,*)
      write(iOut,*)' Hrant - Here is the test code again...'
      do i = 1,nPeaksFC
        write(iOut,*) fcDataTable(1,i),fcDataTable(3,i),  &
          spectrumData_plot_value(fcProgression,fcDataTable(1,i),scaleFactor=scaleFactor)
      endDo
      write(iOut,*)
      write(iOut,*)
!hph-

      write(iOut,*)
      write(iOut,*)
      write(iOut,*)' Hrant - minval = ',minval(fcDataTable(1,:))
      write(iOut,*)' Hrant - maxval = ',maxval(fcDataTable(1,:))
      write(iOut,*)
      write(iOut,*)
      minPeakPosition = minval(fcDataTable(1,:))-5.0*maxFWHM
      maxPeakPosition = maxval(fcDataTable(1,:))+5.0*maxFWHM
      write(iOut,5000) plotStepSize,minPeakPosition,maxPeakPosition
      nPlotPoints = Int((maxPeakPosition-minPeakPosition)/plotStepSize) + 1
      Allocate(plotData(2,nPlotPoints))
      plotData(1,1) = minPeakPosition
      plotData(2,1) = 0.0
      do i = 2,nPlotPoints
        plotData(1,i) = plotData(1,i-1)+plotStepSize
        plotData(2,i) = 0.0
      endDo
      do i = 1,nPlotPoints
        plotData(2,i) = plotData(2,i) +   &
          spectrumData_plot_value(fcProgression,plotData(1,i),scaleFactor=scaleFactor)
      endDo
      write(iOut,5100)
      write(iOut,5110)
      write(simOut,5110)
      do i = 1,nPlotPoints
        write(iOut,5120) plotData(1,i),plotData(2,i)
        write(simOut,5120) plotData(1,i),plotData(2,i)
      endDo
!
!     The end of the program.
!
  999 Continue
      close(unit=fChkUnit)
      close(unit=simOut)
      end program spectrum
