INCLUDE "spectrum_mod.f03"
      program spectrum
!
!     This program reads Franck-Condon progression data from a Gaussian
!     formatted checkpoint file and generates {x,y} plotting data with gaussian
!     enveloping. The program requires three command line arguments:
!           1. the name of the Gaussian fchk file with FC data
!           2. the name of the output text spectrum file
!           3. the name of the output text VMI data file
!
!
!     Hrant P. Hratchian, 2024, 2025.
!     hhratchian@umcerced.edu
!     University of California, Merced
!
!
      USE spectrum_mod
      implicit none
      integer(kind=int64)::i,j,nElements,nPeaksFC,nPlotPoints,  &
        nVMIPlotPoints
      real(kind=real64)::minPeakPosition,maxPeakPosition,  &
        plotStepSize=0.30,maxFWHM,scaleFactor
      real(kind=real64),dimension(:),allocatable::tmpVector,fcDataLinear
      real(kind=real64),dimension(:,:),allocatable::fcDataTable,  &
        plotData,vmiPlotData
      character(len=512)::fchkFileName,spectrumFileName,vmiFileName,  &
        tmpChar
      logical::fail=.false.,OK
      type(MQC_Gaussian_FChk_File)::myFChk
      type(spectrumData)::fcProgression
!
!     Format statements.
!
 1000 format(1x,'Program spectrum.')
 1500 format(/,1x,'FChk File:     ',A,/,  &
               1x,'Spectrum File: ',A,/,  &
               1x,'VMI File:      ',A,/)
 4000 format(/,1x,'Stick Spectrum Data',/,  &
        10x,'Energy (cm^-1)',15x,'Intensity',15x,'Dipole Strength')
 4010 format(/,1x,'Stick Spectrum Data',/,  &
        10x,'Energy (eV)',18x,'Intensity',15x,'Dipole Strength')
 4100 format(5x,f20.5,5x,f20.5,5x,f20.5)
 5000 format(/,1x,'Generating the plotting data.',/,  &
        4x,'Step Size = ',f8.5,3x,'Starting Point = ',f20.5,/,  &
        27x,'Ending   Point = ',f20.5)
 5100 format(/,1x,'Spectral Plotting Data')
 5110 format(20x,'X(cm^-1)',25x,'X(eV)',25x,'Y')
 5120 format(5x,f20.5,5x,f20.5,5x,f20.5)
 5200 format(1x,f10.5)
 5500 format(*(I0,1x))
!
!
!     Begin the program.
!
      write(iOut,1000)
      call get_command_argument(1,fchkFileName)
      call get_command_argument(2,spectrumFileName)
      call get_command_argument(3,vmiFileName)
      call myFChk%openFile(fchkFileName,fChkUnit,OK)
      open(Unit=simOut,file=spectrumFileName)
      open(Unit=vmiOut,file=vmiFileName)
      write(iOut,1500) TRIM(fchkFileName),TRIM(spectrumFileName),  &
        TRIM(vmiFileName)
!
!     Read in the FC data from the FChk file.
!
      call Find_FChk_Entry(myFChk,'FCHT RAssign',OK,tmpChar,nElements,  &
        Vector_Real=fcDataLinear)
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
!     Convert the peak energies to eV. Also, scale the intensities and set
!     maxFWHM.
!
!hph      fcDataTable(1,:) = fcDataTable(1,:)*evPHartree/cmM1PHartree
      Allocate(tmpVector(1))
      tmpVector = MaxLoc(fcDataTable(2,:))
      i = tmpVector(1)
      scaleFactor = fcDataTable(2,i)
      fcDataTable(2,:) = fcDataTable(2,:)/scaleFactor
      write(iOut,4010)
      do i = 1,nPeaksFC
        write(iOut,4100) fcDataTable(:,i)
      endDo

!hph+
!      maxFWHM = mqc_float(2)*evPHartree/cmM1PHartree
      maxFWHM = mqc_float(1)
!hph-

!
!     Initialize the fcProgression object and fill the object with spectral
!     features read from the fchk file. Note that we convert Gaussian's values
!     from wavenumbers to eV.
!
      call fcProgression%init(nPeaksFC,unitsPeakPositions='cm^-1',  &
        unitsPeakIntensities='scaled')
      do i = 1,nPeaksFC
        call fcProgression%addPeak(peakPosition=fcDataTable(1,i),  &
          peakIntensity=fcDataTable(2,i),peakFWHM=maxFWHM,peakBeta=-1.0)
      endDo
!
!     Now, evaluate the plot values for the simulated spectrum.
!
      minPeakPosition = minval(fcDataTable(1,:))-2.0*maxFWHM
      maxPeakPosition = maxval(fcDataTable(1,:))+2.0*maxFWHM
      write(iOut,5000) plotStepSize,minPeakPosition,maxPeakPosition
      nPlotPoints = Int((maxPeakPosition-minPeakPosition)/plotStepSize) + 1
      Allocate(plotData(2,nPlotPoints))
      plotData(1,1) = minPeakPosition
      plotData(2,1) = 0.0
      do i = 2,nPlotPoints
        plotData(1,i) = plotData(1,i-1)+plotStepSize
        plotData(2,i) = 0.0
      endDo
      scaleFactor = mqc_float(1)
      do i = 1,nPlotPoints
        plotData(2,i) = plotData(2,i) +   &
          spectrumData_plot_value(fcProgression,plotData(1,i),scaleFactor=scaleFactor)
      endDo
      tmpVector = MaxLoc(plotData(2,:))
      i = tmpVector(1)
      scaleFactor = plotData(2,i)
      plotData(2,:) = plotData(2,:)/scaleFactor
      if(DEBUG) then
        write(iOut,5100)
        write(iOut,5110)
      endIf
      write(simOut,5110)
      do i = 1,nPlotPoints
        if(DEBUG) write(iOut,5120) plotData(1,i),  &
          plotData(1,i)*evPHartree/cmM1PHartree,plotData(2,i)
        write(simOut,5120) plotData(1,i),  &
          plotData(1,i)*evPHartree/cmM1PHartree,plotData(2,i)
      endDo
!
!     Form a pixel matrix file for simulating a VMI image.
!
      nVMIPlotPoints = 1024
      Allocate(vmiPlotData(nVMIPlotPoints,nVMIPlotPoints))
      call spectrumData_generate_vmi_image(fcProgression,vmiPlotData,  &
        nVMIPlotPoints,nVMIPlotPoints,2.5*cmM1PHartree/evPHartree)
      if(Allocated(tmpVector)) DeAllocate(tmpVector)
      Allocate(tmpVector(2))
      tmpVector = MaxLoc(vmiPlotData)
      i = tmpVector(1)
      j = tmpVector(2)
      scaleFactor = vmiPlotData(i,j)
      vmiPlotData = vmiPlotData/scaleFactor
      vmiPlotData = mqc_float(2500)*vmiPlotData
      tmpVector = MaxLoc(vmiPlotData)
      i = tmpVector(1)
      j = tmpVector(2)
      do i = 1,nVMIPlotPoints
        write(vmiOut,5500) (Int(vmiPlotData(i,j)),j=1,nVMIPlotPoints)
      endDo
!
!     The end of the program.
!
  999 Continue
      call myFChk%closeFile()
      close(unit=simOut)
      close(unit=vmiOut)
      end program spectrum
