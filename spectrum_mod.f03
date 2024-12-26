      module spectrum_mod
!
!     This module is used for building spectra.
!
!     Hrant P. Hratchian, 2024.
!     hhratchian@umcerced.edu
!     University of California, Merced
!
!
      USE iso_fortran_env
      USE MQC_General
      USE MQC_Gaussian
!
!     Set up type/class definitions.
!
!     SpectrumData Type
      type::spectrumData
        logical::initialized=.false.
        integer(kind=int64)::nPeaks,nPeaksAdded=-1
        real(kind=real64),dimension(:),allocatable::peakPositions,peakIntensities,peakFWHMs
        character(len=128)::unitsPeakPositions,unitsPeakIntensities
        logical,dimension(:),allocatable::peakFilled
      Contains
        procedure,public::init     => spectrumData_init
        procedure,public::addPeak  => spectrumData_addPeak
      end type SpectrumData
!
!     Module global variables.
!
      integer(kind=int64),parameter::iOut=6,fChkUnit=25,simOut=26
      real(kind=real64),private::defaultPeakIntensity=1.0,defaultPeakFWHM=0.01
!
!
      CONTAINS
!
!
!PROCEDURE spectrumData_init
      subroutine spectrumData_init(spectrum,nPeaks,unitsPeakPositions,  &
        unitsPeakIntensities)
!
!     This routine initializes a spectrumData object.
!
!
      implicit none
      class(spectrumData)::spectrum
      integer(kind=int64),intent(in)::nPeaks
      character(len=*),intent(in),optional::unitsPeakPositions,unitsPeakIntensities
!
!     Initialize the spectrumData object. If the units aren't sent, default
!     those in the object. The last step here is to set the init flag in the
!     object to TRUE.
!
      if(nPeaks.le.0) call mqc_error('Constructor for spectrumData object sent nPeaks<=0.')
      spectrum%nPeaks = nPeaks
      if(Allocated(spectrum%peakPositions)) DeAllocate(spectrum%peakPositions)
      if(Allocated(spectrum%peakIntensities)) DeAllocate(spectrum%peakIntensities)
      if(Allocated(spectrum%peakFWHMs)) DeAllocate(spectrum%peakFWHMs)
      if(Allocated(spectrum%peakFilled)) DeAllocate(spectrum%peakFilled)
      Allocate(spectrum%peakPositions(nPeaks),  &
        spectrum%peakIntensities(nPeaks),spectrum%peakFWHMs(nPeaks),  &
        spectrum%peakFilled(nPeaks))
      spectrum%peakFilled(:) = .false.
      if(Present(unitsPeakPositions)) then
        spectrum%unitsPeakPositions = TRIM(unitsPeakPositions)
      else
        spectrum%unitsPeakPositions = 'cm-1'
      endIf
      if(Present(unitsPeakIntensities)) then
        spectrum%unitsPeakIntensities = TRIM(unitsPeakIntensities)
      else
        spectrum%unitsPeakIntensities = 'unknown'
      endIf
      spectrum%nPeaksAdded = 0
      spectrum%initialized = .true.
!
      return
      end subroutine spectrumData_init
!
!
!PROCEDURE spectrumData_addPeak
      subroutine spectrumData_addPeak(spectrum,peakPosition,peakIntensity,peakFWHM)
!
!     This routine adds a peak to a spectrumData object (spectrum). The peak
!     position is sent as dummy argument peakPosition, which is a required
!     argument. The peak intensity and peak full-width-at-half-max are sent as
!     dummy arguments peakIntensity and peakFWHM; these are optional arguments.
!     If peakIntensity is NOT sent, it is set to the module default value. If
!     peakFWHM is NOT sent, it is set to the module default value.
!
!
      implicit none
      class(spectrumData)::spectrum
      real(kind=real64),intent(in)::peakPosition
      real(kind=real64),intent(in),optional::peakIntensity,peakFWHM
!
!     Add the peak data to spectrum.
!
      if(.not.spectrum%initialized) call mqc_error('Trying add a peak to an unitialized spectrumData object.')
      if(spectrum%nPeaksAdded.ge.spectrum%nPeaks) call mqc_error('Trying add a peak to a full spectrumData object.')
      spectrum%nPeaksAdded = spectrum%nPeaksAdded+1
      spectrum%peakPositions(spectrum%nPeaksAdded) = peakPosition
      if(Present(peakIntensity)) then
        spectrum%peakIntensities(spectrum%nPeaksAdded) = peakIntensity
      else
        spectrum%peakIntensities(spectrum%nPeaksAdded) = defaultPeakIntensity
      endIf
      if(Present(peakFWHM)) then
        spectrum%peakFWHMs(spectrum%nPeaksAdded) = peakFWHM
      else
        spectrum%peakFWHMs(spectrum%nPeaksAdded) = defaultPeakFWHM
      endIf
      spectrum%peakFilled(spectrum%nPeaksAdded) = .true.
!
      return
      end subroutine spectrumData_addPeak
!
!
!PROCEDURE spectrumData_plot_value
      function spectrumData_plot_value(spectrum,x,scaleFactor) result(plot)
!
!     This routine is used to evaluate the value of the spectrum plot at
!     position x. The optional dummy argument scaleValue is used to scale the
!     plot value.
!
!
      implicit none
      class(spectrumData),intent(in)::spectrum
      real(kind=real64),intent(in)::x
      real(kind=real64),intent(in),optional::scaleFactor
      real(kind=real64)::plot
      integer(kind=int64)::i
      real(kind=real64)::prefactorDenominator,sigma,sigmaDenominator,myScaleFactor
!
!     Figure out the value of myScaleFactor.
!
      myScaleFactor = 1.0
      if(Present(scaleFactor)) myScaleFactor = scaleFactor
!
!     Loop over the peaks in spectrum and build the plot value at x.
!
      prefactorDenominator = SQRT(2.0*Pi)
      sigmaDenominator = float(2)*SQRT(float(2)*log(float(2)))
      plot = float(0)
      do i = 1,spectrum%nPeaksAdded
        sigma = spectrum%peakFWHMs(i)/sigmaDenominator
        plot = plot + myScaleFactor*  &
          (spectrum%peakIntensities(i)/(prefactorDenominator*sigma))*  &
          Exp(-((x-spectrum%peakPositions(i))**2)/(float(2)*sigma**2))
      endDo
!
      return
      end function spectrumData_plot_value
!
!
!
      end module spectrum_mod
