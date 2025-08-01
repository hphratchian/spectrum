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
        real(kind=real64),dimension(:),allocatable::peakPositions,  &
          peakIntensities,peakFWHMs,peakBeta
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
      real(kind=real64),private::defaultPeakIntensity=1.0,  &
        defaultPeakFWHM=0.01,defaultPeakBeta=0.0
!
!
      CONTAINS
!
!
!PROCEDURE spectrumData_init
      subroutine spectrumData_init(spectrum,nPeaks,unitsPeakPositions,  &
        unitsPeakIntensities,peakBeta)
!
!     This routine initializes a spectrumData object.
!
!
      implicit none
      class(spectrumData)::spectrum
      integer(kind=int64),intent(in)::nPeaks
      character(len=*),intent(in),optional::unitsPeakPositions,unitsPeakIntensities
      real(kind=real64),intent(in),optional::peakBeta
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
      if(Allocated(spectrum%peakBeta)) DeAllocate(spectrum%peakBeta)
      if(Allocated(spectrum%peakFilled)) DeAllocate(spectrum%peakFilled)
      Allocate(spectrum%peakPositions(nPeaks),  &
        spectrum%peakIntensities(nPeaks),spectrum%peakFWHMs(nPeaks),  &
        spectrum%peakBeta(nPeaks),spectrum%peakFilled(nPeaks))
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
      if(Present(peakBeta)) then
        spectrum%peakBeta(:) = peakBeta
      else
        spectrum%peakBeta(:) = mqc_float(0)
      endIf
      spectrum%nPeaksAdded = 0
      spectrum%initialized = .true.
!
      return
      end subroutine spectrumData_init
!
!
!PROCEDURE spectrumData_addPeak
      subroutine spectrumData_addPeak(spectrum,peakPosition,  &
        peakIntensity,peakFWHM,peakBeta)
!
!     This routine adds a peak to a spectrumData object (spectrum). The peak
!     position is sent as dummy argument peakPosition, which is a required
!     argument. The peak intensity, full-width-at-half-max, and anisotropy
!     parameter (beta) are sent as dummy arguments peakIntensity, peakFWHM, and
!     peakBeta; these are optional arguments. If any of these are not sent, they
!     are set to the corresponding module default value.
!
!
      implicit none
      class(spectrumData)::spectrum
      real(kind=real64),intent(in)::peakPosition
      real(kind=real64),intent(in),optional::peakIntensity,peakFWHM,  &
        peakBeta
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
      if(Present(peakBeta)) then
        spectrum%peakBeta(spectrum%nPeaksAdded) = peakBeta
      else
        spectrum%peakBeta(spectrum%nPeaksAdded) = defaultPeakBeta
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
!PROCEDURE spectrumData_generate_vmi_image
      subroutine spectrumData_generate_vmi_image(spectrum,image,M,N,  &
        photonEnergy,radiusScale)
!
!     This routine generates a simulated velocity map image (VMI) from a
!     spectrumData object. The VMI is returned as a two-dimensional array
!     (image), with pixel intensities corresponding to electron emission
!     probability mapped onto a radial and angular grid. Each peak in the
!     spectrum is modeled as a Gaussian function in radial coordinate (based on
!     its FWHM) and modulated in angle according to its anisotropy parameter
!     (beta). The radius of each feature is determined from the square root of
!     the photoelectron kinetic energy (eKE), computed as the difference between
!     the supplied photon energy and the binding energy of each peak.
!
!     If the optional dummy argument radiusScale is not provided, a dynamic
!     value is computed so that the outermost feature and its Gaussian tail are
!     scaled to fit neatly within the image bounds. The image is centered
!     automatically, and the vertical axis corresponds to θ = 0 (aligned with
!     laser polarization).
!
!
      implicit none
      class(spectrumData),intent(in)::spectrum
      real(kind=real64),intent(out)::image(:,:)
      integer(kind=int64),intent(in)::M,N
      real(kind=real64),intent(in)::photonEnergy
      real(kind=real64),intent(in),optional::radiusScale
!
!     Local variables
      integer::i,j,k,nPeaksUsed
      real(real64),parameter::pi=acos(-1.0_real64)
      real(kind=real64)::cx,cy,dx,dy,r,theta,radial_gauss,  &
        angular_factor,sigma,norm,peak_eKE,peak_r,scale,  &
        localRadiusScale,maxSqrtEKE,maxSigma,tempEKE,x,y
      real(real64),dimension(:),allocatable::eKEs,sqrtEKEs,sigmas
!
!     Initialize image
!
      image(:,:) = 0.0_real64
      cx = (real(M,real64) - 1.0) / 2.0
      cy = (real(N,real64) - 1.0) / 2.0
!
!     Determine radiusScale if not provided.
!
      if(.not.Present(radiusScale)) then
        Allocate(eKEs(spectrum%nPeaksAdded),sqrtEKEs(spectrum%nPeaksAdded),sigmas(spectrum%nPeaksAdded))
        nPeaksUsed = 0
        do k = 1,spectrum%nPeaksAdded
          if(.not.spectrum%peakFilled(k)) cycle
          tempEKE = photonEnergy - spectrum%peakPositions(k)
          if(tempEKE.le.MQC_Small) cycle
          nPeaksUsed = nPeaksUsed + 1
          eKEs(nPeaksUsed) = tempEKE
          sqrtEKEs(nPeaksUsed) = sqrt(tempEKE)
          sigmas(nPeaksUsed) = spectrum%peakFWHMs(k)/(mqc_float(2)*sqrt(mqc_float(2)*log(mqc_float(2))))
        endDo
        if(nPeaksUsed.eq.0) return
        maxSqrtEKE = maxval(sqrtEKEs(1:nPeaksUsed))
        maxSigma = maxval(sigmas(1:nPeaksUsed))
        localRadiusScale = (min(mqc_float(M),mqc_float(N))/mqc_float(2))/  &
          (maxSqrtEKE+mqc_float(3)*maxSigma)
        Deallocate(eKEs,sqrtEKEs,sigmas)
      else
        localRadiusScale = radiusScale
      endIf
      
!hph+
      write(iOut,*)
      write(iOut,*)' Hrant - localRadiusScale = ',localRadiusScale
      call mqc_print(eKEs,iOut,header='eKEs')
      call mqc_print(sigmas,iOut,header='sigmas')
!hph-

!
!     Main loop over all peaks and build the image matrix.
!
      do k = 1,spectrum%nPeaksAdded
        if(.not.spectrum%peakFilled(k)) cycle
        peak_eKE = photonEnergy - spectrum%peakPositions(k)
        if(peak_eKE.le.MQC_Small) cycle
        peak_r = localRadiusScale*sqrt(peak_eKE)
        sigma = spectrum%peakFWHMs(k)/(mqc_float(2)*sqrt(mqc_float(2)*log(mqc_float(2))))
        norm = spectrum%peakIntensities(k)/(mqc_float(2)*pi*sigma*sigma)
        do i = 1,M
          dx = mqc_float(i-1)-cx
          do j = 1,N
            dy = mqc_float(j-1)-cy
            r = sqrt(dx*dx + dy*dy)
            if(abs(r).lt.MQC_Small) then
              theta = mqc_float(0)
            else
              theta = acos(-dy/r)
            endIf
            angular_factor = mqc_float(1) + spectrum%peakBeta(k)*(1.5*cos(theta)**2-0.5)
            radial_gauss = exp(-(r-peak_r)**2/(mqc_float(2)*sigma**2))
            image(i,j) = image(i,j) + norm*angular_factor*radial_gauss
          endDo
        endDo
      endDo
!
      return
      end subroutine spectrumData_generate_vmi_image
!
!
!
      end module spectrum_mod
