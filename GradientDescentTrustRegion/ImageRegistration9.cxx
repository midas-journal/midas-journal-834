/*=========================================================================

  This has been modified from the original included with ITK 3.20 to show
  the performance of the GradientDescentTrustRegionOptimizer.

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: ImageRegistration9.cxx,v $
  Language:  C++
  Date:      $Date: 2009-06-24 12:09:28 $
  Version:   $Revision: 1.38 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

//  Software Guide : BeginCommandLineArgs
//    INPUTS: {BrainProtonDensitySliceBorder20.png}
//    INPUTS: {BrainProtonDensitySliceR10X13Y17.png}
//    OUTPUTS: {ImageRegistration9Output.png}
//    OUTPUTS: {ImageRegistration9DifferenceBefore.png}
//    OUTPUTS: {ImageRegistration9DifferenceAfter.png}
//    1.0 300
//  Software Guide : EndCommandLineArgs

// Software Guide : BeginLatex
//
// This example illustrates the use of the \doxygen{AffineTransform}
// for performing registration in $2D$. The example code is, for the most part,
// identical to that in \ref{sec:InitializingRegistrationWithMoments}.  
// The main difference is the use of the AffineTransform here instead of the
// \doxygen{CenteredRigid2DTransform}. We will focus on the most
// relevant changes in the current code and skip the basic elements already
// explained in previous examples.
//
// \index{itk::AffineTransform}
//
// Software Guide : EndLatex 

#include "itkImageRegistrationMethod.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkGradientDescentTrustRegionOptimizer.h"
#include "itkImage.h"


#include "itkCenteredTransformInitializer.h"

//  Software Guide : BeginLatex
//  
//  Let's start by including the header file of the AffineTransform.
//
//  \index{itk::AffineTransform!header}
// 
//  Software Guide : EndLatex 

// Software Guide : BeginCodeSnippet
#include "itkAffineTransform.h"
// Software Guide : EndCodeSnippet


#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"


//
//  The following piece of code implements an observer
//  that will monitor the evolution of the registration process.
//
#include "itkCommand.h"
#include "CommandObserver.h"


int main( int argc, char *argv[] )
{
  if( argc < 5 )
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << "   rs|tr fixedImageFile  movingImageFile " << std::endl;
    std::cerr << "   outputImagefile  [differenceBeforeRegistration] " << std::endl;
    std::cerr << "   [differenceAfterRegistration] " << std::endl;
    std::cerr << "   [stepLength] [maxNumberOfIterations] "<< std::endl;
    return EXIT_FAILURE;
    }

  bool mode=true;
  {
     std::string argstr(argv[1]);
     if(argstr=="rs") mode=true;
     else if(argstr=="tr") mode=false;
     else {
        std::cerr << "Could not parse optimizer method, must be rs, or tr.  Got " << argstr << std::endl;
        return EXIT_FAILURE;
     }
  } 
 
  //  Software Guide : BeginLatex
  //
  //  We define then the types of the images to be registered.
  //
  //  Software Guide : EndLatex 

  // Software Guide : BeginCodeSnippet
  const    unsigned int    Dimension = 2;
  typedef  float           PixelType;

  typedef itk::Image< PixelType, Dimension >  FixedImageType;
  typedef itk::Image< PixelType, Dimension >  MovingImageType;
  // Software Guide : EndCodeSnippet


  //  Software Guide : BeginLatex
  //  
  //  The transform type is instantiated using the code below. The template
  //  parameters of this class are the representation type of the space
  //  coordinates and the space dimension.
  //
  //  \index{itk::AffineTransform!Instantiation}
  //
  //  Software Guide : EndLatex 

  // Software Guide : BeginCodeSnippet
  typedef itk::AffineTransform< 
                                  double, 
                                  Dimension  >     TransformType;
  // Software Guide : EndCodeSnippet


  typedef itk::RegularStepGradientDescentOptimizer       OptimizerType1;
  typedef itk::GradientDescentTrustRegionOptimizer       OptimizerType2;
  typedef itk::MeanSquaresImageToImageMetric< 
                                    FixedImageType, 
                                    MovingImageType >    MetricType;
  typedef itk:: LinearInterpolateImageFunction< 
                                    MovingImageType,
                                    double          >    InterpolatorType;
  typedef itk::ImageRegistrationMethod< 
                                    FixedImageType, 
                                    MovingImageType >    RegistrationType;

  MetricType::Pointer         metric        = MetricType::New();
  OptimizerType1::Pointer      optimizer1     = OptimizerType1::New();
  OptimizerType2::Pointer      optimizer2     = OptimizerType2::New();
  InterpolatorType::Pointer   interpolator  = InterpolatorType::New();
  RegistrationType::Pointer   registration  = RegistrationType::New();

  registration->SetMetric(        metric        );
  if(mode) {
     registration->SetOptimizer(     optimizer1     );
  } else {
     registration->SetOptimizer(     optimizer2     );
  }
  registration->SetInterpolator(  interpolator  );


  //  Software Guide : BeginLatex
  //
  //  The transform object is constructed below and passed to the registration
  //  method.
  //
  //  \index{itk::AffineTransform!New()}
  //  \index{itk::AffineTransform!Pointer}
  //  \index{itk::RegistrationMethod!SetTransform()}
  //
  //  Software Guide : EndLatex 

  // Software Guide : BeginCodeSnippet
  TransformType::Pointer  transform = TransformType::New();
  registration->SetTransform( transform );
  // Software Guide : EndCodeSnippet


  typedef itk::ImageFileReader< FixedImageType  > FixedImageReaderType;
  typedef itk::ImageFileReader< MovingImageType > MovingImageReaderType;
  FixedImageReaderType::Pointer  fixedImageReader  = FixedImageReaderType::New();
  MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();
  fixedImageReader->SetFileName(  argv[2] );
  movingImageReader->SetFileName( argv[3] );


  registration->SetFixedImage(    fixedImageReader->GetOutput()    );
  registration->SetMovingImage(   movingImageReader->GetOutput()   );
  fixedImageReader->Update();

  registration->SetFixedImageRegion( 
     fixedImageReader->GetOutput()->GetBufferedRegion() );


  //  Software Guide : BeginLatex
  //  
  //  In this example, we again use the
  //  \doxygen{CenteredTransformInitializer} helper class in order to compute
  //  a reasonable value for the initial center of rotation and the
  //  translation. The initializer is set to use the center of mass of each
  //  image as the initial correspondence correction.
  //
  //  Software Guide : EndLatex 

  // Software Guide : BeginCodeSnippet
  typedef itk::CenteredTransformInitializer< 
                                    TransformType, 
                                    FixedImageType, 
                                    MovingImageType >  TransformInitializerType;
  TransformInitializerType::Pointer initializer = TransformInitializerType::New();
  initializer->SetTransform(   transform );
  initializer->SetFixedImage(  fixedImageReader->GetOutput() );
  initializer->SetMovingImage( movingImageReader->GetOutput() );
  initializer->MomentsOn();
  initializer->InitializeTransform();
  // Software Guide : EndCodeSnippet


  //  Software Guide : BeginLatex
  //  
  //  Now we pass the parameters of the current transform as the initial
  //  parameters to be used when the registration process starts.
  //
  //  Software Guide : EndLatex 

  // Software Guide : BeginCodeSnippet
  registration->SetInitialTransformParameters( 
                                 transform->GetParameters() );
  // Software Guide : EndCodeSnippet


  //  Software Guide : BeginLatex
  //  
  //  Keeping in mind that the scale of units in scaling, rotation and
  //  translation are quite different, we take advantage of the scaling
  //  functionality provided by the optimizers. We know that the first $N
  //  \times N$ elements of the parameters array correspond to the rotation
  //  matrix factor, and the last $N$ are the components of the translation to
  //  be applied after multiplication with the matrix is performed.
  //
  //  Software Guide : EndLatex 


  double translationScale = 1.0 / 1000.0;
  if( argc > 9 )
    {
    translationScale = atof( argv[9] );
    }


  // Software Guide : BeginCodeSnippet
  typedef OptimizerType1::ScalesType       OptimizerScalesType;
  OptimizerScalesType optimizerScales( transform->GetNumberOfParameters() );

  optimizerScales[0] =  1.0;
  optimizerScales[1] =  1.0;
  optimizerScales[2] =  1.0;
  optimizerScales[3] =  1.0;
  optimizerScales[4] =  translationScale;
  optimizerScales[5] =  translationScale;

  optimizer1->SetScales( optimizerScales );
  optimizer2->SetScales( optimizerScales );
  // Software Guide : EndCodeSnippet


  //  Software Guide : BeginLatex
  //  
  //  We also set the usual parameters of the optimization method. In this
  //  case we are using an
  //  \doxygen{RegularStepGradientDescentOptimizer}. Below, we define the
  //  optimization parameters like initial step length, minimal step length
  //  and number of iterations. These last two act as stopping criteria for
  //  the optimization.
  //
  //  Software Guide : EndLatex 

  double steplength = 0.1;

  if( argc > 7 )
    {
    steplength = atof( argv[7] );
    }


  unsigned int maxNumberOfIterations = 300;

  if( argc > 8 )
    {
    maxNumberOfIterations = atoi( argv[8] );
    }


  // Software Guide : BeginCodeSnippet
  optimizer1->SetMaximumStepLength( steplength ); 
  optimizer1->SetMinimumStepLength( 0.0001 );
  optimizer1->SetNumberOfIterations( maxNumberOfIterations );
  optimizer2->SetInitialStepLength( steplength ); 
  optimizer2->SetMinimumStepLength( 0.0001 );
  optimizer2->SetNumberOfIterations( maxNumberOfIterations );
  // Software Guide : EndCodeSnippet


  //  Software Guide : BeginLatex
  //
  //  We also set the optimizer to do minimization by calling the
  //  \code{MinimizeOn()} method.
  //
  //  \index{itk::Regular\-Step\-Gradient\-Descent\-Optimizer!MinimizeOn()}
  //
  //  Software Guide : EndLatex 

  // Software Guide : BeginCodeSnippet
  optimizer1->MinimizeOn();
  optimizer2->MinimizeOn();
  // Software Guide : EndCodeSnippet


  // Create the Command observer and register it with the optimizer.
  //
  CommandIterationUpdate<OptimizerType1>::Pointer observer1 = CommandIterationUpdate<OptimizerType1>::New();
  optimizer1->AddObserver( itk::IterationEvent(), observer1 );
  CommandIterationUpdate<OptimizerType2>::Pointer observer2 = CommandIterationUpdate<OptimizerType2>::New();
  optimizer2->AddObserver( itk::IterationEvent(), observer2 );



  //  Software Guide : BeginLatex
  //
  //  Finally we trigger the execution of the registration method by calling
  //  the \code{Update()} method. The call is placed in a \code{try/catch}
  //  block in case any exceptions are thrown.
  //
  //  Software Guide : EndLatex 

  // Software Guide : BeginCodeSnippet
  try 
    { 
    registration->StartRegistration(); 
    std::cout << "Optimizer stop condition: "
              << registration->GetOptimizer()->GetStopConditionDescription()
              << std::endl;
    } 
  catch( itk::ExceptionObject & err ) 
    { 
    std::cerr << "ExceptionObject caught !" << std::endl; 
    std::cerr << err << std::endl; 
    return EXIT_FAILURE;
    } 
  // Software Guide : EndCodeSnippet


  //  Software Guide : BeginLatex
  //
  //  Once the optimization converges, we recover the parameters from the
  //  registration method. This is done with the
  //  \code{GetLastTransformParameters()} method. We can also recover the
  //  final value of the metric with the \code{GetValue()} method and the
  //  final number of iterations with the \code{GetCurrentIteration()}
  //  method.
  //
  //  \index{itk::RegistrationMethod!GetValue()}
  //  \index{itk::RegistrationMethod!GetCurrentIteration()}
  //  \index{itk::RegistrationMethod!GetLastTransformParameters()}
  //
  //  Software Guide : EndLatex 

  // Software Guide : BeginCodeSnippet
  OptimizerType1::ParametersType finalParameters = 
                    registration->GetLastTransformParameters();

  const double finalRotationCenterX = transform->GetCenter()[0];
  const double finalRotationCenterY = transform->GetCenter()[1];
  const double finalTranslationX    = finalParameters[4];
  const double finalTranslationY    = finalParameters[5];

  const unsigned int numberOfIterations = mode ? optimizer1->GetCurrentIteration() : optimizer2->GetCurrentIteration();
  const double bestValue = mode ? optimizer1->GetValue() : optimizer2->GetValue();

  // Software Guide : EndCodeSnippet


  // Print out results
  //
  std::cout << "Result = " << std::endl;
  std::cout << " Center X      = " << finalRotationCenterX  << std::endl;
  std::cout << " Center Y      = " << finalRotationCenterY  << std::endl;
  std::cout << " Translation X = " << finalTranslationX  << std::endl;
  std::cout << " Translation Y = " << finalTranslationY  << std::endl;
  std::cout << " Iterations    = " << numberOfIterations << std::endl;
  std::cout << " Metric value  = " << bestValue          << std::endl;
  
  //Compute the rotation angle and scaling from SVD of the matrix
  // \todo Find a way to figure out if the scales are along X or along Y.
  // VNL returns the eigenvalues ordered from largest to smallest.
  
  vnl_matrix<double> p(2, 2);
  p[0][0] = (double) finalParameters[0];
  p[0][1] = (double) finalParameters[1];
  p[1][0] = (double) finalParameters[2];
  p[1][1] = (double) finalParameters[3];
  vnl_svd<double> svd(p);
  vnl_matrix<double> r(2, 2);
  r = svd.U() * vnl_transpose(svd.V());
  double angle = vcl_asin(r[1][0]);
  
  const double angleInDegrees = angle * 180.0 / vnl_math::pi;

  std::cout << " Scale 1         = " << svd.W(0)        << std::endl;
  std::cout << " Scale 2         = " << svd.W(1)        << std::endl;
  std::cout << " Angle (degrees) = " << angleInDegrees  << std::endl;
  

  //  Software Guide : BeginLatex
  //  
  //  Let's execute this example over two of the images provided in
  //  \code{Examples/Data}:
  //  
  //  \begin{itemize}
  //  \item \code{BrainProtonDensitySliceBorder20.png} 
  //  \item \code{BrainProtonDensitySliceR10X13Y17.png}
  //  \end{itemize}
  //
  //  The second image is the result of intentionally rotating the first
  //  image by $10$ degrees and then translating by $(-13,-17)$.  Both images
  //  have unit-spacing and are shown in Figure
  //  \ref{fig:FixedMovingImageRegistration9}. We execute the code using the
  //  following parameters: step length=1.0, translation scale= 0.0001 and
  //  maximum number of iterations = 300. With these images and parameters
  //  the registration takes $98$ iterations and produces
  //
  //  \begin{center}
  //  \begin{verbatim}
  //   96 58.09 [0.986481, -0.169104, 0.166411, 0.986174, 12.461, 16.0754]
  //  \end{verbatim}
  //  \end{center}
  //
  //  These results are interpreted as
  //
  //  \begin{itemize}
  //  \item Iterations   = 98
  //  \item Final Metric = 58.09
  //  \item Center       = $( 111.204,   131.6   )$ millimeters
  //  \item Translation  = $(   12.461,  16.0754 )$ millimeters
  //  \item Affine scales = $(1.00185, .999137)$
  //  \end{itemize}
  //  
  //  The second component of the matrix values is usually associated with
  //  $\sin{\theta}$. We obtain the rotation through SVD of the affine
  //  matrix. The value is $9.6526$ degrees, which is approximately the 
  //  intentional misalignment of $10.0$ degrees.
  //
  // \begin{figure}
  // \center
  // \includegraphics[width=0.44\textwidth]{BrainProtonDensitySliceBorder20.eps}
  // \includegraphics[width=0.44\textwidth]{BrainProtonDensitySliceR10X13Y17.eps}
  // \itkcaption[AffineTransform registration]{Fixed and moving images
  // provided as input to the registration method using the AffineTransform.}
  // \label{fig:FixedMovingImageRegistration9}
  // \end{figure}
  //
  //
  // \begin{figure}
  // \center
  // \includegraphics[width=0.32\textwidth]{ImageRegistration9Output.eps}
  // \includegraphics[width=0.32\textwidth]{ImageRegistration9DifferenceBefore.eps}
  // \includegraphics[width=0.32\textwidth]{ImageRegistration9DifferenceAfter.eps} 
  // \itkcaption[AffineTransform output images]{The resampled moving image
  // (left), and the difference between the fixed and moving images before (center) 
  // and after (right) registration with the 
  // AffineTransform transform.}
  // \label{fig:ImageRegistration9Outputs}
  // \end{figure}
  //
  // Figure \ref{fig:ImageRegistration9Outputs} shows the output of the
  // registration. The right most image of this figure shows the squared
  // magnitude difference between the fixed image and the resampled
  // moving image. 
  //
  // \begin{figure}
  // \center
  // \includegraphics[height=0.32\textwidth]{ImageRegistration9TraceMetric.eps}
  // \includegraphics[height=0.32\textwidth]{ImageRegistration9TraceAngle.eps}
  // \includegraphics[height=0.32\textwidth]{ImageRegistration9TraceTranslations.eps} 
  // \itkcaption[AffineTransform output plots]{Metric values,
  // rotation angle and translations during the registration using the 
  // AffineTransform transform.}
  // \label{fig:ImageRegistration9Plots}
  // \end{figure}
  //
  //  Figure \ref{fig:ImageRegistration9Plots} shows the plots of the main
  //  output parameters of the registration process. The metric values at every
  //  iteration are shown on the top plot. The angle values are shown on the bottom left plot,
  //  while the translation components of the registration are presented
  //  on the bottom right plot. Note that the final total offset of the transform
  //  is to be computed as a combination of the shift due rotation plus the
  //  explicit translation set on the transform.
  //
  //  Software Guide : EndLatex 


  //  The following code is used to dump output images to files.
  //  They illustrate the final results of the registration.
  //  We will resample the moving image and write out the difference image
  //  before and after registration. We will also rescale the intensities of the
  //  difference images, so that they look better!
  typedef itk::ResampleImageFilter< 
                            MovingImageType, 
                            FixedImageType >    ResampleFilterType;

  TransformType::Pointer finalTransform = TransformType::New();

  finalTransform->SetParameters( finalParameters );
  finalTransform->SetFixedParameters( transform->GetFixedParameters() );

  ResampleFilterType::Pointer resampler = ResampleFilterType::New();

  resampler->SetTransform( finalTransform );
  resampler->SetInput( movingImageReader->GetOutput() );

  FixedImageType::Pointer fixedImage = fixedImageReader->GetOutput();

  resampler->SetSize(    fixedImage->GetLargestPossibleRegion().GetSize() );
  resampler->SetOutputOrigin(  fixedImage->GetOrigin() );
  resampler->SetOutputSpacing( fixedImage->GetSpacing() );
  resampler->SetOutputDirection( fixedImage->GetDirection() );
  resampler->SetDefaultPixelValue( 100 );
  
  typedef  unsigned char  OutputPixelType;

  typedef itk::Image< OutputPixelType, Dimension > OutputImageType;
  
  typedef itk::CastImageFilter< 
                        FixedImageType,
                        OutputImageType > CastFilterType;
                    
  typedef itk::ImageFileWriter< OutputImageType >  WriterType;


  WriterType::Pointer      writer =  WriterType::New();
  CastFilterType::Pointer  caster =  CastFilterType::New();


  writer->SetFileName( argv[4] );


  caster->SetInput( resampler->GetOutput() );
  writer->SetInput( caster->GetOutput()   );
  writer->Update();


  typedef itk::SubtractImageFilter< 
                                  FixedImageType, 
                                  FixedImageType, 
                                  FixedImageType > DifferenceFilterType;

  DifferenceFilterType::Pointer difference = DifferenceFilterType::New();

  difference->SetInput1( fixedImageReader->GetOutput() );
  difference->SetInput2( resampler->GetOutput() );

  WriterType::Pointer writer2 = WriterType::New();
  
  typedef itk::RescaleIntensityImageFilter< 
                                  FixedImageType, 
                                  OutputImageType >   RescalerType;

  RescalerType::Pointer intensityRescaler = RescalerType::New();

  intensityRescaler->SetInput( difference->GetOutput() );
  intensityRescaler->SetOutputMinimum(   0 );
  intensityRescaler->SetOutputMaximum( 255 );
  
  writer2->SetInput( intensityRescaler->GetOutput() );  
  resampler->SetDefaultPixelValue( 1 );
  
  // Compute the difference image between the 
  // fixed and resampled moving image.
  if( argc > 6 )
    {
    writer2->SetFileName( argv[6] );
    writer2->Update();
    }


  typedef itk::IdentityTransform< double, Dimension > IdentityTransformType;
  IdentityTransformType::Pointer identity = IdentityTransformType::New();

  // Compute the difference image between the 
  // fixed and moving image before registration.
  if( argc > 5 )
    {
    resampler->SetTransform( identity );
    writer2->SetFileName( argv[5] );
    writer2->Update();
    }

  return EXIT_SUCCESS;
}
