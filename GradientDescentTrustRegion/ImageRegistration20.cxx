/*=========================================================================
 *
  This has been modified from the original included with ITK 3.20 to show
  the performance of the GradientDescentTrustRegionOptimizer.

 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

//  Software Guide : BeginCommandLineArgs
//    INPUTS: {brainweb1e1a10f20.mha}
//    INPUTS: {brainweb1e1a10f20Rot10Tx15.mha}
//    ImageRegistration20Output.mhd
//  Software Guide : EndCommandLineArgs

// Software Guide : BeginLatex
//
// This example illustrates the use of the \doxygen{AffineTransform}
// for performing registration in $3D$.
//
// \index{itk::AffineTransform}
//
// Software Guide : EndLatex

#include "itkImageRegistrationMethod.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkCenteredTransformInitializer.h"
#include "itkGradientDescentTrustRegionOptimizer.h"

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
  const    unsigned int    Dimension = 3;
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
  optimizerScales[4] =  1.0;
  optimizerScales[5] =  1.0;
  optimizerScales[6] =  1.0;
  optimizerScales[7] =  1.0;
  optimizerScales[8] =  1.0;
  optimizerScales[9]  =  translationScale;
  optimizerScales[10] =  translationScale;
  optimizerScales[11] =  translationScale;

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

  const unsigned int numberOfIterations = mode ? optimizer1->GetCurrentIteration() : optimizer2->GetCurrentIteration();
  const double bestValue = mode ? optimizer1->GetValue(): optimizer2->GetValue();
  // Software Guide : EndCodeSnippet


  // Print out results
  //
  std::cout << "Result = " << std::endl;
  std::cout << " Iterations    = " << numberOfIterations << std::endl;
  std::cout << " Metric value  = " << bestValue          << std::endl;


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
