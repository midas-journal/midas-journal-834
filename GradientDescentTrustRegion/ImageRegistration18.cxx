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

// Software Guide : BeginLatex
//
// This example illustrates how to use the
// \doxygen{GradientDifferenceImageToImageMetric}.
//
// This metric is particularly useful for registration scenarios where fitting
// the edges of both images is the most relevant criteria for registration
// success.
//
// \index{itk::ImageRegistrationMethod!Monitoring}
//
//
// Software Guide : EndLatex


#include "itkImageRegistrationMethod.h"
#include "itkTranslationTransform.h"
#include "itkGradientDifferenceImageToImageMetric.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkGradientDescentTrustRegionOptimizer.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkCommand.h"
#include "CommandObserver.h"


int main( int argc, char *argv[] )
{
  if( argc < 4 )
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " rs|tr fixedImageFile  movingImageFile ";
    std::cerr << "outputImagefile" << std::endl;
    std::cerr << "[initialTx] [initialTy]" << std::endl;
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


  const    unsigned int    Dimension = 2;
  typedef  unsigned short  PixelType;

  typedef itk::Image< PixelType, Dimension >  FixedImageType;
  typedef itk::Image< PixelType, Dimension >  MovingImageType;

  typedef itk::TranslationTransform< double, Dimension > TransformType;

  typedef itk::RegularStepGradientDescentOptimizer       OptimizerType1;
  typedef itk::GradientDescentTrustRegionOptimizer       OptimizerType2;
 
  typedef itk::LinearInterpolateImageFunction<
                                    MovingImageType,
                                    double             > InterpolatorType;

  typedef itk::ImageRegistrationMethod<
                                    FixedImageType,
                                    MovingImageType   >  RegistrationType;

  typedef itk::GradientDifferenceImageToImageMetric<
                                      FixedImageType,
                                      MovingImageType >  MetricType;

  TransformType::Pointer      transform     = TransformType::New();
  OptimizerType1::Pointer      optimizer1     = OptimizerType1::New();
  OptimizerType2::Pointer      optimizer2     = OptimizerType2::New();
   InterpolatorType::Pointer   interpolator  = InterpolatorType::New();
  RegistrationType::Pointer   registration  = RegistrationType::New();

  if(mode) {
     registration->SetOptimizer(     optimizer1     );
  } else {
     registration->SetOptimizer(     optimizer2     );
  }
   registration->SetTransform(     transform     );
  registration->SetInterpolator(  interpolator  );

  MetricType::Pointer         metric        = MetricType::New();

  metric->SetDerivativeDelta( 0.5 );

  registration->SetMetric( metric  );

  typedef itk::ImageFileReader< FixedImageType  > FixedImageReaderType;
  typedef itk::ImageFileReader< MovingImageType > MovingImageReaderType;

  FixedImageReaderType::Pointer  fixedImageReader  = FixedImageReaderType::New();
  MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();

  fixedImageReader->SetFileName(  argv[2] );
  movingImageReader->SetFileName( argv[3] );

  registration->SetFixedImage(    fixedImageReader->GetOutput()    );
  registration->SetMovingImage(   movingImageReader->GetOutput()   );

  fixedImageReader->Update(); // This is needed to make the BufferedRegion below valid.

  registration->SetFixedImageRegion(
       fixedImageReader->GetOutput()->GetBufferedRegion() );

  typedef RegistrationType::ParametersType ParametersType;
  ParametersType initialParameters( transform->GetNumberOfParameters() );

  initialParameters[0] = 0.0;  // Initial offset in mm along X
  initialParameters[1] = 0.0;  // Initial offset in mm along Y

  if( argc > 5 )
    {
    initialParameters[0] = atof( argv[5] );
    }

  if( argc > 6 )
    {
    initialParameters[1] = atof( argv[6] );
    }

  std::cout << "Initial parameters = " << initialParameters << std::endl;

  registration->SetInitialTransformParameters( initialParameters );

  optimizer1->SetMaximumStepLength( 4.00 );
  optimizer1->SetMinimumStepLength( 0.01 );
  optimizer1->SetNumberOfIterations( 200 );
  optimizer1->SetGradientMagnitudeTolerance( 1e-40 );

  optimizer1->MaximizeOn();

  optimizer2->SetInitialStepLength( 4.00 );
  optimizer2->SetMinimumStepLength( 0.01 );
  optimizer2->SetNumberOfIterations( 200 );
  optimizer2->SetGradientMagnitudeTolerance( 1e-40 );

  optimizer2->MaximizeOn();


  CommandIterationUpdate<OptimizerType1>::Pointer observer1 = CommandIterationUpdate<OptimizerType1>::New();
  optimizer1->AddObserver( itk::IterationEvent(), observer1 );
  CommandIterationUpdate<OptimizerType2>::Pointer observer2 = CommandIterationUpdate<OptimizerType2>::New();
  optimizer2->AddObserver( itk::IterationEvent(), observer2 );


  try
    {
    registration->StartRegistration();
    std::cout << "Optimizer stop condition: "
              << registration->GetOptimizer()->GetStopConditionDescription()
              << std::endl;
    }
  catch( itk::ExceptionObject & err )
    {
    std::cout << "ExceptionObject caught !" << std::endl;
    std::cout << err << std::endl;
    return EXIT_FAILURE;
    }


  ParametersType finalParameters = registration->GetLastTransformParameters();

  const double TranslationAlongX = finalParameters[0];
  const double TranslationAlongY = finalParameters[1];

  const unsigned int numberOfIterations = mode ? optimizer1->GetCurrentIteration() : optimizer2->GetCurrentIteration();
  const double bestValue = mode ? optimizer1->GetValue(): optimizer2->GetValue();

  std::cout << "Registration done !" << std::endl;
  std::cout << "Optimizer stop condition = "
              << registration->GetOptimizer()->GetStopConditionDescription()
              << std::endl;
  std::cout << "Number of iterations = " << numberOfIterations << std::endl;
  std::cout << "Translation along X  = " << TranslationAlongX << std::endl;
  std::cout << "Translation along Y  = " << TranslationAlongY << std::endl;
  std::cout << "Optimal metric value = " << bestValue << std::endl;


  // Prepare the resampling filter in order to map the moving image.
  //
  typedef itk::ResampleImageFilter<
                            MovingImageType,
                            FixedImageType >    ResampleFilterType;

  TransformType::Pointer finalTransform = TransformType::New();

  finalTransform->SetParameters( finalParameters );
  finalTransform->SetFixedParameters( transform->GetFixedParameters() );

  ResampleFilterType::Pointer resample = ResampleFilterType::New();

  resample->SetTransform( finalTransform );
  resample->SetInput( movingImageReader->GetOutput() );

  FixedImageType::Pointer fixedImage = fixedImageReader->GetOutput();

  resample->SetSize(    fixedImage->GetLargestPossibleRegion().GetSize() );
  resample->SetOutputOrigin(  fixedImage->GetOrigin() );
  resample->SetOutputSpacing( fixedImage->GetSpacing() );
  resample->SetOutputDirection( fixedImage->GetDirection() );
  resample->SetDefaultPixelValue( 100 );


  // Prepare a writer and caster filters to send the resampled moving image to
  // a file
  //
  typedef  unsigned char  OutputPixelType;

  typedef itk::Image< OutputPixelType, Dimension > OutputImageType;

  typedef itk::CastImageFilter<
                        FixedImageType,
                        OutputImageType > CastFilterType;

  typedef itk::ImageFileWriter< OutputImageType >  WriterType;

  WriterType::Pointer      writer =  WriterType::New();
  CastFilterType::Pointer  caster =  CastFilterType::New();


  writer->SetFileName( argv[4] );

  caster->SetInput( resample->GetOutput() );
  writer->SetInput( caster->GetOutput()   );
  writer->Update();


  return EXIT_SUCCESS;
}
