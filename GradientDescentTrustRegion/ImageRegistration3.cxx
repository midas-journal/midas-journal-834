/*=========================================================================

  This has been modified from the original included with ITK 3.20 to show
  the performance of the GradientDescentTrustRegionOptimizer.

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: ImageRegistration3.cxx,v $
  Language:  C++
  Date:      $Date: 2009-06-24 12:09:09 $
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

// Software Guide : BeginLatex
//
// Given the numerous parameters involved in tuning a registration method for
// a particular application, it is not uncommon for a registration process to
// run for several minutes and still produce a useless result.  To avoid
// this situation it is quite helpful to track the evolution of the
// registration as it progresses. The following section illustrates the
// mechanisms provided in ITK for monitoring the activity of the
// ImageRegistrationMethod class.
//
// Insight implements the \emph{Observer/Command} design pattern
// \cite{Gamma1995}. (See Section~\ref{sec:EventHandling} for an overview.) 
// The classes involved in this implementation are the \doxygen{Object},
// \doxygen{Command} and \doxygen{EventObject} classes. The Object
// is the base class of most ITK objects. This class maintains a linked
// list of pointers to event observers. The role of observers is played by
// the Command class.  Observers register themselves with an
// Object, declaring that they are interested in receiving
// notification when a particular event happens. A set of events is
// represented by the hierarchy of the Event class. Typical events
// are \code{Start}, \code{End}, \code{Progress} and \code{Iteration}.
//
// Registration is controlled by an \doxygen{Optimizer}, which generally
// executes an iterative process. Most Optimizer classes invoke an
// \doxygen{IterationEvent} at the end of each iteration. When an event is
// invoked by an object, this object goes through its list of registered
// observers (Commands) and checks whether any one of them has expressed
// interest in the current event type. Whenever such an observer is found,
// its corresponding \code{Execute()} method is invoked.  In this context,
// \code{Execute()} methods should be considered \emph{callbacks}.  As such,
// some of the common sense rules of callbacks should be respected.  For
// example, \code{Execute()} methods should not perform heavy computational
// tasks.  They are expected to execute rapidly, for example, printing out a
// message or updating a value in a GUI.
//
// \index{itk::ImageRegistrationMethod!Monitoring}
//
//
// Software Guide : EndLatex 


#include "itkImageRegistrationMethod.h"
#include "itkTranslationTransform.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkGradientDescentTrustRegionOptimizer.h"
#include "itkImage.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"


//  Software Guide : BeginLatex
//
//  The following code illustrates a simple way of creating a 
//  Observer/Command to monitor a registration process. This new
//  class derives from the Command class and provides a specific
//  implementation of the \code{Execute()} method.  First, the header file of
//  the Command class must be included.
//
//  Software Guide : EndLatex 

// Software Guide : BeginCodeSnippet
#include "itkCommand.h"
// Software Guide : EndCodeSnippet


//  Software Guide : BeginLatex
//
//  Our custom command class is called \code{CommandIterationUpdate}. It
//  derives from the Command class and declares for convenience the
//  types \code{Self} and \code{Superclass}. This facilitate the use of
//  standard macros later in the class implementation.
//
//  Software Guide : EndLatex 

// Software Guide : BeginCodeSnippet
#include "CommandObserver.h"


int main( int argc, char *argv[] )
{
  if( argc < 4 )
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " rs|tr fixedImageFile  movingImageFile ";
    std::cerr << "outputImagefile [differenceImage]" << std::endl;
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

  typedef itk::MeanSquaresImageToImageMetric< 
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
  
  registration->SetInitialTransformParameters( initialParameters );

  optimizer1->SetMaximumStepLength( 4.00 );  
  optimizer1->SetMinimumStepLength( 0.01 );
  optimizer1->SetNumberOfIterations( 200 );

  optimizer1->MaximizeOff();

  optimizer2->SetInitialStepLength( 4.00 );  
  optimizer2->SetMinimumStepLength( 0.01 );
  optimizer2->SetNumberOfIterations( 200 );

  optimizer2->MaximizeOff();


  //  Software Guide : BeginLatex
  //
  //  Once all the registration components are in place we can create one
  //  instance of our observer. This is done with the standard \code{New()}
  //  method and assigned to a SmartPointer.
  //
  //  Software Guide : EndLatex 

  // Software Guide : BeginCodeSnippet
  CommandIterationUpdate<OptimizerType1>::Pointer observer1 = CommandIterationUpdate<OptimizerType1>::New();
  CommandIterationUpdate<OptimizerType2>::Pointer observer2 = CommandIterationUpdate<OptimizerType2>::New();
  // Software Guide : EndCodeSnippet


  //  Software Guide : BeginLatex
  //  
  // \begin{figure}
  // \center
  // \includegraphics[width=\textwidth]{ImageRegistration3Observer.eps}
  // \itkcaption[Command/Observer and the Registration Framework]{Interaction
  // between the Command/Observer and the Registration Method.}
  // \label{fig:ImageRegistration3Observer}
  // \end{figure}
  //
  //
  //  Software Guide : EndLatex 


  //  Software Guide : BeginLatex
  //
  //  The newly created command is registered as observer on the
  //  optimizer, using the \code{AddObserver()} method. Note
  //  that the event type is provided as the first argument to this
  //  method. In order for the RTTI mechanism to work correctly, a newly
  //  created event of the desired type must be passed as the first
  //  argument. The second argument is simply the smart pointer to the
  //  optimizer. Figure \ref{fig:ImageRegistration3Observer} illustrates the
  //  interaction between the Command/Observer class and the registration
  //  method.
  //
  //  Software Guide : EndLatex 

  // Software Guide : BeginCodeSnippet
  optimizer1->AddObserver( itk::IterationEvent(), observer1 );
  optimizer2->AddObserver( itk::IterationEvent(), observer2 );
  // Software Guide : EndCodeSnippet


  //  Software Guide : BeginLatex
  //
  //  At this point, we are ready to execute the registration. The
  //  typical call to \code{StartRegistration()} will do it. Note again the
  //  use of the \code{try/catch} block around the \code{StartRegistration()}
  //  method in case an exception is thrown.
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
    std::cout << "ExceptionObject caught !" << std::endl; 
    std::cout << err << std::endl; 
    return EXIT_FAILURE;
    } 
  // Software Guide : EndCodeSnippet


  //  Software Guide : BeginLatex
  //
  //  The registration process is applied to the following images in \code{Examples/Data}:
  //  
  //  \begin{itemize}
  //  \item \code{BrainProtonDensitySliceBorder20.png} 
  //  \item \code{BrainProtonDensitySliceShifted13x17y.png}
  //  \end{itemize}
  //
  //  It produces the following output.
  //
  //  \begin{verbatim}
  //   0 = 4499.45 : [2.9287, 2.72447]
  //   1 = 3860.84 : [5.62751, 5.67683]
  //   2 = 3450.68 : [8.85516, 8.03952]
  //   3 = 3152.07 : [11.7997, 10.7469]
  //   4 = 2189.97 : [13.3628, 14.4288]
  //   5 = 1047.21 : [11.292, 17.851]
  //   6 = 900.189 : [13.1602, 17.1372]
  //   7 = 19.6301 : [12.3268, 16.5846]
  //   8 = 237.317 : [12.7824, 16.7906]
  //   9 = 38.1331 : [13.1833, 17.0894]
  //   10 = 18.9201 : [12.949, 17.002]
  //   11 = 1.15456 : [13.074, 16.9979]
  //   12 = 2.42488 : [13.0115, 16.9994]
  //   13 = 0.0590549 : [12.949, 17.002]
  //   14 = 1.15451 : [12.9803, 17.001]
  //   15 = 0.173731 : [13.0115, 16.9997]
  //   16 = 0.0586584 : [12.9959, 17.0001]
  //  \end{verbatim}
  //  You can verify from the code in the \code{Execute()} method that the first
  //  column is the iteration number, the second column is the metric value and
  //  the third and fourth columns are the parameters of the transform, which
  //  is a $2D$ translation transform in this case. By tracking these values as
  //  the registration progresses, you will be able to determine whether the
  //  optimizer is advancing in the right direction and whether the step-length
  //  is reasonable or not.  That will allow you to interrupt the registration
  //  process and fine-tune parameters without having to wait until the
  //  optimizer stops by itself.
  //
  //  Software Guide : EndLatex 


  ParametersType finalParameters = registration->GetLastTransformParameters();
  
  const double TranslationAlongX = finalParameters[0];
  const double TranslationAlongY = finalParameters[1];
  
  const unsigned int numberOfIterations = (mode) ? optimizer1->GetCurrentIteration() : optimizer2->GetCurrentIteration();
  
  const double bestValue = mode ? optimizer1->GetValue() : optimizer2->GetValue();

  std::cout << "Registration done !" << std::endl;
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
