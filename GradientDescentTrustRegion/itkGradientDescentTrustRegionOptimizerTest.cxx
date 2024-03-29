/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkRegularStepGradientDescentOptimizerTest.cxx,v $
  Language:  C++
  Date:      $Date: 2007-09-10 16:32:01 $
  Version:   $Revision: 1.23 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include "itkGradientDescentTrustRegionOptimizer.h"
#include <vnl/vnl_math.h>

/** 
 *  The objectif function is the quadratic form:
 *
 *  1/2 x^T A x - b^T x
 *
 *  Where A is a matrix and b is a vector
 *  The system in this example is:
 *
 *     | 3  2 ||x|   | 2|   |0|
 *     | 2  6 ||y| + |-8| = |0|
 *
 *
 *   the solution is the vector | 2 -2 |
 *
 */ 
class RSGCostFunction : public itk::SingleValuedCostFunction 
{
public:

  typedef RSGCostFunction                     Self;
  typedef itk::SingleValuedCostFunction      Superclass;
  typedef itk::SmartPointer<Self>            Pointer;
  typedef itk::SmartPointer<const Self>      ConstPointer;
  itkNewMacro( Self );

  enum { SpaceDimension=2 };
  
  typedef Superclass::ParametersType      ParametersType;
  typedef Superclass::DerivativeType      DerivativeType;
  typedef Superclass::MeasureType         MeasureType ;


  RSGCostFunction() 
  {
     m_Factor=1.0;
  }


  MeasureType  GetValue( const ParametersType & parameters ) const
  { 
    
    double x = parameters[0];
    double y = parameters[1];

    std::cout << "GetValue( " ;
    std::cout << x << " ";
    std::cout << y << ") = ";

    MeasureType measure = (0.5*(3*x*x+4*x*y+6*y*y) - 2*x + 8*y)*m_Factor;

    std::cout << measure << std::endl; 
    return measure;

  }

  void GetDerivative( const ParametersType & parameters,
                            DerivativeType  & derivative ) const
  {

    double x = parameters[0];
    double y = parameters[1];

    std::cout << "GetDerivative( " ;
    std::cout << x << " ";
    std::cout << y << ") = ";

    derivative = DerivativeType( SpaceDimension ); 
    derivative[0] = (3 * x + 2 * y -2)*m_Factor;
    derivative[1] = (2 * x + 6 * y +8)*m_Factor;
    std::cout << derivative << std::endl;

  }

 
  unsigned int GetNumberOfParameters(void) const
    {
    return SpaceDimension;
    }

  itkSetMacro(Factor,double);

private:
  double m_Factor;


};



int itkGradientDescentTrustRegionOptimizerTest(int, char* [] ) 
{
  std::cout << "GradientDescentTrustRegionOptimizer Test ";
  std::cout << std::endl << std::endl;

  typedef  itk::GradientDescentTrustRegionOptimizer  OptimizerType;

  typedef  OptimizerType::ScalesType            ScalesType;
  
  
  // Declaration of a itkOptimizer
  OptimizerType::Pointer  itkOptimizer = OptimizerType::New();


  // Declaration of the CostFunction 
  RSGCostFunction::Pointer costFunction = RSGCostFunction::New();


  itkOptimizer->SetCostFunction( costFunction.GetPointer() );

  
  typedef RSGCostFunction::ParametersType    ParametersType;


  const unsigned int spaceDimension = 
                      costFunction->GetNumberOfParameters();

  // We start not so far from  | 2 -2 |
  ParametersType  initialPosition( spaceDimension );
  initialPosition[0] =  100;
  initialPosition[1] = -100;
  
  ScalesType    parametersScale( spaceDimension );
  parametersScale[0] = 1.0;
  parametersScale[1] = 1.0;

  itkOptimizer->MinimizeOn();
  itkOptimizer->SetScales( parametersScale );
  itkOptimizer->SetGradientMagnitudeTolerance( 1e-6 );
  itkOptimizer->SetMaximumStepLength( 30.0 );
  itkOptimizer->SetMinimumStepLength( 1e-6 );
  itkOptimizer->SetNumberOfIterations( 900 );

  itkOptimizer->SetInitialPosition( initialPosition );

  try 
    {
    itkOptimizer->StartOptimization();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cout << "Exception thrown ! " << std::endl;
    std::cout << "An error ocurred during Optimization" << std::endl;
    std::cout << "Location    = " << e.GetLocation()    << std::endl;
    std::cout << "Description = " << e.GetDescription() << std::endl;
    return EXIT_FAILURE;
    }


  ParametersType finalPosition = itkOptimizer->GetCurrentPosition();
  std::cout << "Solution        = (";
  std::cout << finalPosition[0] << "," ;
  std::cout << finalPosition[1] << ")" << std::endl;  
  std::cout << "Number of iterations = " << itkOptimizer->GetCurrentIteration() << std::endl;
  //
  // check results to see if it is within range
  //
  bool pass = true;
  double trueParameters[2] = { 2, -2 };
  for( unsigned int j = 0; j < 2; j++ )
    {
    if( vnl_math_abs( finalPosition[j] - trueParameters[j] ) > 0.01 )
      {
      pass = false;
      }
    }

  if( !pass )
    {
    std::cout << "Test failed." << std::endl;
    return EXIT_FAILURE;
    }


  // Run now with a different config
  itkOptimizer->SetLowerDecreaseRatio(0.01);
  itkOptimizer->SetMiddleDecreaseRatio(0.3);
  itkOptimizer->SetUpperDecreaseRatio(0.5);
 
  {
  itkOptimizer->SetInitialPosition( initialPosition );

  try 
    {
    itkOptimizer->StartOptimization();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cout << "Exception thrown ! " << std::endl;
    std::cout << "An error ocurred during Optimization" << std::endl;
    std::cout << "Location    = " << e.GetLocation()    << std::endl;
    std::cout << "Description = " << e.GetDescription() << std::endl;
    return EXIT_FAILURE;
    }

  finalPosition = itkOptimizer->GetCurrentPosition();
  std::cout << "Solution        = (";
  std::cout << finalPosition[0] << "," ;
  std::cout << finalPosition[1] << ")" << std::endl;  
  std::cout << "Number of iterations = " << itkOptimizer->GetCurrentIteration() << std::endl;

  //
  // check results to see if it is within range
  //
  pass = true;
  for( unsigned int j = 0; j < 2; j++ )
    {
    if( vnl_math_abs( finalPosition[j] - trueParameters[j] ) > 0.01 )
      {
      pass = false;
      }
    }

  if( !pass )
    {
    std::cout << "Test failed." << std::endl;
    return EXIT_FAILURE;
    }

  }

  // Make the function -1 of itself and check maximization
  {
     std::cout << "Trying maximization" << std::endl;
  itkOptimizer->SetInitialPosition( initialPosition );
  costFunction->SetFactor(-1.0);
  itkOptimizer->MaximizeOn();


  try 
    {
    itkOptimizer->StartOptimization();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cout << "Exception thrown ! " << std::endl;
    std::cout << "An error ocurred during Optimization" << std::endl;
    std::cout << "Location    = " << e.GetLocation()    << std::endl;
    std::cout << "Description = " << e.GetDescription() << std::endl;
    return EXIT_FAILURE;
    }

  finalPosition = itkOptimizer->GetCurrentPosition();
  std::cout << "Solution        = (";
  std::cout << finalPosition[0] << "," ;
  std::cout << finalPosition[1] << ")" << std::endl;  
  std::cout << "Number of iterations = " << itkOptimizer->GetCurrentIteration() << std::endl;

  //
  // check results to see if it is within range
  //
  pass = true;
  for( unsigned int j = 0; j < 2; j++ )
    {
    if( vnl_math_abs( finalPosition[j] - trueParameters[j] ) > 0.01 )
      {
      pass = false;
      }
    }

  if( !pass )
    {
    std::cout << "Test failed." << std::endl;
    return EXIT_FAILURE;
    }

  }

  //
  // Verify that the optimizer doesn't run if the 
  // number of iterations is set to zero.
  //
  {
  itkOptimizer->SetNumberOfIterations( 0 );
  itkOptimizer->SetInitialPosition( initialPosition );

  try 
    {
    itkOptimizer->StartOptimization();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cout << excp << std::endl;
    return EXIT_FAILURE;
    }

  if( itkOptimizer->GetCurrentIteration() > 0 )
    {
    std::cerr << "The optimizer is running iterations despite of ";
    std::cerr << "having a maximum number of iterations set to zero" << std::endl;
    return EXIT_FAILURE;
    }
  }

  std::cout << "TEST PASSED !" << std::endl;
  return EXIT_SUCCESS;

}

int main(int argc , char ** argv)
{
 return itkGradientDescentTrustRegionOptimizerTest(argc,argv); 
}
