/*=========================================================================

  This class has been extracted from the TAPIR registration software suite
  Copyright Rupert Brooks, McGill University 2006-2011

  It should be replaced with an officially published variant.

=========================================================================*/
#ifndef _tapirTrustRegionOptimizer_txx
#define _tapirTrustRegionOptimizer_txx

#include "tapirTrustRegionOptimizer.h"
#include "itkCommand.h"
#include "itkEventObject.h"
#include "vnl/vnl_math.h"
#include "itkImage.h"
//#include "../metrics/itkAnytimeImageToImageMetric.h"
//#include "../transforms/itkTransformFacade.h"

//#define TELLALL

namespace tapir
{

/**
 * Constructor
 */
TrustRegionOptimizer
::TrustRegionOptimizer()
{

  itkDebugMacro("Constructor");
      
  m_MaximumStepLength = 1.0;
  m_MinimumStepLength = 1e-3;
  m_GradientMagnitudeTolerance = 1e-6;
  m_NumberOfIterations = 100;
  m_CurrentIteration   =   0;
  m_Value = 0;
  m_Maximize = false;
  m_CostFunction = 0;
  m_CurrentStepLength   =   0;
  m_StopCondition = MaximumNumberOfIterations;
  m_Gradient.Fill( 0.0f );
  m_PreviousGradient.Fill( 0.0f );
  m_IdentityParameters.Fill(0.0);
  m_Accuracy=1.0;
  m_LowerDecreaseRatio=0.001;
  m_MiddleDecreaseRatio=0.1;
  m_UpperDecreaseRatio=0.75;
}


/**
 * Start the optimization
 */
void
TrustRegionOptimizer
::StartOptimization( void )
{

  itkDebugMacro("StartOptimization");

  m_CurrentStepLength         = m_MaximumStepLength;
  m_CurrentIteration          = 0;

    const unsigned int spaceDimension = 
    m_CostFunction->GetNumberOfParameters();


  m_Gradient = DerivativeType( spaceDimension );
  m_PreviousGradient = DerivativeType( spaceDimension );
  m_Gradient.Fill( 0.0f );
  m_PreviousGradient.Fill( 0.0f );

  ParametersType currentPosition = GetInitialPosition();
  this->SetCurrentPosition( currentPosition );

    if(!m_ScalesInitialized)
    {
    const unsigned int numberOfParameters = 
      m_CostFunction->GetNumberOfParameters();

    ScalesType scales( numberOfParameters );
    scales.Fill( 1.0f );
    SetScales( scales );
    m_ScalesInitialized = true;
	} else {
      scales = this->GetScales();
	}
  // Make sure the scales have been set properly
  if (scales.size() != spaceDimension)
    {
    itkExceptionMacro(<< "The size of Scales is "
                      << scales.size()
                      << ", but the NumberOfParameters for the CostFunction is "
                      << spaceDimension
                      << ".");
    }

	  if( this->m_Maximize ) 
    {
    direction = 1.0;
    }
  else 
    {
    direction = -1.0;
    }

  transformedGradient.SetSize( spaceDimension );
  
  this->EvaluateCostFunction();

	force=0;
  this->ResumeOptimization();

}





/**
 * Resume the optimization
 */
void
TrustRegionOptimizer
::ResumeOptimization( void )
{
  
  itkDebugMacro("ResumeOptimization");

  m_Stop = false;

  this->InvokeEvent( itk::StartEvent() );

  while( !m_Stop ) 
    {

    

    if( m_Stop )
      {
      break;
      }

   if( scaledGradientMagnitude < m_GradientMagnitudeTolerance ) 
    {
    m_StopCondition = GradientMagnitudeTolerance;
    this->StopOptimization();
    break;
    }
    

  if( m_CurrentStepLength < m_MinimumStepLength/10 )
    {
    m_StopCondition = StepTooSmall;
    this->StopOptimization();
    break;
    }

  if( m_CurrentIteration == m_NumberOfIterations )
      {
      m_StopCondition = MaximumNumberOfIterations;
      this->StopOptimization();
      break;
      }

	  this->AdvanceOneStep();

    
    }
    

}





/**
 * Stop optimization
 */
void
TrustRegionOptimizer
::StopOptimization( void )
{
  itkDebugMacro("StopOptimization");
  m_Stop = true;
  this->InvokeEvent( itk::EndEvent() );
}




/**
 * Advance one Step following the gradient direction
 */
void
TrustRegionOptimizer
::AdvanceOneStep( void )
{ 

	itkDebugMacro("AdvanceOneStep");
	m_PreviousGradient=m_Gradient;
	m_PreviousValue=m_Value;
		previousTransformedGradient=transformedGradient;

	ParametersType oldPosition = this->GetCurrentPosition();
	const unsigned int  spaceDimension =
		m_CostFunction->GetNumberOfParameters();

	DerivativeType step( spaceDimension );


	// MAY HAVE TO CHANGE BACK
	//const double factor = 
	//	direction * m_CurrentStepLength / gradientMagnitude;
	double transformedGradientMagnitude=0.0;
	for(unsigned int i=0;i<spaceDimension;i++) transformedGradientMagnitude+=transformedGradient[i]*transformedGradient[i]*scales[i];
	const double factor = m_CurrentStepLength / sqrt(transformedGradientMagnitude);

	//std::cout<<"factor : "<<factor<<std::endl;
	double pred,cred,rho;
	pred=0;
		// MIGHT HAVE TO CHANGE BACK
		double fct=factor;
		fct*=direction;
		for(unsigned int i=0; i<spaceDimension; i++)
		{
			step[i]=transformedGradient[i]*fct;
		}
		pred=0.0;
		for(unsigned int i=0;i<m_Gradient.Size();++i) {
			pred+=step[i]*m_Gradient[i];
		}


	// This method StepAlongGradient() will 
	// be overloaded in non-vector spaces
	this->StepAlongGradient( step );
	this->EvaluateCostFunction();

	cred=(m_PreviousValue-m_Value)*direction;
	rho=cred/pred;
   bool atLimit=true;
	if(rho<m_LowerDecreaseRatio) {
		// Unacceptable decrease.  Roll back the step, Shrink stepsize and retry

		this->SetCurrentPosition( oldPosition );
		m_Gradient=m_PreviousGradient;
		m_Value=m_PreviousValue;
			transformedGradient=previousTransformedGradient;

			m_CurrentStepLength*=0.1;
			force=0;
		double magnitudeSquare = 0;
		for(unsigned int dim=0; dim<spaceDimension; dim++)
		{
			const double weighted = transformedGradient[dim];
			magnitudeSquare += weighted * weighted;
		}

		gradientMagnitude = vcl_sqrt( magnitudeSquare );

		magnitudeSquare = 0;
		for(unsigned int dim=0; dim<spaceDimension; dim++)
		{
			const double weighted = m_Gradient[dim]/scales[dim];
			magnitudeSquare += weighted * weighted;
		}

		scaledGradientMagnitude = vcl_sqrt( magnitudeSquare );

	} else {
		// this step will be accepted
		
		if (rho<m_MiddleDecreaseRatio) {
		// kind of crappy, shrink step by half
		m_CurrentStepLength*=0.5;
		force=0;
      this->InvokeEvent( itk::IterationEvent() );
		m_CurrentIteration++;

		} else if(rho>m_UpperDecreaseRatio && atLimit ) {
		m_CurrentStepLength*=2;
		force=0;
      this->InvokeEvent( itk::IterationEvent() );
		m_CurrentIteration++;

	} else {
		// Everything is Cool... carry on
		force=0;
      this->InvokeEvent( itk::IterationEvent() );
		m_CurrentIteration++;
	}
	}
}

void
TrustRegionOptimizer
::EvaluateCostFunction() {
	ParametersType currentPosition = this->GetCurrentPosition();
	const unsigned int  spaceDimension =
		m_CostFunction->GetNumberOfParameters();

			m_CostFunction->GetValueAndDerivative(currentPosition, 
				m_Value, m_Gradient);
			for(unsigned int i = 0;  i < spaceDimension; i++)
			{
				transformedGradient[i]  = m_Gradient[i] / scales[i];    
			}

	double magnitudeSquare = 0;
	for(unsigned int dim=0; dim<spaceDimension; dim++)
	{
		const double weighted = transformedGradient[dim];
		magnitudeSquare += weighted * weighted;
	}

	gradientMagnitude = vcl_sqrt( magnitudeSquare );
		magnitudeSquare = 0;
		for(unsigned int dim=0; dim<spaceDimension; dim++)
		{
			const double weighted = m_Gradient[dim]/scales[dim];
			magnitudeSquare += weighted * weighted;
		}

		scaledGradientMagnitude = vcl_sqrt( magnitudeSquare );
}

void
TrustRegionOptimizer
::PrintSelf( std::ostream& os, itk::Indent indent ) const
{
  Superclass::PrintSelf(os,indent);
  os << indent << "MaximumStepLength: "
     << m_MaximumStepLength << std::endl;
  os << indent << "MinimumStepLength: "
     << m_MinimumStepLength << std::endl;
  os << indent << "GradientMagnitudeTolerance: "
     << m_GradientMagnitudeTolerance << std::endl;
  os << indent << "NumberOfIterations: "
     << m_NumberOfIterations << std::endl;
  os << indent << "CurrentIteration: "
     << m_CurrentIteration   << std::endl;
  os << indent << "Value: "
     << m_Value << std::endl;
  os << indent << "Maximize: "
     << m_Maximize << std::endl;
  if (m_CostFunction)
    {
    os << indent << "CostFunction: "
       << &m_CostFunction << std::endl;
    }
  else
    {
    os << indent << "CostFunction: "
       << "(None)" << std::endl;
    }
  os << indent << "CurrentStepLength: "
     << m_CurrentStepLength << std::endl;
  os << indent << "StopCondition: "
     << m_StopCondition << std::endl;
  os << indent << "Gradient: "
     << m_Gradient << std::endl;
}


void
TrustRegionOptimizer
::StepAlongGradient( const DerivativeType & step )
{ 

  itkDebugMacro(<<step );

  const unsigned int spaceDimension =
    m_CostFunction->GetNumberOfParameters();
  DerivativeType innerStep(spaceDimension);
  ParametersType newPosition( spaceDimension );
  ParametersType currentPosition = this->GetCurrentPosition();

  for(unsigned int j=0; j<spaceDimension; j++)
    {
    newPosition[j] = currentPosition[j] + step[j];
    }
  itkDebugMacro(<<"new position = " << newPosition );

  this->SetCurrentPosition( newPosition );

}



} // end namespace itk

#endif
