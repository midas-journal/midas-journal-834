/*=========================================================================

  This code is part of the Insight Journal submission
  "A Gradient Descent Trust-Region Optimizer"

  Copyright (c) 2005-2011 Rupert Brooks.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkGradientDescentTrustRegionOptimizer_txx
#define _itkGradientDescentTrustRegionOptimizer_txx

#include "itkGradientDescentTrustRegionOptimizer.h"
#include "itkCommand.h"
#include "itkEventObject.h"
#include "vnl/vnl_math.h"
#include "itkImage.h"


namespace itk {

/**
 * Constructor
 */
GradientDescentTrustRegionOptimizer
::GradientDescentTrustRegionOptimizer() {

   itkDebugMacro( "Constructor" );

   m_MaximumStepLength          = 64.0;
   m_InitialStepLength          = 1.0;
   m_MinimumStepLength          = 1e-3;
   m_GradientMagnitudeTolerance = 1e-6;
   m_NumberOfIterations         = 100;
   m_CurrentIteration           = 0;
   m_Value                      = 0;
   m_Maximize                   = false;
   m_CostFunction               = 0;
   m_CurrentStepLength          = 0;
   m_StopCondition              = Unknown;
   m_Gradient.Fill( 0.0f );
   m_LowerDecreaseRatio=0.001;
   m_MiddleDecreaseRatio=0.1;
   m_UpperDecreaseRatio=0.75;
   m_RejectedStepDecreaseFactor = 0.1;
   m_PoorStepDecreaseFactor = 0.5;
   m_StepIncreaseFactor = 2.0;


}


/**
 * Start the optimization
 */
void
GradientDescentTrustRegionOptimizer
::StartOptimization( void ) {

   itkDebugMacro( "StartOptimization" );

   m_CurrentStepLength         = m_InitialStepLength;
   m_CurrentIteration          = 0;

   m_StopCondition = Unknown;
   m_StopConditionDescription.str( "" );
   m_StopConditionDescription << this->GetNameOfClass() << ": ";

   const unsigned int spaceDimension = m_CostFunction->GetNumberOfParameters();

   m_Gradient = DerivativeType( spaceDimension );
   m_Gradient.Fill( 0.0f );

   ParametersType currentPosition = GetInitialPosition();
   this->SetCurrentPosition( currentPosition );

   if ( !m_ScalesInitialized ) {
      const unsigned int numberOfParameters =
         m_CostFunction->GetNumberOfParameters();

      ScalesType scales( numberOfParameters );
      scales.Fill( 1.0f );
      SetScales( scales );
      m_ScalesInitialized = true;
   } 
   // Make sure the scales have been set properly
   const unsigned int sz=this->GetScales().size();
   if ( sz != spaceDimension ) {
      itkExceptionMacro( << "The size of Scales is "
                            << sz
                            << ", but the NumberOfParameters for the CostFunction is "
                            << spaceDimension
                            << "." );
   }

   m_TransformedGradient.SetSize( spaceDimension );

   this->EvaluateCostFunction();

   this->ResumeOptimization();

}





/**
 * Resume the optimization
 */
void
GradientDescentTrustRegionOptimizer
::ResumeOptimization( void ) {

   itkDebugMacro( "ResumeOptimization" );

   m_Stop = false;

   this->InvokeEvent( StartEvent() );

   while ( !m_Stop ) {

      if ( m_ScaledGradientMagnitude < m_GradientMagnitudeTolerance ) {
         m_StopCondition = GradientMagnitudeTolerance;
         m_StopConditionDescription << "Gradient magnitude tolerance met after "
            << m_CurrentIteration
            << " iterations. Gradient magnitude (scaled) ("
            << m_ScaledGradientMagnitude
            << ") is less than gradient magnitude tolerance ("
            << m_GradientMagnitudeTolerance
            << ").";

         this->StopOptimization();
         break;
      }

      if ( m_CurrentStepLength < m_MinimumStepLength/10 ) {
         m_StopCondition = StepTooSmall;
         m_StopConditionDescription << "Trust region too small after "
            << m_CurrentIteration
            << " iterations. Current size ("
            << m_CurrentStepLength
            << ") is less than minimum size ("
            << m_MinimumStepLength
            << ").";
         this->StopOptimization();
         break;
      }

      if ( m_CurrentIteration == m_NumberOfIterations ) {
         m_StopCondition = MaximumNumberOfIterations;
         m_StopConditionDescription << "Maximum number of iterations ("
            << m_NumberOfIterations
            << ") exceeded.";
         this->StopOptimization();
         break;
      }

      if ( m_Stop ) { break;  }

      this->AdvanceOneStep();
   }
}





/**
 * Stop optimization
 */
void
GradientDescentTrustRegionOptimizer
::StopOptimization( void ) {
   itkDebugMacro( "StopOptimization" );
   m_Stop = true;
   this->InvokeEvent( EndEvent() );
}




/**
 * Advance one Step following the gradient direction
 */
void
GradientDescentTrustRegionOptimizer
::AdvanceOneStep( void ) {

   const ScalesType & scales=this->GetScales();

   itkDebugMacro( "AdvanceOneStep" );
   const DerivativeType   previousGradient=m_Gradient;
   const MeasureType      previousValue=m_Value;
   const DerivativeType   previousTransformedGradient=m_TransformedGradient;
   const MeasureType      previousScaledGradientMagnitude=m_ScaledGradientMagnitude;

   ParametersType oldPosition = this->GetCurrentPosition();
   const unsigned int  spaceDimension =
      m_CostFunction->GetNumberOfParameters();

   DerivativeType step( spaceDimension );


   double transformedGradientMagnitude=0.0;
   // The trust region radius is computed in the scaled parameter space - if the step is divided by scales, then the parameters are scaled
   // by the square root of scales.
   for ( unsigned int i=0;i<spaceDimension;i++ ) { 
      transformedGradientMagnitude+=m_TransformedGradient[i]*m_TransformedGradient[i]*scales[i]; 
   }
   const double factor = m_CurrentStepLength / sqrt( transformedGradientMagnitude )*( this->m_Maximize ? 1 : -1);

   for ( unsigned int i=0; i<spaceDimension; i++ ) {
      step[i]=m_TransformedGradient[i]*factor;
   }
   double predictedChange=0.0;
   for ( unsigned int i=0;i<m_Gradient.Size();++i ) {
      predictedChange+=step[i]*m_Gradient[i];
   }

   // This method StepAlongGradient() will
   // be overloaded in non-vector spaces
   this->StepAlongGradient( step );
   this->EvaluateCostFunction();

   const double actualChange=( m_Value-previousValue );
   const double rho=actualChange/predictedChange;
   if ( rho<m_LowerDecreaseRatio ) {
      // Unacceptable decrease.  Roll back the step, Shrink stepsize and retry

      this->SetCurrentPosition( oldPosition );
      m_Gradient=previousGradient;
      m_Value=previousValue;
      m_TransformedGradient=previousTransformedGradient;
      m_ScaledGradientMagnitude=previousScaledGradientMagnitude;

      m_CurrentStepLength*=m_RejectedStepDecreaseFactor;
      m_CurrentIteration++;

   } else {
      // this step will be accepted

      if ( rho<m_MiddleDecreaseRatio ) {
         // kind of poor result, shrink step by half
         m_CurrentStepLength*=m_PoorStepDecreaseFactor;
         this->InvokeEvent( IterationEvent() );
         m_CurrentIteration++;

      } else if ( rho>m_UpperDecreaseRatio ) {
         // great result, expand the step
         m_CurrentStepLength=m_CurrentStepLength*m_StepIncreaseFactor;
         m_CurrentStepLength=std::min(m_MaximumStepLength,m_CurrentStepLength*m_StepIncreaseFactor);
         this->InvokeEvent( IterationEvent() );
         m_CurrentIteration++;

      } else {
         // Everything is Cool... carry on
         this->InvokeEvent( IterationEvent() );
         m_CurrentIteration++;
      }
   }
}

void
GradientDescentTrustRegionOptimizer
::EvaluateCostFunction() {
   ParametersType currentPosition = this->GetCurrentPosition();
   const ScalesType & scales=this->GetScales();
   const unsigned int  spaceDimension =
      m_CostFunction->GetNumberOfParameters();

   try {
      m_CostFunction->GetValueAndDerivative( currentPosition,
                                             m_Value, m_Gradient );
   } catch ( ExceptionObject & excp ) {
      m_StopCondition = CostFunctionError;
      m_StopConditionDescription << "Cost function error after "
         << m_CurrentIteration
         << " iterations. "
         << excp.GetDescription();
      this->StopOptimization();
      throw excp;
   }


   for ( unsigned int i = 0;  i < spaceDimension; i++ ) {
      m_TransformedGradient[i]  = m_Gradient[i] / scales[i];
   }
   m_ScaledGradientMagnitude = m_TransformedGradient.magnitude();

}

void
GradientDescentTrustRegionOptimizer
::PrintSelf( std::ostream& os, Indent indent ) const {
   Superclass::PrintSelf( os,indent );
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
   os << indent << "LowerDecreaseRatio: "
      << m_LowerDecreaseRatio << std::endl;
   os << indent << "MiddleDecreaseRatio: "
      << m_MiddleDecreaseRatio << std::endl;
   os << indent << "UpperDecreaseRatio: "
      << m_UpperDecreaseRatio << std::endl;
   os << indent << "Value: "
      << m_Value << std::endl;
   os << indent << "Maximize: "
      << m_Maximize << std::endl;
   if ( m_CostFunction ) {
      os << indent << "CostFunction: "
         << &m_CostFunction << std::endl;
   } else {
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
GradientDescentTrustRegionOptimizer
::StepAlongGradient( const DerivativeType & step ) {

   itkDebugMacro( <<step );

   const unsigned int spaceDimension =
      m_CostFunction->GetNumberOfParameters();
   DerivativeType innerStep( spaceDimension );
   ParametersType newPosition( spaceDimension );
   ParametersType currentPosition = this->GetCurrentPosition();

   for ( unsigned int j=0; j<spaceDimension; j++ ) {
      newPosition[j] = currentPosition[j] + step[j];
   }
   itkDebugMacro( <<"new position = " << newPosition );

   this->SetCurrentPosition( newPosition );

}

const std::string
GradientDescentTrustRegionOptimizer
::GetStopConditionDescription() const {
   return m_StopConditionDescription.str();
}



} // end namespace itk

#endif
