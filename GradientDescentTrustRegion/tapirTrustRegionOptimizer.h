/*=========================================================================

  This class has been extracted from the TAPIR registration software suite
  Copyright Rupert Brooks, McGill University 2006-2011

  It should be replaced with an officially published variant.

=========================================================================*/
#ifndef __itkTrustRegionOptimizer_h
#define __itkTrustRegionOptimizer_h

#include "itkSingleValuedNonLinearOptimizer.h"

namespace tapir
{
  
/** \class TrustRegionOptimizer
 * \brief Implement a gradient descent optimizer
 *
 * \ingroup Numerics Optimizers
 */
class ITK_EXPORT TrustRegionOptimizer : 
   public itk::SingleValuedNonLinearOptimizer
{
public:
  /** Standard "Self" typedef. */
  typedef TrustRegionOptimizer      Self;
  typedef itk::SingleValuedNonLinearOptimizer               Superclass;
  typedef itk::SmartPointer<Self>                           Pointer;
  typedef itk::SmartPointer<const Self>                     ConstPointer;
  
  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  
  /** Run-time type information (and related methods). */
  itkTypeMacro( TrustRegionOptimizer, 
                SingleValuedNonLinearOptimizer );
  

  /** Codes of stopping conditions. */
  typedef enum {
    GradientMagnitudeTolerance=1,
    StepTooSmall,
    ImageNotAvailable,
    SamplesNotAvailable,
    MaximumNumberOfIterations
  } StopConditionType;

  /** Specify whether to minimize or maximize the cost function. */
  itkSetMacro( Maximize, bool );
  itkGetConstReferenceMacro( Maximize, bool );
  itkBooleanMacro( Maximize );

  bool GetMinimize( ) const
  { return !m_Maximize; }

  void SetMinimize(bool v)
  { this->SetMaximize(!v); }

  void    MinimizeOn(void) 
  { SetMaximize( false ); }

  void    MinimizeOff(void) 
  { SetMaximize( true ); }

  /** Start optimization. */
  void    StartOptimization( void );

  /** Resume previously stopped optimization with current parameters.
   * \sa StopOptimization */
  void    ResumeOptimization( void );

  /** Stop optimization.
   * \sa ResumeOptimization */
  void    StopOptimization( void );

  /** Set/Get parameters to control the optimization process. */
  itkSetMacro( MaximumStepLength, double );
  itkSetMacro( MinimumStepLength, double );
  itkSetMacro( NumberOfIterations, unsigned long );
  itkSetMacro( GradientMagnitudeTolerance, double );
  itkGetConstReferenceMacro( CurrentStepLength, double);
  itkGetConstReferenceMacro( MaximumStepLength, double );
  itkGetConstReferenceMacro( MinimumStepLength, double );
  itkGetConstReferenceMacro( NumberOfIterations, unsigned long );
  itkGetConstReferenceMacro( GradientMagnitudeTolerance, double );
  itkGetConstMacro( CurrentIteration, unsigned int );
  itkGetConstReferenceMacro( StopCondition, StopConditionType );
  itkGetConstReferenceMacro( Value, MeasureType );
  itkGetConstReferenceMacro( Gradient, DerivativeType );
  
  itkSetClampMacro( Accuracy, double, 0.0, 1.0 );
  itkGetConstMacro( Accuracy, double); 
  itkSetClampMacro( LowerDecreaseRatio, double, 0.0, 1.0 );
  itkGetConstMacro( LowerDecreaseRatio, double); 
  itkSetClampMacro( MiddleDecreaseRatio, double, 0.0, 1.0 );
  itkGetConstMacro( MiddleDecreaseRatio, double); 
  itkSetClampMacro( UpperDecreaseRatio, double, 0.0, 1.0 );
  itkGetConstMacro( UpperDecreaseRatio, double); 

protected:
  TrustRegionOptimizer();
  virtual ~TrustRegionOptimizer() {};
  void PrintSelf(std::ostream& os, itk::Indent indent) const;
  ScalesType scales;
  double direction,gradientMagnitude,scaledGradientMagnitude;
  DerivativeType transformedGradient;
  DerivativeType previousTransformedGradient;
  double lastComputationLevel;
  int force;

  /** Advance one step following the gradient direction
   * This method verifies if a change in direction is required
   * and if a reduction in steplength is required. */
  virtual void AdvanceOneStep( void );
  virtual void EvaluateCostFunction( void );

  /** Advance one step along the corrected gradient taking into
   * account the steplength represented by factor.
   * This method is invoked by AdvanceOneStep. It is expected
   * to be overrided by optimization methods in non-vector spaces
   * \sa AdvanceOneStep */
  virtual void StepAlongGradient( const DerivativeType&);
  

private:  
  TrustRegionOptimizer(const Self&); //purposely not implemented
  void operator=(const Self&);//purposely not implemented

protected:
  DerivativeType                m_Gradient; 
  DerivativeType                m_PreviousGradient; 
  DerivativeType                m_IdentityParameters; 

  bool                          m_Stop;
  bool                          m_Maximize;
  MeasureType                   m_Value;
  MeasureType                   m_PreviousValue;
  double                        m_LowerDecreaseRatio;
  double                        m_MiddleDecreaseRatio;
  double                        m_UpperDecreaseRatio;

  double                        m_GradientMagnitudeTolerance;
  double                        m_MaximumStepLength;
  double                        m_MinimumStepLength;
  double                        m_CurrentStepLength;
  StopConditionType             m_StopCondition;
  unsigned long                 m_NumberOfIterations;
  unsigned long                 m_CurrentIteration;
  double                        m_Accuracy;

};

} // end namespace itk



#endif



