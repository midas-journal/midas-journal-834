/*=========================================================================

  Copyright (c) 2005-2011 Rupert Brooks. 

  This code is part of the insight journal submission 
  "A Gradient Descent Trust Region Optimizer"

  This software is distributed WITHOUT ANY WARRANTY; without even 
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
  PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkGradientDescentTrustRegionOptimizer_h
#define __itkGradientDescentTrustRegionOptimizer_h

#include "itkSingleValuedNonLinearOptimizer.h"
#include <sstream>

namespace itk
{
  
/** \class GradientDescentTrustRegionOptimizer
 * \brief Performs gradient descent using trust region principles to control 
 *        the step size
 *
 * \ingroup Numerics Optimizers
 */
class ITK_EXPORT GradientDescentTrustRegionOptimizer : 
    public SingleValuedNonLinearOptimizer
{
public:
  /** Standard "Self" typedef. */
  typedef GradientDescentTrustRegionOptimizer          Self;
  typedef SingleValuedNonLinearOptimizer               Superclass;
  typedef SmartPointer<Self>                           Pointer;
  typedef SmartPointer<const Self>                     ConstPointer;
  
  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  
  /** Run-time type information (and related methods). */
  itkTypeMacro( GradientDescentTrustRegionOptimizer, 
                SingleValuedNonLinearOptimizer );
  

  /** Codes of stopping conditions. 
   *  for compatibility, uses same codes as RegularStepGradientDescent
   *  but not all can possibly happen which is why some numbers are skipped
   */
  
  typedef enum {
    GradientMagnitudeTolerance=1,
    StepTooSmall=2,
    CostFunctionError = 4,
    MaximumNumberOfIterations = 5,
    Unknown = 6
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

  /** Set/Get The maximum allowable trust region size. */
  itkSetMacro( MaximumStepLength, double );
  itkGetConstReferenceMacro( MaximumStepLength, double );
  /** Set/Get The starting trust region size. */
  itkSetMacro( InitialStepLength, double );
  itkGetConstReferenceMacro( InitialStepLength, double );
  /** Set/Get smallest allowable trust region size. This is a termination criterion. */
  itkSetMacro( MinimumStepLength, double );
  itkGetConstReferenceMacro( MinimumStepLength, double );
  /** Set/Get Maximum allowable number of iterations. This is a termination criterion. */
  itkSetMacro( NumberOfIterations, unsigned long );
  itkGetConstReferenceMacro( NumberOfIterations, unsigned long );
  /** Set/Get Minimum allowable scaled gradient magnitude. This is a termination criterion. */
  itkSetMacro( GradientMagnitudeTolerance, double );
  itkGetConstReferenceMacro( GradientMagnitudeTolerance, double );
  /** Get the current trust region size. */
  itkGetConstReferenceMacro( CurrentStepLength, double);
  /** Get the current iteration */
  itkGetConstMacro( CurrentIteration, unsigned int );
  /** Get the Stop condition enum */
  itkGetConstReferenceMacro( StopCondition, StopConditionType );
  /** Get the Current value */
  itkGetConstReferenceMacro( Value, MeasureType );
  /** Get the Current Gradient */
  itkGetConstReferenceMacro( Gradient, DerivativeType );
  /** Set/Get the improvement ratio below which the step is rejected */
  itkSetClampMacro( LowerDecreaseRatio, double, 0.0, 1.0 );
  itkGetConstMacro( LowerDecreaseRatio, double); 
  /** Set/Get the improvement ratio below which the step is accepted but the trust region is reduced*/
  itkSetClampMacro( MiddleDecreaseRatio, double, 0.0, 1.0 );
  itkGetConstMacro( MiddleDecreaseRatio, double); 
  /** Set/Get the improvement ratio above which the trust region is enlarged */
  itkSetClampMacro( UpperDecreaseRatio, double, 0.0, 1.0 );
  itkGetConstMacro( UpperDecreaseRatio, double); 
  /** Set/Get the factor by which the trust region is reduced when the step is rejected */
  itkSetClampMacro( RejectedStepDecreaseFactor, double, 0.0, 1.0 );
  itkGetConstMacro( RejectedStepDecreaseFactor, double); 
  /** Set/Get the factor by which the trust region is reduced when the step is poor */
  itkSetClampMacro( PoorStepDecreaseFactor, double, 0.0, 1.0 );
  itkGetConstMacro( PoorStepDecreaseFactor, double); 
  /** Set/Get the factor by which the trust region is enlarged when the step is good */
  itkSetClampMacro( StepIncreaseFactor, double, 1.0, itk::NumericTraits<double>::max() );
  itkGetConstMacro( StepIncreaseFactor, double); 
 
  /** Get the reason for termination as text */
  virtual const std::string GetStopConditionDescription() const;
 
protected:
  GradientDescentTrustRegionOptimizer();
  virtual ~GradientDescentTrustRegionOptimizer() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** Continue solving the trust region subproblem until one step is accepted or termination */
  virtual void AdvanceOneStep( void );
  /** Evaluate cost function and fill internal variables */
  virtual void EvaluateCostFunction( void );
  /** Update the parameters by the step */
  virtual void StepAlongGradient( const DerivativeType&);
  

protected:
  DerivativeType                m_Gradient; 
  double                        m_ScaledGradientMagnitude;
  DerivativeType                m_TransformedGradient;

  bool                          m_Stop;
  bool                          m_Maximize;
  MeasureType                   m_Value;
  double                        m_LowerDecreaseRatio;
  double                        m_MiddleDecreaseRatio;
  double                        m_UpperDecreaseRatio;

  double                        m_GradientMagnitudeTolerance;
  double                        m_MaximumStepLength;
  double                        m_InitialStepLength;
  double                        m_MinimumStepLength;
  double                        m_CurrentStepLength;
  StopConditionType             m_StopCondition;
  unsigned long                 m_NumberOfIterations;
  unsigned long                 m_CurrentIteration;

  double                        m_RejectedStepDecreaseFactor;
  double                        m_PoorStepDecreaseFactor;
  double                        m_StepIncreaseFactor;
  
  std::ostringstream            m_StopConditionDescription;

private:  
  GradientDescentTrustRegionOptimizer(const Self&); //purposely not implemented
  void operator=(const Self&);//purposely not implemented


};

} // end namespace itk



#endif



