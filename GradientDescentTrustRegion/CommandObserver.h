/*=========================================================================

  This is the command observer, common to all the examples, which has been 
  copied out into a separate header for brevity.

=========================================================================*/

template <class TOptimizer>
class CommandIterationUpdate : public itk::Command 
{
public:
  typedef  CommandIterationUpdate<TOptimizer>   Self;
  typedef  itk::Command             Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  itkNewMacro( Self );

protected:
  CommandIterationUpdate() {};

public:

  typedef TOptimizer                          OptimizerType;
  typedef const OptimizerType                         *OptimizerPointer;

  void Execute(itk::Object *caller, const itk::EventObject & event)
  {
    Execute( (const itk::Object *)caller, event);
  }

  void Execute(const itk::Object * object, const itk::EventObject & event)
  {
    OptimizerPointer optimizer = 
                         dynamic_cast< OptimizerPointer >( object );

    if( ! itk::IterationEvent().CheckEvent( &event ) )
      {
      return;
      }
    typename OptimizerType::DerivativeType gradient = optimizer->GetGradient();
    typename OptimizerType::ScalesType     scales   = optimizer->GetScales();

    double magnitude2 = 0.0;

    for(unsigned int i=0; i<gradient.size(); i++)
      {
      const double fc = gradient[i] / scales[i];
      magnitude2 += fc * fc;
      }

    const double gradientMagnitude = vcl_sqrt( magnitude2 );


    std::cout << optimizer->GetCurrentIteration() << " = ";
    std::cout << optimizer->GetValue() << " : ";
    std::cout << gradient << " ( " << gradientMagnitude << " ) : "; 
    std::cout << optimizer->GetCurrentPosition() << std::endl;
  }
   
};


