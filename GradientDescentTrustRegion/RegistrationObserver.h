
//  The following section of code implements a Command observer
//  that will control the modification of optimizer parameters
//  at every change of resolution level.
//
template <typename TRegistration>
class RegistrationInterfaceCommand : public itk::Command
{
public:
  typedef  RegistrationInterfaceCommand   Self;
  typedef  itk::Command                   Superclass;
  typedef  itk::SmartPointer<Self>        Pointer;
  itkNewMacro( Self );
protected:
  RegistrationInterfaceCommand() {};
public:
  typedef   TRegistration                              RegistrationType;
  typedef   RegistrationType *                         RegistrationPointer;
  typedef   itk::RegularStepGradientDescentOptimizer   OptimizerType1;
  typedef   OptimizerType1 *                           Optimizer1Pointer;
  typedef   itk::GradientDescentTrustRegionOptimizer   OptimizerType2;
  typedef   OptimizerType2 *                           Optimizer2Pointer;

  void Execute(itk::Object * object, const itk::EventObject & event)
  {
    if( !(itk::IterationEvent().CheckEvent( &event )) )
      {
      return;
      }
    RegistrationPointer registration =
                        dynamic_cast<RegistrationPointer>( object );
    Optimizer1Pointer optimizer1 = dynamic_cast< Optimizer1Pointer >(
                       registration->GetOptimizer() );
    Optimizer2Pointer optimizer2 = dynamic_cast< Optimizer2Pointer >(
                       registration->GetOptimizer() );

    std::cout << "-------------------------------------" << std::endl;
    std::cout << "MultiResolution Level : "
              << registration->GetCurrentLevel()  << std::endl;
    std::cout << std::endl;

    if(optimizer1) {
       if ( registration->GetCurrentLevel() == 0 )
         {
         optimizer1->SetMaximumStepLength( 16.00 );
         optimizer1->SetMinimumStepLength(  0.01 );
         }
       else
         {
         optimizer1->SetMaximumStepLength( optimizer1->GetMaximumStepLength() / 4.0 );
         optimizer1->SetMinimumStepLength( optimizer1->GetMinimumStepLength() / 10.0 );
         std::cout << "Iterations: " << optimizer1->GetCurrentIteration() << std::endl;
         }
    }
    if(optimizer2) {
       if ( registration->GetCurrentLevel() == 0 )
         {
         optimizer2->SetInitialStepLength( 16.00 );
         optimizer2->SetMinimumStepLength(  0.01 );
         }
       else
         {
         optimizer2->SetInitialStepLength( optimizer2->GetInitialStepLength() / 4.0 );
         optimizer2->SetMinimumStepLength( optimizer2->GetMinimumStepLength() / 10.0 );
         std::cout << "Iterations: " << optimizer2->GetCurrentIteration() << std::endl;
         }
    }
  }
  void Execute(const itk::Object * , const itk::EventObject & )
    { return; }
};
