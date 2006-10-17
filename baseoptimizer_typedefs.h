// Since in C++ classes do not inherit base class typedefs, let's reuse by using the preprocessor
    typedef typename BaseOptimizer::Function Function;
    typedef typename BaseOptimizer::Value Value;
    typedef typename BaseOptimizer::Tolerance Tolerance;
    typedef typename BaseOptimizer::PSet PSet;
    typedef typename BaseOptimizer::Parameter Parameter;
    typedef typename BaseOptimizer::Parameter::Type ParamType;
