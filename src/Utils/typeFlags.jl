
# Flags representing function types to dispatch on
abstract type FunctionType end
struct Dynamics <: FunctionType end
struct Cost <: FunctionType end 
struct Algebraic <: FunctionType end

# Flags representing Jacobian types to dispatch on
abstract type JacobianType end 
struct State <: JacobianType end 
struct Control <: JacobianType end 
struct Time <: JacobianType end 
struct Static <: JacobianType end

# Transcription types 
abstract type TranscriptionType end
struct ImplicitRK <: TranscriptionType end

# NLP Solver types
abstract type NLPSolver end
struct IPOPT <: NLPSolver end