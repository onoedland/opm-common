/*
Adapts the RegulaFalsi solver to work with Evaluations (instead of doubles)

At least during the testing stage, we keep this in a separate file
from RootFinders.hpp!

TO DO (maybe): Also include version with initial guess??
*/

#ifndef OPM_ROOTFINDERS_EVALUATION_HEADER
#define OPM_ROOTFINDERS_EVALUATION_HEADER

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>

#include <string>
#include <algorithm>
#include <limits>
#include <cmath>
#include <iostream>

// ============================================================================
#include <exception>
namespace Opm{

    class ParameterOutsideOfValidRangeException : public std::domain_error{
        public:

            ParameterOutsideOfValidRangeException() : std::domain_error("Parameter is outside of valid range...") { };
            ParameterOutsideOfValidRangeException(const std::string& error_msg) : std::domain_error(error_msg) { };

    };

    class TooManyIterationsException : public std::domain_error{
        public:

            TooManyIterationsException() : std::domain_error("Number of iterations exceeded the maximum number without finding a solution...") { };
            TooManyIterationsException(const std::string& error_msg) : std::domain_error(error_msg) { };

    };
}
// ============================================================================

namespace Opm
{

    template <typename Evaluation>
    class EvaluationRegulaFalsi
    {
        public:

        /// Implements a modified regula falsi method as described in
        /// "Improved algorithms of Illinois-type for the numerical
        /// solution of nonlinear equations"
        /// by J. A. Ford.
        /// Current variant is the 'Pegasus' method.
        template <class Functor>
        inline static Evaluation solve(const Functor& f,
                                       const Evaluation& a,
                                       const Evaluation& b,
                                       const int max_iter,
                                       const double tolerance,
                                       int& iterations_used)
        {
            using namespace std;
            const double macheps = numeric_limits<double>::epsilon();
            const double max_ab = max(fabs(a.value()), fabs(b.value()));
            const double eps = tolerance + macheps * max(max_ab, 1.0);

            Evaluation x0(a);
            Evaluation x1(b);
            Evaluation f0 = f(x0);

            const double epsF = tolerance + macheps * max(abs(f0.value()), 1.0);
            if (Opm::abs(f0) < epsF) {
                return x0;
            }

            Evaluation f1 = f(x1);
            if (abs(f1.value()) < epsF) {
                return x1;
            }
            if (f0*f1 > 0.0) {
                throw ParameterOutsideOfValidRangeException("Zero not bracketed!!");
                //return ErrorPolicy::handleBracketingFailure(a, b, f0, f1);
            }
            iterations_used = 0;
            // In every iteraton, x1 is the last point computed, and x0 is
            // the last point computed that makes it a bracket.
            while (abs(x1 - x0).value() >= eps) {

                Evaluation xnew = regulaFalsiStep(x0, x1, f0, f1);
                Evaluation fnew = f(xnew);
                ++iterations_used;
                if (iterations_used > max_iter) {
                    throw TooManyIterationsException("");
                    //return ErrorPolicy::handleTooManyIterations(x0, x1, max_iter);
                }

                if (abs(fnew.value()) < epsF) {
                    return xnew;
                }

                // Now we must check which point we must replace.
                if ((fnew > 0.0) == (f0 > 0.0)) {
                    // We must replace x0.
                    x0 = x1;
                    f0 = f1;
                }
                else {
                // We must replace x1, this is the case where
                // the modification to regula falsi kicks in,
                // by modifying f0.
                // 1. The classic Illinois method
                //                  const double gamma = 0.5;
                // @afr: The next two methods do not work??!!?
                // 2. The method called 'Method 3' in the paper.
                //                  const double phi0 = f1/f0;
                //                  const double phi1 = fnew/f1;
                //                  const double gamma = 1.0 - phi1/(1.0 - phi0);
                // 3. The method called 'Method 4' in the paper.
                //                  const double phi0 = f1/f0;
                //                  const double phi1 = fnew/f1;
                //                  const double gamma = 1.0 - phi0 - phi1;
                //                  cout << "phi0 = " << phi0 <<" phi1 = " << phi1 <<
                //                  " gamma = " << gamma << endl;
                // 4. The 'Pegasus' method
                    const Evaluation gamma = f1 / (f1 + fnew);
                    f0 *= gamma;
                }
                x1 = xnew;
                f1 = fnew;
            } // end while loop

            return 0.5*(x0 + x1);

        } // end solve()


    private:

        inline static Evaluation regulaFalsiStep(const Evaluation a,
                                                 const Evaluation b,
                                                 const Evaluation fa,
                                                 const Evaluation fb)
        {
            assert(fa*fb < 0.0);
            return (b*fa - a * fb) / (fa - fb);
        }
    };

} // end namespace Opm

#endif
