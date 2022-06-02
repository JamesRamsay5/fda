\name{NEWS}
\alias{NEWS}
\title{Recent changes to the fda package}
\description{
  \itemize{
    \item{Changes in version fda_6.0.3 2022-05-02:}{
      \itemize{
        \item{Landmark registration:}{Landmark registration using function 
          \code{landmarkreg} can no longer be done by using
          function \code{smooth.basis} instead of function \code{smooth.morph}.  The 
          warping function must be strictly monotonic, and we have found that using 
          \code{smooth.basis} too often violates this monotonicity constraint.  Function 
          \code{smooth.morph} ensures monotonicity and in most applications takes negligible
          computer time to do so.
        }
        \item{PACE in fda:}{
          Function \code{pcaPACE} arries out a functional PCA with regularization from the 
          estimate of the covariance surface.

          Function \code{scoresPACE} estimates functional Principal Component 
          scores through Conditional Expectation (PACE).
        }
        \item{Further changes to \code{smooth.morph} and \code{landmarkreg:}}{
          \code{Smooth.morph} estimates a warping function when the target of the fit by
          registration is a functional data object.  This function has been extended 
          to work when the target for the fit and the fitted functions have different
          ranges or domains.  The warping also maps each boundary into its target
          boundary.  Simiarly \code{landmarkreg} uses a small number of discrete
          values to define the warping, and how has an extra argument, \code{x0lim},
          that defines the range of the target domain.  Since it defaults to the 
          range of the warped domain, it continues to work if not used and the
          domains have the same range.
        }
        \item{Surprisal smoothing:}{This function works with multinomial data that
          evolve over a continuum, such as the value of a latent variable in 
          psychometrics.  A multinomial observation consists of a set of 
          probabilities that are in the open interval (0,1) and sum to one.
          The surprisal value S(P_m) corresponding to a probabity P_m is 
          -log_M(P_m), where M is the number of probabities and is the base of
          the logarithm.  The inverse function is P(S_m) = M^(-S_m).
          
          Surprisal is also known as "self-information" in the field of information
          theory.  It has the characteristics of a true metric:  Surprisals can be 
          added, multiplied by positive numbers, and the difference between two 
          surprisal values mean the same thing everywhere along the information.
          continuum. The unit of the metric is called the "M-bit", the
          generalization of the familiar "bit" or "2-bit" for binary data.
          The metric property is not possessed by so-called latent  
          variables because they can be arbitrarily monotonically transformed.
          
          Smoothing surprisal data is much easier and faster than smoothing
          probabilities since surprisal values are only constrained to be 
          non-negative and are otherwise unbounded.  
          
          The function \code{smooth.surp} estimates smooth curves which fit a set of 
          surprisal values and which also satisfy the constraint that their 
          probability versions sum to one.
        }
        \item{Improvements in iterative optimisation:}{
          Many functions in the fda package optimize a fitting criterion
          iteratively.  Function \code{smooth.monotone} is an example.
          The optimisation algorithm used was a rather early design,
          and many improvements have since been made.  In most of our
          optimisations, we have switched to the algorithm to be found
          in Press, Teukolsky, Vetterling and Flannery Numerical Recipes
          volumes.  We have noticed a bit improvement in speed, are in 
          the process of upgrading all of our optimisers using this 
          approach.
        }
      }
    }
    \item{Changes in version fda_5.5.0 2021-10-28:}{
      \itemize{
        \item{Smooth and constrained curves:}{
          Many data smoothing situations require that the smooth curves satisfy some constraints. 
    
          Take function \code{smooth.monotone.R} for example. Its curves are either strictly 
          increasing or strictly decreasing, even though the data are not.  This is the case in 
          modelling human growth, where we can reasonably assume that daily or monthly measurements 
          will reflect a trend that increases everywhere.  
    
          Function \code{smooth.morph.R}, which plays an important role in curve registration, 
          adds the additional constraint that the domain limits are mapped exactly into the range 
          limits.
    
          In this version two new constrained curves are introduced.  Nonsingular multinomial 
          probability vectors contain nonzero probabilities that sum to zero.  A simple 
          transformation of these probabilities, $S = -log(P)$, converts probabilities into what 
          is often called surprisal.  Surprisal is a measure of information where the unit of 
          measurement is the M-bit, where $M$ is the length of the multinomial vector.  
          Information measured in this way can be added and subtracted, and fixed differences 
          mean the same thing anywhere along the surprisal continuum, which is positive with 
          an origin at 0.  Probability 1 corresponds to surprisal 0, and a very small 
          probability produces a very large positive surprisal.  Probabilities 0.05 and 0.01 
          correspond to 2-bit surprisals 4.3 and 6.1, respectively.
    
          Probability curves result if the probabilities change with over continuous scale, 
          often called a latent variable in statistics.  The corresponding surprisal curves 
          satisfy the constraint at any index value $log(sum(M^-S)) = 0.$  The unbounded nature 
          of surprisal curves plus their metric property render them much easier to work with    
          computationally, as well having the metric property.
    
          Functions smooth.surp.R and error sum of squares fit function surp.fit.R are added 
          in this version in order to support a package \code{TestGardener} that analyzes 
          choice or psychometric data.
        }
      }
    }
  }
}
