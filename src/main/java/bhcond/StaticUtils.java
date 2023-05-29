package bhcond;

import static bayonet.math.SpecialFunctions.*;

public class StaticUtils
{
  public static double betaBinomialLogDensity(double alpha, double beta, int Lr, int Rr, int Er) {
    return 
        lnGamma(alpha + beta) 
      - lnGamma(alpha) 
      - lnGamma(beta) 
      + lnGamma(Er + alpha) 
      + lnGamma( (Lr * Rr) + beta - Er) 
      - lnGamma( (Lr * Rr) + alpha + beta);
  }
  
  public static double truncatedGeometricLogDensity(double rho, int Lr, int Rr, int Er) {
    return 
        lnGamma(Er + 1) 
      + lnGamma( (Lr * Rr) - Er + 1) 
      - lnGamma( (Lr * Rr) + 1)
      + Er * Math.log(rho) + Math.log(1 - rho)
      - Math.log(1 - Math.pow(rho, (Lr*Rr + 1)) );
  }
}
