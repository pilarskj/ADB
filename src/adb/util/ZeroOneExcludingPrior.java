package adb.util;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.inference.distribution.Prior;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;

@Description("Like Prior, but excludes any zero and one elements of the input parameter from the calculation.")
public class ZeroOneExcludingPrior extends Prior {

    @Override
    public double calculateLogP() {
        Function x = m_x.get();
        if (x instanceof RealParameter || x instanceof IntegerParameter) {
            // test that parameter is inside its bounds
            double l;
            double h;
            if (x instanceof RealParameter) {
                l = ((RealParameter)x).getLower();
                h = ((RealParameter)x).getUpper();
            } else {
                l = ((IntegerParameter)x).getLower();
                h = ((IntegerParameter)x).getUpper();
            }
            for (int i = 0; i < x.getDimension(); i++) {
                double value = x.getArrayValue(i);
                if (value < l || value > h) {
                    logP = Double.NEGATIVE_INFINITY;
                    return Double.NEGATIVE_INFINITY;
                }
            }
        }

        logP = 0.0;
        for (double elementVal : x.getDoubleValues()) {
            if (elementVal == 0.0 || elementVal == 1.0)
                continue;

            logP += dist.logDensity(elementVal);
        }

        if (logP == Double.POSITIVE_INFINITY) {
            logP = Double.NEGATIVE_INFINITY;
        }

        return logP;
    }

}
