/*
 * Created on Apr 17, 2014
 *
 */
package org.reactome.pathway.factorgraph;

import java.util.ArrayList;
import java.util.List;

/**
 * This class is used to caluclate IPAs.
 * @author gwu
 *
 */
public class IPACalculator {
    
    /**
     * Calculate IPA values for factor graph based data analysis.
     * @param priorProbs
     * @param postProbs
     * @return
     */
    public static double calculateIPA(List<Double> priorProbs,
                                      List<Double> postProbs) {
        if (priorProbs == null || postProbs == null || priorProbs.size() < 3 || postProbs.size() < 3)
            return 0.0d;
        List<Double> ratios = new ArrayList<Double>(3);
        for (int i = 0; i < 3; i++) {
            double ratio = calculateLogRatio(priorProbs.get(i),
                                             postProbs.get(i));
            ratios.add(ratio);
        }
        return calculateIPA(ratios);
    }
    
    public static double calculateIPA(double[] prior,
                                      double[] posterior) {
        List<Double> priorList = new ArrayList<Double>();
        for (int i = 0; i < prior.length; i++)
            priorList.add(prior[i]);
        List<Double> posteriorList = new ArrayList<Double>();
        for (int i = 0; i < posterior.length; i++)
            posteriorList.add(posterior[i]);
        return calculateIPA(priorList, posteriorList);
    }
    
    private static double calculateIPA(List<Double> ratios) {
        double down = ratios.get(0);
        double normal = ratios.get(1);
        double up = ratios.get(2);
        if (up > down && up > normal)
            return up;
        if (down > up && down > normal)
            return -down;
        return 0;
    }
    
    public static double calculateLogRatio(double priorValue,
                                           double postValue) {
        double tmp = postValue / (1.0d - postValue) * (1.0d - priorValue) / priorValue;
        // The following check is used to avoid infinity
        if (tmp < Double.MIN_EXPONENT)
            tmp = Double.MIN_EXPONENT; 
        if (tmp > Double.MAX_EXPONENT)
            tmp = Double.MAX_EXPONENT;
        double value = Math.log10(tmp);
        return value;
    }
    
}
