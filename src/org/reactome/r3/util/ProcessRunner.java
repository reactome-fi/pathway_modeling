/*
 * Created on Feb 22, 2011
 *
 */
package org.reactome.r3.util;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

import org.junit.Test;

/**
 * This class is used to do a process call (e.g. R, Perl or MCL call).
 * @author wgm
 */
public class ProcessRunner {
    
    public ProcessRunner() {
    }
    
    public String[] runRScript(String[] parameters) throws IOException {
        String[] extraParameters = new String[parameters.length + 1];
        extraParameters[0] = "Rscript";
        for (int i = 0; i < parameters.length; i++)
            extraParameters[i + 1] = parameters[i];
        return runScript(extraParameters);
    }
    
    /**
     * Run a script. The actual script name should be the in the passed
     * parameters.
     * @param parameters
     * @return two element String array. The fist element should be generated from
     * System.out, while the second from System.err.
     * @throws IOException
     */
    public String[] runScript(String[] parameters) throws IOException {
        Runtime runtime = Runtime.getRuntime();
        Process process = runtime.exec(parameters);
        InputStream is = process.getInputStream();
        String output = output(is);
        is.close();
        is = process.getErrorStream();
        String error = output(is);
        is.close();
        String[] rtn = new String[2];
        rtn[0] = output;
        rtn[1] = error;
        return rtn;
    }
    
    /**
     * Helper method to output the result
     * @param is
     * @return
     * @throws IOException
     */
    private String output(InputStream is) throws IOException {
        InputStreamReader isr = new InputStreamReader(is);
        BufferedReader br = new BufferedReader(isr);
        String line = null;
        StringBuilder builder = new StringBuilder();
        while ((line = br.readLine()) != null) {
            if (line.startsWith("Loading required package"))
                continue; // Escape these information lines
            builder.append(line).append("\n");
        }
        br.close();
        isr.close();
        return builder.toString();
    }
    
    /**
     * This method is used to test run a R script.
     */
    @Test
    public void testRscript() {
        String dirName = "/Users/wgm/Documents/EclipseWorkspace/caBigR3/RSource/";
        String rscriptName = dirName + "CGISurvivalAnalysis.R";
        String dataDirName = "datasets/TCGA/OvarianCancer/";
        String scoreFileName = dataDirName + "SamplesToExomeMutationModules_3More_091310.txt";
        String clinFileName = dataDirName + "data_031910/Batches9-22_tcga_OV_clinical_csv.2010-01-25-noDates.txt";
        String plotFileName = dirName + "JavaTest.pdf";
        String method = "coxph";
        method = "kaplan-meier";
        Runtime runtime = Runtime.getRuntime();
        try {
            // Try coxph analysis
            String[] parameters = new String[] {
                    "Rscript",
                    rscriptName,
                    scoreFileName,
                    clinFileName,
                    method,
                    "4",
                    plotFileName,
                    "PDF"
            };
            String[] output = new ProcessRunner().runScript(parameters);
            System.out.println("Output: \n" + output[0]);
            if (output[1].length() > 0)
                System.out.println("\nError: \n" + output[1]);
        }
        catch(Exception e) {
            e.printStackTrace();
        }
    }
}
