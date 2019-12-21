import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

public class HMM {
    // Starting
    private final static double START_PROB = 0.5;

    // H
    private final static double H_TO_H = 0.5;
    private final static double H_TO_L = 0.5;

    private final static double A_H = 0.2;
    private final static double C_H = 0.3;
    private final static double G_H = 0.3;
    private final static double T_H = 0.2;

    // L
    private final static double L_TO_L = 0.6;
    private final static double L_TO_H = 0.4;

    private final static double A_L = 0.3;
    private final static double C_L = 0.2;
    private final static double G_L = 0.2;
    private final static double T_L = 0.3;

    void ViterbiAlgorithm(Sequence seq){


        // Split sequence
        String sequence = seq.getSequence();

        // Create matrix for scores
        double[][] matrix = new double[2][sequence.length()];

        // Fill out matrix[1][1] and matrix[1][2]
        if (sequence.charAt(0) == 'A'){
            matrix[0][0] = START_PROB * A_H;
        } else if (sequence.charAt(0) == 'C'){
            matrix[0][0] = START_PROB * C_H;
        } else if (sequence.charAt(0) == 'G'){
            matrix[0][0] = START_PROB * G_H;
        } else if (sequence.charAt(0) == 'T'){
            matrix[0][0] = START_PROB * T_H;
        }

        if (sequence.charAt(0) == 'A'){
            matrix[1][0] = START_PROB * A_L;
        } else if (sequence.charAt(0) == 'C'){
            matrix[1][0] = START_PROB * C_L;
        } else if (sequence.charAt(0) == 'G'){
            matrix[1][0] = START_PROB * G_L;
        } else if (sequence.charAt(0) == 'T'){
            matrix[1][0] = START_PROB * T_L;
        }

        // Fill out the rest of the table
        for (int i = 1; i < sequence.length(); i++){
            // Find values for next rows
            if (sequence.charAt(i) == 'A'){
                matrix[0][i] = A_H * (Math.max(matrix[0][i-1] * H_TO_H, matrix[1][i-1] * L_TO_H));
                matrix[1][i] = A_L * (Math.max(matrix[0][i-1] * H_TO_L, matrix[1][i-1] * L_TO_L));
            } else if (sequence.charAt(i) == 'C'){
                matrix[0][i] = C_H * (Math.max(matrix[0][i-1] * H_TO_H, matrix[1][i-1] * L_TO_H));
                matrix[1][i] = C_L * (Math.max(matrix[0][i-1] * H_TO_L, matrix[1][i-1] * L_TO_L));
            } else if (sequence.charAt(i) == 'G'){
                matrix[0][i] = G_H * (Math.max(matrix[0][i-1] * H_TO_H, matrix[1][i-1] * L_TO_H));
                matrix[1][i] = G_L * (Math.max(matrix[0][i-1] * H_TO_L, matrix[1][i-1] * L_TO_L));
            } else if (sequence.charAt(i) == 'T'){
                matrix[0][i] = T_H * (Math.max(matrix[0][i-1] * H_TO_H, matrix[1][i-1] * L_TO_H));
                matrix[1][i] = T_L * (Math.max(matrix[0][i-1] * H_TO_L, matrix[1][i-1] * L_TO_L));
            }
        }

        /**
         * FILE OUTPUT PART B
         */
        try{
            //Create writer
            BufferedWriter writer = new BufferedWriter(new FileWriter("4.o2"));

            // Start
            writer.write("- 0 ");

            // Write sequence
            for (int i = 0; i < sequence.length(); i++){
                writer.write(sequence.charAt(i) + " ");
            }

            // next line
            writer.write('\n');

            // Print viterbi matrix
            for (int i = 0; i < 3; i++){
                if (i == 0){
                    writer.write("0 ");
                } else if (i == 1){
                    writer.write("H ");
                } else if (i == 2){
                    writer.write("L ");
                }
                for (int j = 0; j < sequence.length() + 1; j++){
                    if (j == 0 || i == 0){
                        writer.write("0 ");
                    } else if (i == 1){
                        writer.write(Double.toString(matrix[0][j-1]) + " ");
                    } else if (i == 2){
                        writer.write(Double.toString(matrix[1][j-1]) + " ");
                    }
                }
                writer.write('\n');
            }

            // close writer
            writer.close();
        } catch (IOException e){
            System.out.println("Error writing to file 4.o2");
        }
        /**
         * END FILE OUTPUT PART B
         */

        // Find traceback
        String traceback = ViterbiTraceback(matrix);

        /**
         * FILE OUTPUT PART C
         */
        try{
            BufferedWriter writer = new BufferedWriter(new FileWriter("4.o3"));

            writer.write(traceback);

            writer.close();
        } catch (IOException e){
            System.out.println("Error writing to file 4.o3");
        }
        /**
         * END FILE OUTPUT PART C
         */

        /**
         * FILE OUTPUT PART D
         */
        try{
            BufferedWriter writer = new BufferedWriter(new FileWriter("4.o4"));

            if (traceback.charAt(sequence.length() - 1) == 'H'){
                writer.write(Double.toString(matrix[0][sequence.length() - 1]));
            } else {
                writer.write(Double.toString(matrix[1][sequence.length() - 1]));
            }

            writer.close();
        } catch (IOException e){
            System.out.println("Error writing to file 4.o4");
        }

        /**
         * END FILE OUTPUT PART D
         */
    }

    String ViterbiTraceback(double[][] matrix){
        // Initialize trace
        StringBuilder trace = new StringBuilder();

        // Default no multiple paths
        boolean multiplePaths = false;

        // Appends to the beginning
        // Traceback
        for (int i = matrix[0].length - 1; i >= 0; i--) {
            // If starting point
            if (i == matrix[0].length - 1){
                // If we're starting with H
                if (matrix[0][i] > matrix[1][i]){
                    // Append H
                    trace.append("H");

                // Not starting with H
                } else {
                    // Check if multiple
                    if(matrix[0][i] == matrix[1][i]){
                        // Set multiple to true
                        multiplePaths = true;
                    }

                    // Append L
                    trace.append("L");
                }
            }

            // If we are at H
            if (trace.charAt(0) == 'H') {
                // If path from H is higher prob
                if ((matrix[0][i] * H_TO_H) > (matrix[1][i] * L_TO_H)){
                    // Append H
                    trace.insert(0, "H");

                // Else path from L is higher prob (or equal)
                } else {
                    // Check for multiple
                    if ((matrix[0][i] * H_TO_H) == (matrix[1][i] * L_TO_H)){
                        // Set multiple to true
                        multiplePaths = true;
                    }

                    // Append L
                    trace.insert(0, "L");
                }
            // Else we are at L
            } else {
                // If path from L is higher prob
                if ((matrix[1][i] * L_TO_L) > (matrix[0][i] * H_TO_L)){
                    // Append L
                    trace.insert(0, "L");

                // Else path from H is higher prob
                } else {
                    // Check for multiple
                    if((matrix[1][i] * L_TO_L) == (matrix[0][i] * H_TO_L)){
                        // Set multiple to true
                        multiplePaths = true;
                    }

                    // Append H
                    trace.insert(0, "H");
                }
            }
        }

        /**
         * FILE OUTPUT PART E
         */

        try{
            BufferedWriter writer = new BufferedWriter(new FileWriter("4.o5"));

            if (multiplePaths){
                writer.write("YES");
            } else {
                writer.write("NO");
            }

            writer.close();
        } catch (IOException e){
            System.out.println("Error writing to file 4.o5");
        }

        /**
         * END FILE OUTPUT PART E
         */

        // Create string path to return from the trace
        String path = trace.toString();

        // Return path taken
        return path;
    }

    double ForwardAlgorithm(Sequence seq){
        double probability;

        // Split sequence
        String sequence = seq.getSequence();

        // Create matrix for scores
        double[][] matrix = new double[2][sequence.length()];

        // Fill out matrix[0][0] and matrix[1][0]
        if (sequence.charAt(0) == 'A'){
            matrix[0][0] = START_PROB * A_H;
        } else if (sequence.charAt(0) == 'C'){
            matrix[0][0] = START_PROB * C_H;
        } else if (sequence.charAt(0) == 'G'){
            matrix[0][0] = START_PROB * G_H;
        } else if (sequence.charAt(0) == 'T'){
            matrix[0][0] = START_PROB * T_H;
        }

        if (sequence.charAt(0) == 'A'){
            matrix[1][0] = START_PROB * A_L;
        } else if (sequence.charAt(0) == 'C'){
            matrix[1][0] = START_PROB * C_L;
        } else if (sequence.charAt(0) == 'G'){
            matrix[1][0] = START_PROB * G_L;
        } else if (sequence.charAt(0) == 'T'){
            matrix[1][0] = START_PROB * T_L;
        }

        // Fill out the rest of the table
        for (int i = 1; i < sequence.length(); i++){
            // Find values for next rows
            if (sequence.charAt(i) == 'A'){
                matrix[0][i] = A_H * Double.sum(matrix[0][i-1] * H_TO_H, matrix[1][i-1] * L_TO_H);
                matrix[1][i] = A_L * Double.sum(matrix[0][i-1] * H_TO_L, matrix[1][i-1] * L_TO_L);
            } else if (sequence.charAt(i) == 'C'){
                matrix[0][i] = C_H * Double.sum(matrix[0][i-1] * H_TO_H, matrix[1][i-1] * L_TO_H);
                matrix[1][i] = C_L * Double.sum(matrix[0][i-1] * H_TO_L, matrix[1][i-1] * L_TO_L);
            } else if (sequence.charAt(i) == 'G'){
                matrix[0][i] = G_H * Double.sum(matrix[0][i-1] * H_TO_H, matrix[1][i-1] * L_TO_H);
                matrix[1][i] = G_L * Double.sum(matrix[0][i-1] * H_TO_L, matrix[1][i-1] * L_TO_L);
            } else if (sequence.charAt(i) == 'T'){
                matrix[0][i] = T_H * Double.sum(matrix[0][i-1] * H_TO_H, matrix[1][i-1] * L_TO_H);
                matrix[1][i] = T_L * Double.sum(matrix[0][i-1] * H_TO_L, matrix[1][i-1] * L_TO_L);
            }
        }

        // Calculate probability
        probability = matrix[0][sequence.length() - 1] + matrix[1][sequence.length() - 1];

        // Make it log
        //probability = Math.log(probability)/Math.log(2.0);

        // run backwards algorithm
        BackwardAlgorithm(seq, matrix, probability);

        // Return probability
        return probability;
    }

    void BackwardAlgorithm(Sequence seq, double[][] forwardMatrix, double probability){
        // sequence
        String sequence = seq.getSequence();

        // initialize matrix
        double[][] matrix = new double[2][sequence.length()];

        // last column is 1
        matrix[0][sequence.length()-1] = 1;
        matrix[1][sequence.length()-1] = 1;

        // calculate the rest of the table
        for (int i = sequence.length() - 2; i >= 0; i--){
            // Find values for next rows
            if (sequence.charAt(i) == 'A'){
                matrix[0][i] = Double.sum(A_H * matrix[0][i+1] * H_TO_H, A_H * matrix[1][i+1] * L_TO_H);
                matrix[1][i] = Double.sum(A_L * matrix[0][i+1] * H_TO_L, A_L * matrix[1][i+1] * L_TO_L);
            } else if (sequence.charAt(i) == 'C'){
                matrix[0][i] = Double.sum(C_H * matrix[0][i+1] * H_TO_H, C_H * matrix[1][i+1] * L_TO_H);
                matrix[1][i] = Double.sum(C_L * matrix[0][i+1] * H_TO_L, C_L * matrix[1][i+1] * L_TO_L);
            } else if (sequence.charAt(i) == 'G'){
                matrix[0][i] = Double.sum(G_H * matrix[0][i+1] * H_TO_H, G_H * matrix[1][i+1] * L_TO_H);
                matrix[1][i] = Double.sum(G_L * matrix[0][i+1] * H_TO_L, G_L * matrix[1][i+1] * L_TO_L);
            } else if (sequence.charAt(i) == 'T'){
                matrix[0][i] = Double.sum(T_H * matrix[0][i+1] * H_TO_H, T_H * matrix[1][i+1] * L_TO_H);
                matrix[1][i] = Double.sum(T_L * matrix[0][i+1] * H_TO_L, T_L * matrix[1][i+1] * L_TO_L);
            }
        }

        /**
         * FILE OUTPUT PART F
         */
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter("4.o6"));

            // Calculate posterior probabilities and write to file
            writer.write(Double.toString((forwardMatrix[0][sequence.length() - 1] * matrix[0][sequence.length() - 1]) / probability));
            writer.write("\n");
            writer.write(Double.toString((forwardMatrix[1][sequence.length() - 1] * matrix[1][sequence.length() - 1]) / probability));

            // close writer
            writer.close();
        } catch (IOException e) {
            System.out.println("Error writing to file 4.o6");
        }

        /**
         * END FILE OUTPUT PART F
         */

    }

}
