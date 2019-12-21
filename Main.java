import java.io.*;

public class Main {

    public static void main(String[] args) {
        // Create file
        File file = new File("4.in");

        String name = "";
        String seq = "";

        // Try reading from file
        try{
            // Scanner to read form file
            BufferedReader reader = new BufferedReader(new FileReader(file));

            // Set values
            name = reader.readLine().substring(1);
            seq = reader.readLine();

            // close reader
            reader.close();
        } catch (IOException e) {
            // Tell user file is not found
            System.out.println("Error file not found.");
        }

        // Create sequence
        Sequence sequence = new Sequence(name, seq);

        // Create out model
        HMM model = new HMM();

        // Get probability of producing X given the model
        double probProdX = model.ForwardAlgorithm(sequence);

        /**
         * FILE OUTPUT PART A
         */

        // Try writing to file
        try{
            // Scanner to read form file
           BufferedWriter writer = new BufferedWriter(new FileWriter("4.o1"));

            // Write probability of sequence x given model
            writer.write(Double.toString(probProdX));

            // close writer
            writer.close();
        } catch (IOException e) {
            System.out.println("Error writing to file 4.o1");
        }

        /**
         * END FILE OUTPUT PART A
         */

        // Run viterbi algorithm
        model.ViterbiAlgorithm(sequence);

    }
}
