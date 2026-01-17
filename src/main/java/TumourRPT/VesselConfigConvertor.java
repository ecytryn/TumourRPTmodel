package TumorRPT;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

/**
 * Loads vessel configuration from CSV files
 * 
 * Supports loading from:
 * - Filesystem paths (for local development)
 * - Classpath resources (for packaged jars)
 */
public class VesselConfigConvertor {
    public double[][] grid;

    /**
     * Initialize grid from CSV file
     * Tries resource loading first, falls back to filesystem
     * 
     * @param filePath Path to CSV file (e.g., "vasculature/uniform.csv")
     */
    public VesselConfigConvertor(String filePath) throws IOException {
        // Try loading as resource first (from src/main/resources)
        InputStream resourceStream = getClass().getClassLoader().getResourceAsStream(filePath);
        
        if (resourceStream != null) {
            // Load from classpath resource
            loadFromStream(resourceStream);
        } else {
            // Fall back to filesystem path
            loadFromFile(filePath);
        }
    }
    
    /**
     * Load grid from an input stream (for resources)
     */
    private void loadFromStream(InputStream stream) throws IOException {
        List<double[]> tempGrid = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new InputStreamReader(stream))) {
            String line;
            while ((line = br.readLine()) != null) {
                String[] values = line.split(",");
                double[] row = new double[values.length];
                for (int i = 0; i < values.length; i++) {
                    row[i] = Double.parseDouble(values[i].trim());
                }
                tempGrid.add(row);
            }
        }
        grid = tempGrid.toArray(new double[0][]);
    }
    
    /**
     * Load grid from filesystem (for local files)
     */
    private void loadFromFile(String filePath) throws IOException {
        List<double[]> tempGrid = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(filePath))) {
            String line;
            while ((line = br.readLine()) != null) {
                String[] values = line.split(",");
                double[] row = new double[values.length];
                for (int i = 0; i < values.length; i++) {
                    row[i] = Double.parseDouble(values[i].trim());
                }
                tempGrid.add(row);
            }
        }
        grid = tempGrid.toArray(new double[0][]);
    }

    /**
     * Print grid to console (for debugging)
     */
    public void printGrid() {
        for (double[] row : grid) {
            for (double val : row) {
                System.out.print(val + " ");
            }
            System.out.println();
        }
    }
    
    /**
     * Get grid dimensions
     */
    public int getRows() {
        return grid.length;
    }
    
    public int getCols() {
        return grid.length > 0 ? grid[0].length : 0;
    }
}
