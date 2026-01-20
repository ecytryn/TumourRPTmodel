package TumorRPT;

import javax.imageio.IIOException;
import javax.imageio.ImageIO;
import javax.swing.*;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

public class DataLogger {

    public DataLogger() {
    }

    public void log(ArrayList<double[]> list, String fileName) {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(fileName, false))) {
            for (double[] element : list) {
                writer.write(convertToCSV(element));
                writer.newLine();
            }
            writer.flush();
        } catch (IOException e) {
            System.err.println("Error writing to log file: " + e.getMessage());
        }
    }

    /**
     * Save figure with timestamp and optional legend
     * @param fileName Output filename
     * @param drawer DaVinci object with visualization
     * @param dayCount Current day number for timestamp
     * @param showLegend If true, adds a color legend to the image
     */
    public void saveFigureTotal(String fileName, DaVinci drawer, int dayCount, boolean showLegend) {   
        // Works in both live and headless modes
        BufferedImage img = new BufferedImage(drawer.xDim, drawer.yDim, BufferedImage.TYPE_INT_RGB);
        
        // Copy pixels using DaVinci's GetPix method
        for (int x = 0; x < drawer.xDim; x++) {
            for (int y = 0; y < drawer.yDim; y++) {
                int color = drawer.GetPix(x, y);
                img.setRGB(x, y, color);
            }
        }
        
        // Add timestamp text overlay
        Graphics2D g2d = img.createGraphics();
        
        // Enable anti-aliasing for smoother text
        g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        g2d.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
        
        // Create timestamp text
        String timestamp = String.format("Day %d", dayCount);
        
        // Font settings - adjust size based on image dimensions
        int fontSize = Math.max(12, drawer.xDim / 20);  // Scale with image size
        Font font = new Font("Arial", Font.BOLD, fontSize);
        g2d.setFont(font);
        
        // Get text dimensions for positioning
        FontMetrics fm = g2d.getFontMetrics();
        int textWidth = fm.stringWidth(timestamp);
        int textHeight = fm.getHeight();
        
        // Position in top-left corner with padding
        int padding = 10;
        int x = padding;
        int y = padding + fm.getAscent();
        
        // Draw black outline for contrast (makes text readable on any background)
        g2d.setColor(Color.BLACK);
        for (int dx = -2; dx <= 2; dx++) {
            for (int dy = -2; dy <= 2; dy++) {
                if (dx != 0 || dy != 0) {
                    g2d.drawString(timestamp, x + dx, y + dy);
                }
            }
        }
        
        // Draw white text on top
        g2d.setColor(Color.WHITE);
        g2d.drawString(timestamp, x, y);
        
        // Add legend if requested
        if (showLegend) {
            drawLegend(g2d, img.getWidth(), img.getHeight());
        }
        
        // Always add scale bar
        drawScaleBar(g2d, img.getWidth(), img.getHeight());
        
        g2d.dispose();
        
        // Save image
        try{
            ImageIO.write(img,"png",new File(fileName));
        }catch (IIOException e){
            e.printStackTrace();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        System.out.println("Image Saved with timestamp: " + fileName);
    }
    
    /**
     * Save figure with timestamp but no legend
     * @param fileName Output filename
     * @param drawer DaVinci object with visualization  
     * @param dayCount Current day number for timestamp
     */
    public void saveFigureTotal(String fileName, DaVinci drawer, int dayCount) {
        saveFigureTotal(fileName, drawer, dayCount, false);
    }
    
    /**
     * Draw legend showing cell types and oxygen colors
     */
    private void drawLegend(Graphics2D g2d, int imgWidth, int imgHeight) {
        // Legend configuration
        int legendWidth = 180;
        int legendHeight = 200;
        int padding = 10;
        int x = imgWidth - legendWidth - padding;  // Top-right corner
        int y = padding;
        
        // Semi-transparent background
        g2d.setColor(new Color(0, 0, 0, 200));  // Black with alpha
        g2d.fillRoundRect(x, y, legendWidth, legendHeight, 10, 10);
        
        // White border
        g2d.setColor(Color.WHITE);
        g2d.drawRoundRect(x, y, legendWidth, legendHeight, 10, 10);
        
        // Legend entries
        String[] labels = {
            "Oxygen (low→high)",
            "Normoxic tumor",
            "Hypoxic tumor", 
            "Necrotic",
            "Apoptotic",
            "Vessel"
        };
        
        // Colors from SimParams.COLORLIST
        // We'll create representative colors for each
        Color[] colors = new Color[6];
        colors[0] = new Color(0, 0, 128);      // Dark blue for oxygen (representative)
        colors[1] = new Color(SimParams.COLORLIST[SimParams.NORMAL]);   // Pink/magenta
        colors[2] = new Color(SimParams.COLORLIST[SimParams.HYPOXIC]);     // Purple
        colors[3] = new Color(SimParams.COLORLIST[SimParams.NECROTIC]);    // Gray
        colors[4] = new Color(SimParams.COLORLIST[SimParams.APOPTOTIC]);     // Yellow
        colors[5] = new Color(SimParams.COLORLIST[SimParams.VESSEL]);         // Cyan
        
        // Draw legend items
        int itemHeight = 25;
        int colorBoxSize = 15;
        int textOffset = 25;
        Font legendFont = new Font("Arial", Font.PLAIN, 12);
        g2d.setFont(legendFont);
        
        for (int i = 0; i < labels.length; i++) {
            int itemY = y + 20 + i * itemHeight;
            
            // Draw color box
            g2d.setColor(colors[i]);
            g2d.fillRect(x + 10, itemY, colorBoxSize, colorBoxSize);
            
            // Draw white border around color box
            g2d.setColor(Color.WHITE);
            g2d.drawRect(x + 10, itemY, colorBoxSize, colorBoxSize);
            
            // Draw label text
            g2d.setColor(Color.WHITE);
            g2d.drawString(labels[i], x + 10 + textOffset, itemY + colorBoxSize - 2);
        }
    }
    
    /**
     * Original method without timestamp (kept for backward compatibility)
     */
    public void saveFigureTotal(String fileName, DaVinci drawer) {
        // Works in both live and headless modes
        BufferedImage img = new BufferedImage(drawer.xDim, drawer.yDim, BufferedImage.TYPE_INT_RGB);
        for (int x = 0; x < drawer.xDim; x++) {
            for (int y = 0; y < drawer.yDim; y++) {
                int color = drawer.GetPix(x, y);
                img.setRGB(x, y, color);
            }
        }
        try{
            ImageIO.write(img,"png",new File(fileName));
        }catch (IIOException e){
            e.printStackTrace();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        System.out.println("Image Saved ...");
    }

    public void saveVessel(String fileName, DaVinci drawer) {
        // Works in both live and headless modes
        BufferedImage img = new BufferedImage(drawer.xDim, drawer.yDim, BufferedImage.TYPE_INT_RGB);
        for (int x = 0; x < drawer.xDim; x++) {
            for (int y = 0; y < drawer.yDim; y++) {
                int color = drawer.GetPix(x, y);
                img.setRGB(x, y, color);
            }
        }
        try{
            ImageIO.write(img,"png",new File(fileName));
        }catch (IIOException e){
            e.printStackTrace();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        System.out.println("Image Saved ...");
    }

    private String convertToCSV(double[] element) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < element.length; i++) {
            sb.append(element[i]);
            if (i < element.length - 1) {
                sb.append(",");
            }
        }
        return sb.toString();
    }

    private void drawScaleBar(Graphics2D g2d, int imgWidth, int imgHeight) {
        // Scale bar parameters
        int barLength_cells = 50;  // 50 cells = 500 um
        double barLength_um = barLength_cells * SimParams.CELL_LENGTH * 1e6;
        
        int barHeight = 8;
        int padding = 15;
        int textGap = 5;
        
        // Position in bottom-right
        int x = imgWidth - barLength_cells - padding;
        int y = imgHeight - padding - barHeight - 20;
        
        // Draw white bar with black outline
        g2d.setColor(Color.WHITE);
        g2d.fillRect(x, y, barLength_cells, barHeight);
        g2d.setColor(Color.BLACK);
        g2d.drawRect(x, y, barLength_cells, barHeight);
        
        // Draw scale bar label
        String label = String.format("%.0f μm", barLength_um);
        Font scaleFont = new Font("Arial", Font.BOLD, 12);
        g2d.setFont(scaleFont);
        FontMetrics fm = g2d.getFontMetrics();
        int textWidth = fm.stringWidth(label);
        
        int textX = x + (barLength_cells - textWidth) / 2;
        int textY = y + barHeight + textGap + fm.getAscent();
        
        // Draw text with black outline for visibility
        g2d.setColor(Color.BLACK);
        for (int dx = -1; dx <= 1; dx++) {
            for (int dy = -1; dy <= 1; dy++) {
                if (dx != 0 || dy != 0) {
                    g2d.drawString(label, textX + dx, textY + dy);
                }
            }
        }
        g2d.setColor(Color.WHITE);
        g2d.drawString(label, textX, textY);
    }

/*
	private void drawScaleBar(Graphics2D g2d, int width, int height) {
		// Choose scale bar size in pixels
		int barLength = 10;  // adjust to match physical scale
		int barHeight = 3;
	
		// Placement: bottom-left corner with some padding
		int padding = 20;
		int x = padding;
		int y = height - padding - barHeight;
	
		// Draw background (optional shadow)
		g2d.setColor(new Color(0,0,0,150)); // translucent black
		g2d.fillRect(x - 2, y - 2, barLength + 4, barHeight + 4);
	
		// Draw white bar
		g2d.setColor(Color.WHITE);
		g2d.fillRect(x, y, barLength, barHeight);
	
		// Draw label (e.g., "100 μm")
		g2d.setFont(new Font("SansSerif", Font.PLAIN, 16));
		g2d.drawString("Scale: 100 μm", x, y - 5);
	}
*/

}
