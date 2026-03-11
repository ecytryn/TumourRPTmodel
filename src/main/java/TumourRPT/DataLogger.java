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
import java.io.PrintWriter;

public class DataLogger {

    public DataLogger() {
    }

	public void log(ArrayList<double[]> list, String fileName) {
		try (BufferedWriter writer = new BufferedWriter(new FileWriter(fileName, false))) {
			// Write header based on file type
			if (fileName.contains("populations.csv")) {
				writer.write("healthy,normal,hypoxic,necrotic,apoptotic,vessel");
				writer.newLine();
			} else if (fileName.contains("pkStateVariables.csv")) {
				writer.write("N_cen_hot,N_cen_cold,N_v_hot,N_v_cold,N_ec_hot,N_ec_cold,N_b_hot,N_b_cold,N_ic_hot,N_ic_cold,A_blob");
				writer.newLine();
			} else if (fileName.contains("dose.csv")) {
				writer.write("dose_rate_Gy_per_hr");
				writer.newLine();
			}
			
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
    public void saveFigureTotal(String fileName, DaVinci drawer, int dayCount, boolean showLegend, boolean showScaleBar) {   
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
        if ( MyUtils.isElementPresent(SimParams.INJECTION_SCHEDULE, dayCount) ) {
        	timestamp = timestamp + " (injection)";
        }
        
        // Font settings - adjust size based on image dimensions
        int fontSize = (int)(20 * SimParams.FONT_SCALE_FACTOR);  // Scale with image size
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
        if (showScaleBar) {
			drawScaleBar(g2d, img.getWidth(), img.getHeight());
		}        
        g2d.dispose();
        
        // Save image
        try{
            ImageIO.write(img,"png",new File(fileName));
        }catch (IIOException e){
            e.printStackTrace();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
//        System.out.println("Image Saved with timestamp: " + fileName);
    }
        
    /**
     * Save figure with timestamp but no legend
     * @param fileName Output filename
     * @param drawer DaVinci object with visualization  
     * @param dayCount Current day number for timestamp
     */
    public void saveFigureTotal(String fileName, DaVinci drawer, int dayCount) {
        saveFigureTotal(fileName, drawer, dayCount, false, false);
    }

	/**
	 * Save a zoomed, upscaled crop of the tumour centred on the grid centre.
	 * Zoom parameters are controlled by SimParams.ZOOM_HALF_WIDTH and
	 * SimParams.ZOOM_PIXEL_SCALE. Draws day label and a scale bar that
	 * reflects the physical size of the zoomed region.
	 *
	 * @param fileName  Output PNG path
	 * @param drawer    DaVinci with rendered pixel buffer
	 * @param grid      Grid (for centre coordinates)
	 * @param dayCount  Current simulation day (for label)
	 */
	public void saveFigureZoomed(String fileName, DaVinci drawer, Grid grid, 
						int dayCount, boolean showLegend, boolean showScaleBar) {
	
		// ---------------------------------------------------------------
		// Zoom geometry — derived from SimParams tuning constants
		// ---------------------------------------------------------------
		final int hw   = SimParams.ZOOM_HALF_WIDTH;   // half-width in cells
		final int scale = SimParams.ZOOM_PIXEL_SCALE; // pixels per cell
	
		final int cx = grid.xDim / 2;  // grid centre x (cell coords)
		final int cy = grid.yDim / 2;  // grid centre y
	
		final int outPx = (2 * hw) * scale;  // output image side in pixels
	
		BufferedImage img = new BufferedImage(outPx, outPx, BufferedImage.TYPE_INT_RGB);
	
		// ---------------------------------------------------------------
		// Copy + upscale: each source cell -> scale x scale pixel block
		// ---------------------------------------------------------------
		for (int dy = -hw; dy < hw; dy++) {
			for (int dx = -hw; dx < hw; dx++) {
				int srcX = cx + dx;
				int srcY = cy + dy;
	
				// Clamp to grid bounds (shouldn't trigger for reasonable hw)
				srcX = Math.max(0, Math.min(grid.xDim - 1, srcX));
				srcY = Math.max(0, Math.min(grid.yDim - 1, srcY));
	
				int color = drawer.GetPix(srcX, srcY);
	
				// Write to all pixels in the scaled block
				int destX0 = (dx + hw) * scale;
				int destY0 = (dy + hw) * scale;
				for (int py = 0; py < scale; py++) {
					for (int px = 0; px < scale; px++) {
						img.setRGB(destX0 + px, destY0 + py, color);
					}
				}
			}
		}
	
		// ---------------------------------------------------------------
		// Overlays: day label + scale bar
		// ---------------------------------------------------------------
		Graphics2D g2d = img.createGraphics();
		g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING,  RenderingHints.VALUE_ANTIALIAS_ON);
		g2d.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
	
		// --- Day label (top-left) ---
		String timestamp = String.format("Day %d", dayCount);
		if (MyUtils.isElementPresent(SimParams.INJECTION_SCHEDULE, dayCount)) {
			timestamp += " (injection)";
		}
		int fontSize = SimParams.ZOOM_FONT_SIZE;
		Font font = new Font("Arial", Font.BOLD, fontSize);
		g2d.setFont(font);
		FontMetrics fm = g2d.getFontMetrics();
		int labelPad = 8;
		int labelX = labelPad;
		int labelY = labelPad + fm.getAscent();
		// Black outline
		g2d.setColor(Color.BLACK);
		for (int ddx = -2; ddx <= 2; ddx++)
			for (int ddy = -2; ddy <= 2; ddy++)
				if (ddx != 0 || ddy != 0)
					g2d.drawString(timestamp, labelX + ddx, labelY + ddy);
		g2d.setColor(Color.WHITE);
		g2d.drawString(timestamp, labelX, labelY);
	
		if (showScaleBar) {
			// --- Scale bar (bottom-right) ---
			// Choose a round physical length that fits comfortably in the image.
			// We pick the largest of {50, 100, 200, 500} um that is <= 40% of the crop width.
			double cropWidth_um = 2.0 * hw * SimParams.CELL_LENGTH * 1e6;  // um
			int[] niceValues_um = {50, 100, 200, 500};
			int barPhys_um = niceValues_um[0];
			for (int v : niceValues_um) {
				if (v <= 0.40 * cropWidth_um) barPhys_um = v;
			}
			// Convert physical length to output pixels
			double umPerCell = SimParams.CELL_LENGTH * 1e6;
			int barPx = (int) Math.round((barPhys_um / umPerCell) * scale);
		
			int barH   = SimParams.ZOOM_SCALE_BAR_HEIGHT;
			int sbPad  = 12;
			int barX   = outPx - barPx - sbPad - 20;
			int barY   = outPx - barH  - sbPad - fontSize - 6;
		
			g2d.setColor(Color.WHITE);
			g2d.fillRect(barX, barY, barPx, barH);
			g2d.setColor(Color.BLACK);
			g2d.drawRect(barX, barY, barPx, barH);
		
			// Scale bar label
			String barLabel = barPhys_um >= 1000
				? String.format("%.0f mm", barPhys_um / 1000.0)
				: String.format("%d \u00b5m", barPhys_um);   // µm using unicode
			Font scaleFont = new Font("Arial", Font.PLAIN, Math.max(9, fontSize - 2));
			g2d.setFont(scaleFont);
			FontMetrics sfm = g2d.getFontMetrics();
			int barLabelX = barX + (barPx - sfm.stringWidth(barLabel)) / 2;
			int barLabelY = barY + barH + 4 + sfm.getAscent();
			g2d.setColor(Color.BLACK);
			for (int ddx = -1; ddx <= 1; ddx++)
				for (int ddy = -1; ddy <= 1; ddy++)
					if (ddx != 0 || ddy != 0)
						g2d.drawString(barLabel, barLabelX + ddx, barLabelY + ddy);
			g2d.setColor(Color.WHITE);
			g2d.drawString(barLabel, barLabelX, barLabelY);
		}	
		g2d.dispose();

		if (showLegend) {
			drawLegend(g2d, outPx, outPx);
		}	

		// ---------------------------------------------------------------
		// Save
		// ---------------------------------------------------------------
		try {
			ImageIO.write(img, "png", new File(fileName));
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

    
    /**
     * Draw legend showing cell types and oxygen colors
     */
    private void drawLegend(Graphics2D g2d, int imgWidth, int imgHeight) {
        // Legend configuration
        int legendWidth = 350;
        int legendHeight = 290;
        int padding = 25;
        int x = imgWidth - legendWidth - padding;  // Top-right corner
        int y = 90;
        
        // Semi-transparent background
        g2d.setColor(new Color(0, 0, 0, 100));  // Black with alpha
        g2d.fillRoundRect(x, y, legendWidth, legendHeight, 10, 10);
        
        // White border
        g2d.setColor(Color.WHITE);
        g2d.drawRoundRect(x, y, legendWidth, legendHeight, 10, 10);
        
        // Legend entries
        String[] labels = {
            "Normoxic cell",
            "Hypoxic cell", 
            "Necrotic cell",
            "Apoptotic cell",
            "Oxygen (graduated)",
            "Open vessel",
            "Occluded vessel"
        };
        
        // Colors from SimParams.COLORLIST
        // We'll create representative colors for each
        Color[] colors = new Color[7];
        colors[0] = new Color(SimParams.COLORLIST[SimParams.NORMAL]);   // light green
        colors[1] = new Color(SimParams.COLORLIST[SimParams.HYPOXIC]);     // mid green
        colors[2] = new Color(SimParams.COLORLIST[SimParams.NECROTIC]);    // dark green
        colors[3] = new Color(SimParams.COLORLIST[SimParams.APOPTOTIC]);     // magenta
        colors[4] = new Color(24, 21, 248);      //  blue for oxygen (representative)
        colors[5] = new Color(SimParams.COLORLIST[SimParams.VESSEL]);         // bright red
        colors[6] = new Color(SimParams.COLORLIST[6]);         // dark maroon for occluded vessel
                
        // Draw legend items
        int itemHeight = 40;
        int colorBoxSize = 25;
        int textOffset = 35;
		int columnWidth = 225;
		int fontSize = (int)(20 * SimParams.FONT_SCALE_FACTOR);
		Font legendFont = new Font("Arial", Font.PLAIN, fontSize);
        g2d.setFont(legendFont);
        
        for (int i = 0; i < labels.length; i++) {
//			int column = i / 4;  // 0 or 1 (4 items per column)
//			int row = i % 4;     // 0-3

//			int itemX = x + 10 + column * columnWidth;
			int itemX = x + 10;
//			int itemY = y + 12 + row * itemHeight;
			int itemY = y + 12 + i * itemHeight;

			// Draw color box
			g2d.setColor(colors[i]);
			g2d.fillRect(itemX, itemY, colorBoxSize, colorBoxSize);
			
			// Draw white border around color box
			g2d.setColor(Color.WHITE);
			g2d.drawRect(itemX, itemY, colorBoxSize, colorBoxSize);

            // Draw label text
            g2d.setColor(Color.WHITE);
            g2d.drawString(labels[i], itemX + textOffset, itemY + colorBoxSize + 1);
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
//        System.out.println("Image Saved ...");
    }
    
	/**
	 * Save survival fraction visualization as heatmap
	 * @param fileName Output filename
	 * @param grid Grid containing cells
	 * @param dayCount Current day number for timestamp
	 */
	public void saveSFVisualization(String fileName, Grid grid, int dayCount) {
		int xDim = grid.xDim;
		int yDim = grid.yDim;
		
		// Create buffered image with extra space for color bar
		int imageWidth = xDim + 120;
		int imageHeight = yDim;
		java.awt.image.BufferedImage img = new java.awt.image.BufferedImage(
			imageWidth, imageHeight, java.awt.image.BufferedImage.TYPE_INT_RGB);
		
		// White background
		java.awt.Graphics2D g2d_bg = img.createGraphics();
		g2d_bg.setColor(java.awt.Color.WHITE);
		g2d_bg.fillRect(0, 0, imageWidth, imageHeight);
		g2d_bg.dispose();

		// Draw main heatmap
		for (int x = 0; x < xDim; x++) {
			for (int y = 0; y < yDim; y++) {
				int idx = grid.I(x, y);
				Cell cell = grid.GetAgent(idx);
				
				int color;
				if (cell == null) {
					color = 0xFFFFFFFF;  // White for empty
				} else if (cell.type == SimParams.VESSEL) {
					color = 0xFFFF0000;  // Red for vessels
				} else {
					// Show SF for ALL living cells, mark dead/problem cells specially
					double sf = cell.survivalProb;
					
					if (cell.type == SimParams.NECROTIC) {
						color = 0xFF404040;  // Dark gray for necrotic
					} else if (cell.type == SimParams.APOPTOTIC) {
						color = 0xFF800080;  // Purple for apoptotic
					} else if (sf < 0.0 || Double.isNaN(sf)) {
						color = 0xFFFF00FF;  // Magenta for invalid SF
					} else if (sf == 0.0) {
						color = 0xFFFFFF00;  // Yellow for SF exactly 0.0 (suspicious!)
					} else {
						// Zoom into high SF range to see variation
						double sf_rescaled = 0.8 * sf; 
						sf_rescaled = Math.max(0.0, Math.min(1.0, sf_rescaled));
						color = HAL.Util.HeatMapBGR(sf_rescaled, 0, 1);
					}
				}
				
				img.setRGB(x, y, color);
			}
		}		
		
		// Draw color bar
		int barX = xDim + 20;
		int barWidth = 30;
		int barY = 50;
		int barHeight = yDim - 100;
		
		for (int i = 0; i < barHeight; i++) {
			double sf = 1.0 - (double)i / barHeight;  // Top = 1.0, bottom = 0.0
			int color = HAL.Util.HeatMapBGR(sf, 0, 1);
			for (int j = 0; j < barWidth; j++) {
				img.setRGB(barX + j, barY + i, color);
			}
		}
		
		// Add text annotations
		java.awt.Graphics2D g2d = img.createGraphics();
		g2d.setRenderingHint(java.awt.RenderingHints.KEY_ANTIALIASING, 
							java.awt.RenderingHints.VALUE_ANTIALIAS_ON);
		g2d.setRenderingHint(java.awt.RenderingHints.KEY_TEXT_ANTIALIASING,
							java.awt.RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
		
		// Title with timestamp
		String title = String.format("Day %d - Survival Fraction", dayCount);
		java.awt.Font titleFont = new java.awt.Font("Arial", java.awt.Font.BOLD, 18);
		g2d.setFont(titleFont);
		
		int titleX = 15;
		int titleY = 25;
		
		// Black outline for contrast
		g2d.setColor(java.awt.Color.BLACK);
		for (int dx = -2; dx <= 2; dx++) {
			for (int dy = -2; dy <= 2; dy++) {
				if (dx != 0 || dy != 0) {
					g2d.drawString(title, titleX + dx, titleY + dy);
				}
			}
		}
		g2d.setColor(java.awt.Color.WHITE);
		g2d.drawString(title, titleX, titleY);
		
		// Color bar border and labels
		g2d.setColor(java.awt.Color.BLACK);
		g2d.drawRect(barX, barY, barWidth, barHeight);
		
		g2d.setFont(new java.awt.Font("Arial", java.awt.Font.PLAIN, 11));
		
		// Max label (SF = 1.0)
		g2d.drawString("1.0", barX + barWidth + 5, barY + 5);
		g2d.setFont(new java.awt.Font("Arial", java.awt.Font.PLAIN, 9));
		g2d.drawString("(survive)", barX + barWidth + 5, barY + 17);
		
		// Mid label (SF = 0.5)
		g2d.setFont(new java.awt.Font("Arial", java.awt.Font.PLAIN, 11));
		g2d.drawString("0.5", barX + barWidth + 5, barY + barHeight / 2);
		
		// Min label (SF = 0.0)
		g2d.drawString("0.0", barX + barWidth + 5, barY + barHeight - 5);
		g2d.setFont(new java.awt.Font("Arial", java.awt.Font.PLAIN, 9));
		g2d.drawString("(die)", barX + barWidth + 5, barY + barHeight + 7);
		
		// Legend
		g2d.setFont(new java.awt.Font("Arial", java.awt.Font.PLAIN, 10));
		int refY = barY + barHeight + 30;
		g2d.drawString("Colors:", barX - 10, refY);
		g2d.drawString("Blue = high SF", barX - 10, refY + 15);
		g2d.drawString("Red = low SF", barX - 10, refY + 30);
		g2d.drawString("Gray = dead", barX - 10, refY + 45);
		
		g2d.dispose();
		
		// Save image
		try {
			javax.imageio.ImageIO.write(img, "png", new java.io.File(fileName));
//			System.out.printf("Saved SF visualization: %s\n", fileName);
		} catch (Exception e) {
			System.err.printf("Warning: Could not save SF image: %s\n", e.getMessage());
		}
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
//        System.out.println("Image Saved ...");
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
        int barLength_cells = 100;  // 100 cells = 1 mm
        double barLength_um = barLength_cells * SimParams.CELL_LENGTH * 1e6;
        
        int barHeight = 8;
        int padding = 25;
        int textGap = 5;
        
        // Position in top-right
        int x = imgWidth - barLength_cells - padding;
        int y = imgWidth - barHeight - padding - 30;
        
        // Draw white bar with black outline
        g2d.setColor(Color.WHITE);
        g2d.fillRect(x, y, barLength_cells, barHeight);
        g2d.setColor(Color.BLACK);
        g2d.drawRect(x, y, barLength_cells, barHeight);
        
        // Draw scale bar label
        String label = String.format("%.0f mm", barLength_um*1e-3);
		int fontSize = (int)(20 * SimParams.FONT_SCALE_FACTOR);
        Font scaleFont = new Font("Arial", Font.PLAIN, fontSize);
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
	
		// Draw label (e.g., "100 um")
		g2d.setFont(new Font("SansSerif", Font.PLAIN, 16));
		g2d.drawString("Scale: 100 um", x, y - 5);
	}
*/

	public void saveOxygenFieldCSV(String fileName, Grid grid) {
		try (PrintWriter out = new PrintWriter(new FileWriter(fileName))) {
			for (int x = 0; x < grid.xDim; x++) {
				for (int y = 0; y < grid.yDim; y++) {
					int idx = grid.I(x, y);
					double o2 = grid.oxygenGrid.Get(idx);
					out.print(o2);
					if (y < grid.yDim - 1) out.print(",");
				}
				out.println();
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

}
