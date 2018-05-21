import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import javax.swing.*;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

public class Paint {
    /**
     * Переменные - пути с файлами данных, использованных для построения графиков
     */
    private static String Filename_graph_x = System.getProperty("user.dir")+"\\data\\x.txt";
    private static String Filename_graph_y = System.getProperty("user.dir")+"\\data\\y.txt";
    private static String Filename_graph = System.getProperty("user.dir")+"\\data\\graph.jpeg";
    private static String Filename_graph_exp = System.getProperty("user.dir")+"\\data\\exp.jpeg";

    /**
     * Демонстрирует содержимое коллекции XYSeries
     * в промежутке от t0 до T
     * @param series коллекция для вывода
     * @param t0 левая граница интервала
     * @param T правая граница интервала
     */
    private static void show (XYSeries series, int t0, int T){
        for (int i = t0; i < T; i++){
            System.out.println("x = "+series.getX(i)+", y = "+series.getY(i));
        }
    }

    /**
     * Демонстрирует содержимое коллекции XYSeries
     * @param series коллекция для вывода
     */
    private static void show (XYSeries series){
        for (int i = 0; i < series.getItemCount(); i++){
            System.out.println("DataItem["+i+"] = "+series.getDataItem(i));
        }
    }

    /**
     * Выводит на фрейм график кривой exp(t)
     * @param t0 левая граница интервала
     * @param T правая граница интервала
     */
    private static void test_exp (double t0, double T) throws IOException {
        int width = 1280;
        int height = 720;

        XYSeries series = new XYSeries("exp(t)");

        for(double i = t0; i < T; i+=0.01){
            series.add(i, Math.exp(i));
        }

        XYDataset xyDataset = new XYSeriesCollection(series);
        JFreeChart chart = ChartFactory
                .createXYLineChart("y = exp(t)", "t", "y",
                        xyDataset,
                        PlotOrientation.VERTICAL,
                        true, true, true);
        JFrame frame =
                new JFrame("MinimalStaticChart");
        // Помещаем график на фрейм
        frame.getContentPane()
                .add(new ChartPanel(chart));
        frame.setSize(width,height);
        frame.setSize(width,height);
        ChartUtilities.saveChartAsJPEG( new File(Filename_graph_exp), chart , width , height );
        frame.show();
    }

    /**
     * Возращает коллекцию XYSeries,
     * сформированную из файлов
     */
    private static XYSeries read_file() throws IOException {
        XYSeries series = new XYSeries("σ(t)");
        BufferedReader reader_x = new BufferedReader(new FileReader(Filename_graph_x));
        BufferedReader reader_y = new BufferedReader(new FileReader(Filename_graph_y));
        String lineX, lineY;
        while (((lineX = reader_x.readLine()) != null)&&(lineY = reader_y.readLine()) != null) {
            series.add(Double.parseDouble(lineX), Double.parseDouble(lineY));
        }
        return series;
    }

    /**
     * Выводит на фрейм график кривой exp(t)
     * @param width ширина окна
     * @param height высота окна
     */
    public static void demonstration(int width, int height) throws IOException {
        //формирование XYSeries - коллекции
        XYSeries series = read_file();
        //формирование набора данных на основе ранее созданной коллекции
        XYDataset xyDataset = new XYSeriesCollection(series);
        JFreeChart chart = ChartFactory
                .createXYLineChart("y = σ(t)", "t", "σ",
                        xyDataset,
                        PlotOrientation.VERTICAL,
                        true, true, true);
        JFrame frame = new JFrame("MinimalStaticChart");
        // Помещаем график на фрейм
        frame.getContentPane().add(new ChartPanel(chart));
        frame.setSize(width,height);
        ChartUtilities.saveChartAsJPEG( new File(Filename_graph), chart , width , height );
        frame.show();
    }
}
