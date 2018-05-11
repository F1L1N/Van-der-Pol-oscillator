import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Map;
import java.util.Scanner;
import java.util.TreeMap;

class SigmaFunction {
    //амплитуда
    private double A;
    //угловая частота
    private double omega;
    //шаг выполнения оценки
    private double h;
    //управляющий параметр
    private double nu;
    //левая граница рабочего промежутка
    private double t0;
    //правая граница рабочего промежутка
    private double T;
    //переменные для хранения максимальных и минимальных значений нормы сигмы - матрицы,
    //а также их индексов
    private double g_max = -Math.pow(10,6);
    private double g_max_t;
    private double g_min = Math.pow(10,6);
    private double g_min_t;
    //пути создания файлов с координатами
    private static String Filename_text = System.getProperty("user.dir")+"\\data\\evaluation.txt";
    private static String Filename_graph_X = System.getProperty("user.dir")+"\\data\\x.txt";
    private static String Filename_graph_Y = System.getProperty("user.dir")+"\\data\\y.txt";
    //словарь для хранения пар "X - Y координаты"
    private static TreeMap<Double, Double> dictionary = new TreeMap<>();

    //констуркторы класса
    SigmaFunction (double s1, double s2, double s3, double s4, double s5, double s6){
        h = s1;
        A = s2;
        omega = s3;
        nu = s4;
        t0 = s5;
        T = s6;
    }

    SigmaFunction (){
        System.out.print("Ввод шага: ");
        h = new Scanner(System.in).nextDouble();
        System.out.print("Ввод амплитуды: ");
        A = new Scanner(System.in).nextDouble();
        System.out.print("Ввод угловой частоты: ");
        omega = new Scanner(System.in).nextDouble();
        System.out.print("Ввод управляющего параметра: ");
        nu = new Scanner(System.in).nextDouble();
        System.out.print("Ввод левой границы интервала: ");
        t0 = new Scanner(System.in).nextDouble();
        System.out.print("Ввод правой границы интервала: ");
        T = new Scanner(System.in).nextDouble();
    }

    //метода для получения nu и последующего подсчета
    private double[] nuAnalysis(){
        double[] lambda = new double[2];
        if ((nu >= 0)&&(nu < 2)) {
            lambda[0] = (nu + Math.sqrt(Math.abs(Math.pow(nu,2)-4))*Math.sqrt(-1))/2;
            lambda[1] = (nu - Math.sqrt(Math.abs(Math.pow(nu,2)-4))*Math.sqrt(-1))/2;
        }else
            if (nu >= 2){
                lambda[0] = (nu + Math.sqrt(Math.pow(nu,2)-4))/2;
                lambda[1] = (nu - Math.sqrt(Math.pow(nu,2)-4))/2;
            }else{
                System.out.print("nu из неверного промежутка");
            }
        System.out.println(Arrays.toString(lambda));
        return lambda;
    }

    //функция формирования фундаментальной матрицы
    private double[][] F (double t){
        double[][] F = new double[2][2];
        F[0][0] = (1-t)*Math.exp(t);
        F[0][1] = t*Math.exp(t);
        F[1][0] = -t*Math.exp(t);
        F[1][1] = (1+t)*Math.exp(t);
        return F;
    }

    //функция формирования матрицы, обратной фундаментальной
    private double[][] F_reverse (double t){
        double[][] F_reverse = new double[2][2];
        F_reverse[0][0] = (1+t)/Math.exp(t);
        F_reverse[0][1] = -t/Math.exp(t);
        F_reverse[1][0] = t/Math.exp(t);
        F_reverse[1][1] = (1-t)/Math.exp(t);
        return F_reverse;
    }

    //функция формирования матрицы Коши
    private double[][] K (double t, double tao){
        return matrix_mult(F(t), F_reverse(tao));
    }

    //набор подынтегральных выражений, используемыхв расчетах
    private double f(double tao, double t, String s){
        switch (s){
            case "a1":
                return Math.exp(tao)*Math.sin(omega*tao)*((1-tao)*(1+t)-t*tao)/Math.exp(t);
            case "a2":
                return Math.exp(tao)*Math.sin(omega*tao)*(-t*(1-tao)+tao*(1-t))/Math.exp(t);
            case "a3":
                return Math.exp(tao)*Math.sin(omega*tao)*(tao*(1+t)-t*(1+tao))/Math.exp(t);
            case "a4":
                return Math.exp(tao)*Math.sin(omega*tao)*(-t*tao+(1+tao)*(1-t))/Math.exp(t);
            case "b1":
                return (1+tao)*Math.sin(omega*tao)/Math.exp(tao);
            case "b2":
                return -tao*Math.sin(omega*tao)/Math.exp(tao);
            case "b3":
                return (1-tao)*Math.sin(omega*tao)/Math.exp(tao);
            default: return 0.0;
        }
    }

    /*
    Метод Симпсона для нахождения определенного интеграла
    */
    private double integral(String s, double a, double b){
        int n = 1000;
        double sum = 0, sum2 = 0;
        double[] x = new double[n];
        double h = (b-a)/n;
        for(int i = 0; i < n; i++){
            x[i] = a+i*h;
        }
        for(int i = 1; i < n; i++){
            sum += f(x[i], b, s);
            sum2 += f((x[i-1]+x[i])/2, b, s);
        }
        return h/6*(f(a, b, s)+f(b, b, s)+2*sum+4*(sum2+b));
    }


    //метод для суммирования матриц
    private double[][] matrix_sum(double[][] a, double[][] b){
        double[][] result = new double[a.length][a[0].length];
        for (int i = 0; i < a.length; i++)
            for (int j = 0; j < a[i].length; j++)
                result[i][j] = a[i][j] + b[i][j];
        return result;
    }

    //метод для получения суммы матрицы с числом
    private double[][] matrix_sum(double[][] a, double b){
        for (int i = 0; i < a.length; i++)
            for (int j = 0; j < a[i].length; j++)
                a[i][j] = a[i][j] + b;
        return a;
    }

    //метод для получения евклидовой нормы матрицы
    private static double euc_norm(double[][] sigma){
        double sum = 0;
        for (double[] h : sigma)
            for (int j = 0; j < sigma[0].length;j++) {
                sum += Math.pow(h[j], 2);
            }
        return Math.sqrt(sum);
    }

    //метод для получения произведения матриц
    private double[][] matrix_mult(double[][] mA, double[][] mB){
        double[][] res = new double[mA.length][mB[0].length];
        for (int i = 0; i < mA.length; i++) {
            for (int j = 0; j < mB[0].length; j++) {
                for (int k = 0; k < mB.length; k++) {
                    res[i][j] += mA[i][k] * mB[k][j];
                }
            }
        }
        return res;
    }

    //метод получения произведения матрицы на число
    private double[][] matrix_mult(double[][] mA, double mB){
        for (int i = 0; i < mA.length; i++) {
            for (int j = 0; j < mA[i].length; j++) {
                mA[i][j] = mA[i][j] * mB;
            }
        }
        return mA;
    }

    //метод для формирования файлов с координатами для графиков
    private void file_graphic(){
        //сформируем первый файл с X - координатами
        try(FileWriter writer = new FileWriter(Filename_graph_X, false))
        {
            System.out.println("Производится запись файла по пути "+Filename_graph_X+"...");
            for (Map.Entry entry : dictionary.entrySet()) {
                double temp = (double) entry.getKey();
                writer.write(temp+"\r\n");
            }
            writer.flush();
            System.out.println("Запись завершена.");
        }
        catch(IOException ex){
            System.out.println(ex.getMessage());
        }
        //сформируем первый файл с Y - координатами
        try(FileWriter writer = new FileWriter(Filename_graph_Y, false))
        {
            System.out.println("Производится запись файла по пути "+Filename_graph_Y+"...");
            for (Map.Entry entry : dictionary.entrySet()) {
                double temp = (double) entry.getValue();
                writer.write(temp+"\r\n");
            }
            writer.flush();
            System.out.println("Запись завершена.");
        }
        catch(IOException ex){
            System.out.println(ex.getMessage());
        }
    }

    //метод для получения сигма - функции
    private double[][] form_sigma(double t0, double T){
        double[][] sigma = new double[2][2];
        int j, k, p = 0;
        //сформируем первое слагаемое из формулы
        double[][] K1 = new double[2][2];
        for (j = 0; j < sigma.length; j++)
            for (k = 0; k < sigma[0].length; k++)
            {
                p++;
                K1[j][k] = A*integral("a"+p, t0, T);
            }
        //сформируем второе слагаемое
        double[][] K2 = new double[2][2];
        K2[0][0] = integral("b1", t0, T);
        K2[0][1] = integral("b2", t0, T);
        K2[1][0] = integral("b2", t0, T);
        K2[1][1] = integral("b3", t0, T);
        //осуществим суммирование
        for (j = 0; j < sigma.length; j++)
            for (k = 0; k < sigma[j].length; k++)
            {
                sigma[j][k] = matrix_sum(matrix_mult(matrix_mult(F(T),2*A),K2),K1)[j][k];
            }
        return sigma;
    }

    //метод для оценки сигма - функции
    public void evaluation(){
        //подготовка к записи в файл
        try(FileWriter writer = new FileWriter(Filename_text, false))
        {
            System.out.println("Производится запись файла по пути "+Filename_text+"...");
            for (double i = t0; i < T; i+=h) {
                //сформируем сигма - матрицу
                double[][] sigma = form_sigma(t0, i);
                //подсчет нормы матрицы
                double norma = euc_norm(sigma);
                //сравнение с целью получения точек экстремума
                if (norma > g_max) {
                    g_max = norma;
                    g_max_t = i;
                }
                if (norma < g_min) {
                    g_min = norma;
                    g_min_t = i;
                }
                dictionary.put(i, norma);
                //запись данных в буфер
                writer.write("t = "+i+", ||δ|| = "+norma+"\r\n");
                writer.write("δ: "+ Arrays.deepToString(sigma)+"\r\n\r\n");
            }
            //непосредственно запись в файл
            writer.flush();
            System.out.println("Запись завершена.");
            file_graphic();
            System.out.println("Минимум достигнут в точке ("+g_min_t+"; "+g_min+")");
            System.out.println("Максимум достигнут в точке ("+g_max_t+"; "+g_max+")");
        }
        catch(IOException ex){
            System.out.println(ex.getMessage());
        }
    }
}
