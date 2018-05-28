import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Scanner;

/**
 * Класс сигма - функции
 * @version 1.0
 */

class SigmaFunction
{
    /** амплитуда */
    private double A;
    /** угловая частота */
    private double omega;
    /** шаг выполнения оценки */
    private double h;
    /** управляющий параметр */
    private double nu;
    /** левая граница рабочего промежутка */
    private double t0;
    /** правая граница рабочего промежутка */
    private double T;
    /** переменные для хранения
     *  максимальных и минимальных значений нормы сигмы - матрицы,
     *  а также их индексов*/
    private double g_max = -Math.pow(10,6);
    private double g_max_t;
    private double g_min = Math.pow(10,6);
    private double g_min_t;
    /**
     * переменная для сравнивания чисел типа double
     */
    private double eps = Math.pow(10,-3);
    /**
     * пути создания файлов
     */
    private static String Filename_text = System.getProperty("user.dir")+"\\data\\evaluation.txt";
    private static String Filename_graph_X = System.getProperty("user.dir")+"\\data\\x.txt";
    private static String Filename_graph_Y = System.getProperty("user.dir")+"\\data\\y.txt";
    /**
     * словарь для хранения пар "X - Y координаты"
     */
    private static HashMap<Double, Double> dictionary = new HashMap<>();

    /**
     * Конструктор класса с уже определенными параметрами
     * @param s1 - амплитуда
     * @param s2 - угловая частота
     * @param s3 - шаг выполнения оценки
     * @param s4 - управляющий параметр
     * @param s5 - левая граница рабочего промежутка
     * @param s6 - правая граница рабочего промежутка
     */
    SigmaFunction (double s1, double s2, double s3, double s4, double s5, double s6)
    {
        h = s1;
        A = s2;
        omega = s3;
        nu = s4;
        t0 = s5;
        T = s6;
    }

    /**
     * Конструктор класса с вводом параметров
     */
    SigmaFunction ()
    {
        System.out.print("Ввод шага: ");
        h = new Scanner(System.in).nextDouble();
        System.out.print("Ввод амплитуды: ");
        A = new Scanner(System.in).nextDouble();
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

    /**
     * Анализирует nu и возвращает собственные значения системы
     */
    private double[] nuAnalysis()
    {
        double[] lambda = new double[2];
        if (nu >= 2) {
            lambda[0] = (nu + Math.sqrt(Math.pow(nu, 2) - 4)) / 2;
            lambda[1] = (nu - Math.sqrt(Math.pow(nu, 2) - 4)) / 2;
        }else{
            System.out.print("nu из неверного промежутка");
        }
        return lambda;
    }

    /**
     * Возращает массив с коэффициентами для получения фундаментальной матрицы
     * @param t текущее значение времени
     */
    private double[] B(double t)
    {
        double[] lambda = nuAnalysis();
        double b0, b1;
        if (Math.abs(lambda[0]-lambda[1]) >= eps){
            b1 = (Math.exp(lambda[0]*t) - Math.exp(lambda[1]*t))/(lambda[0]-lambda[1]);
            b0 = Math.exp(lambda[0]*t) - b1 * lambda[0];
        }else{
            b1 = t * Math.exp(lambda[0] * t);
            b0 = Math.exp(lambda[0] * t) - b1 * lambda[0];
        }
        return new double[]{b0, b1};
    }

    /**
     * Возращает массив с коэффициентами
     * для получения фундаментальной матрицы
     * для последующего парсинга
     */
    private String[] B()
    {
        double[] lambda = nuAnalysis();
        String b0, b1;
        if (Math.abs(lambda[0]-lambda[1]) >= eps){
            b1 = "((e^("+lambda[0]+"*t)-e^("+lambda[1]+"*t))/("+lambda[0]+"-"+lambda[1]+"))";
            b0 = "(e^("+lambda[0]+"*t)-"+b1+"*"+lambda[0]+")";
        }else{
            b1 = "(t*e^("+lambda[0]+"*t))";
            b0 = "(e^("+lambda[0]+"*t)-"+b1+"*"+lambda[0]+")";
        }
        return new String[]{b0, b1};
    }

    /**
     * Возращает фундаментальную матрицу
     * @param t текущее значение времени
     */
    private double[][] F(double t)
    {
        double[] b = B(t);
        return matrixSum(matrixMultiplication(E(),b[0]), matrixMultiplication(A(),b[1]));
    }

    /**
     * Возращает фундаментальную матрицу
     * для последующего парсинга
     */
    private String[][] F()
    {
        String[][] F_reverse = new String[2][2];
        String[] b = B();
        F_reverse[0][0] = "("+b[0]+")";
        F_reverse[0][1] = "("+b[1]+")";
        F_reverse[1][0] = "((-1)*" + b[1]+")";
        F_reverse[1][1] = "("+b[0]+"+"+nu+"*"+b[1]+")";
        return F_reverse;
    }

    /**
     * Возращает матрицу, обратную фундаментальной
     * @param F фундаментальная матрица
     */
    private double[][] FReverse(double[][] F)
    {
        double[][] F_reverse = new double[2][2];
        double det = F[0][0]*F[1][1] - F[0][1]*F[1][0];
        F_reverse[0][0] = (1/det) * F[1][1];
        F_reverse[0][1] = (1/det) * F[1][0];
        F_reverse[1][0] = (1/det) * F[0][1];
        F_reverse[1][1] = (1/det) * F[0][0];
        return F_reverse;
    }

    /**
     * Возращает матрицу, обратную фундаментальной,
     * для последующего парсинга
     */
    private String[][] FReverse()
    {
        String[][] F_reverse = new String[2][2];
        String[][] F_str = F();
        String det ="("+F_str[0][0]+"*"+F_str[1][1]+"-"+F_str[0][1]+"*"+F_str[1][0]+")";
        F_reverse[0][0] = "((1/"+det+")*"+F_str[1][1]+")";
        F_reverse[0][1] = "((-1)*(1/"+det+")*"+F_str[1][0]+")";
        F_reverse[1][0] = "((-1)*(1/"+det+")*"+F_str[0][1]+")";
        F_reverse[1][1] = "((1/"+det+")*"+F_str[0][0]+")";
        return F_reverse;
    }

    /**
     * Возвращает единичную матрицу
     */
    private double[][] E()
    {
        double[][] E = new double[2][2];
        E[0][0] = 1;
        E[0][1] = 0;
        E[1][0] = 0;
        E[1][1] = 1;
        return E;
    }

    /**
     * Возвращает матрицу коэффициентов системы
     */
    private double[][] A()
    {
        double[][] A = new double[2][2];
        A[0][0] = 0;
        A[0][1] = 1;
        A[1][0] = -1;
        A[1][1] = nu;
        return A;
    }

    /**
     * Возвращает матрицу Коши
     * @param t первый параметр Матрицы Коши
     * @param tao второй параметр Матрицы Коши
     */
    private double[][] K (double t, double tao)
    {
        return matrixMultiplication(F(t), FReverse(F(tao)));
    }

    /**
     * Возвращает значение, подсчитанное после
     * синтаксического анализа строки
     * @param t фундаментальная матрица
     * @param expression строка для синтаксического анализа
     */
    private double f(double t, String expression)
    {
        MathParser parser = new MathParser();
        parser.setVariable("omega",omega);
        parser.setVariable("t",t);
        try {
            double result = parser.Parse(expression);
            if (Double.isNaN(result)||Double.isInfinite(result))
                return 0.0;
            else return result;
        } catch (Exception e) {
            e.printStackTrace();
            return 1.0;
        }
    }

    /**
     * Возвращает определенный интеграл,
     * подсчитанный методом Симпсона
     * @param s подынтегральное выражение
     * @param a нижний предел интегрирования
     * @param b верхний предел интегрирования
     */
    private double integral(String s, double a, double b)
    {
        int n = 10;
        double sum = 0, sum2 = 0;
        double[] x = new double[n];
        double h = (b-a)/n;
        for(int i = 0; i < n; i++){
            x[i] = a+i*h;
        }
        for(int i = 1; i < n; i++){
            sum += f(x[i],s);
            sum2 += f((x[i-1]+x[i])/2,s);
        }
        return h/6*(f(a,s)+f(b,s)+2*sum+4*(sum2+b));
    }

    /**
     * Возвращает матрицу определенных интегралов,
     * подсчитанных методом Симпсона
     * @param F матрица подынтегральных выражений
     * @param a нижний предел интегрирования
     * @param b верхний предел интегрирования
     */
    private double[][] integral(String[][] F, double a, double b)
    {
        double[][] temp = new double[2][2];
        for (int i = 0; i < F.length; i++)
            for (int j = 0; j < F[i].length; j++)
                temp[i][j] = integral("("+F[i][j]+")*(sin(omega*t))",a,b);
        return temp;
    }

    /**
     * Возвращает сумму матриц
     * @param a первое слагаемое
     * @param b второе слагаемое
     */
    private double[][] matrixSum(double[][] a, double[][] b)
    {
        double[][] result = new double[a.length][a[0].length];
        for (int i = 0; i < a.length; i++)
            for (int j = 0; j < a[i].length; j++)
                result[i][j] = a[i][j] + b[i][j];
        return result;
    }

    /**
     * Возвращает евклидову норму матрицы
     * @param sigma входная матрица
     */
    private static double eucNorm(double[][] sigma)
    {
        double sum = 0;
        for (double[] h : sigma)
            for (int j = 0; j < sigma[0].length;j++) {
                sum += Math.pow(h[j], 2);
            }
        return Math.sqrt(sum);
    }

    /**
     * Возвращает произведение матриц
     * @param mA первый множитель
     * @param mB второй множитель
     */
    private double[][] matrixMultiplication(double[][] mA, double[][] mB)
    {
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

    /**
     * Возвращает произведение матрицы на число
     * @param mA матрица
     * @param mB число
     */
    private double[][] matrixMultiplication(double[][] mA, double mB)
    {
        for (int i = 0; i < mA.length; i++) {
            for (int j = 0; j < mA[i].length; j++) {
                mA[i][j] = mA[i][j] * mB;
            }
        }
        return mA;
    }

    /**
     * Формирует файлы с координатами для графиков
     */
    private void fileGraphic()
    {
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

    /**
     * Возвращает сигма - функцию
     * @param t0 левая граница рабочего промежутка
     * @param T правая граница рабочего промежутка
     */
    private double[][] formSigma(double t0, double T)
    {
        return matrixMultiplication(matrixMultiplication(F(T),3*A), integral(FReverse(),t0,T));
    }

    /**
     * Формирует оценку погрешности линеаризации
     */
    public void evaluation()
    {
        //подготовка к записи в файл
        try(FileWriter writer = new FileWriter(Filename_text, false))
        {
            System.out.println("Производится запись файла по пути "+Filename_text+"...");
            for (double i = t0; i < T; i+=h) {
                //сформируем сигма - матрицу
                double[][] sigma = formSigma(t0, i);
                //подсчет нормы матрицы
                double norma = eucNorm(sigma);
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
            fileGraphic();
            System.out.println("Минимум достигнут в точке ("+g_min_t+"; "+g_min+")");
            System.out.println("Максимум достигнут в точке ("+g_max_t+"; "+g_max+")");
        }
        catch(IOException ex){
            System.out.println(ex.getMessage());
        }
    }
}
