import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Map;
import java.util.Scanner;
import java.util.TreeMap;

/**
 * Класс сигма - функции
 * @author Филин Павел
 * @version 1.0
 */

class SigmaFunction {
    /** амплитуда */
    private double Ampl;
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
    /**переменные для хранения
     * максимальных и минимальных значений нормы сигмы - матрицы,
     * а также их индексов*/
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
    private static TreeMap<Double, Double> dictionary = new TreeMap<>();

    /**
     * Конструктор класса с уже определенными параметрами
     * @param s1 - амплитуда
     * @param s2 - угловая частота
     * @param s3 - шаг выполнения оценки
     * @param s4 - управляющий параметр
     * @param s5 - левая граница рабочего промежутка
     * @param s6 - правая граница рабочего промежутка
     */
    SigmaFunction (double s1, double s2, double s3, double s4, double s5, double s6){
        h = s1;
        Ampl = s2;
        omega = s3;
        nu = s4;
        t0 = s5;
        T = s6;
    }

    /**
     * Конструктор класса с вводом параметров
     */
    SigmaFunction (){
        System.out.print("Ввод шага: ");
        h = new Scanner(System.in).nextDouble();
        System.out.print("Ввод амплитуды: ");
        Ampl = new Scanner(System.in).nextDouble();
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
    private double[] nuAnalysis(){
        double[] lambda = new double[2];
        if ((nu >= 0)&&(nu < 2)) {
            Complex lambda1 = new Complex(nu/2, Math.sqrt(Math.abs(Math.pow(nu,2)-4))/2);
            Complex lambda2 = new Complex(nu/2, -Math.sqrt(Math.abs(Math.pow(nu,2)-4))/2);
        }else
            if (nu >= 2){
                lambda[0] = (nu + Math.sqrt(Math.pow(nu,2)-4))/2;
                lambda[1] = (nu - Math.sqrt(Math.pow(nu,2)-4))/2;
            }else{
                System.out.print("nu из неверного промежутка");
            }
        //System.out.println(Arrays.toString(lambda));
        return lambda;
    }

    /**
     * Возращает массив с коэффициентами для получения фундаментальной матрицы
     * @param t текущее значение времени
     */
    private double[] formB (double t){
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
    private String[] formBStr (){
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
    private double[][] F (double t){
        double[] b = formB(t);
        return matrix_sum(matrix_mult(E(),b[0]),matrix_mult(A(),b[1]));
    }

    /**
     * Возращает фундаментальную матрицу
     * для последующего парсинга
     */
    private String[][] FStr(){
        String[][] F_reverse = new String[2][2];
        String[] b = formBStr();
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
    private double[][] F_reverse (double[][] F){
        double[][] F_reverse = new double[2][2];
        double det = F[0][0]*F[1][1] - F[0][1]*F[1][0];
        F_reverse[0][0] = (1/det) * F[1][1];
        F_reverse[0][1] = (1/det) * F[1][0];
        F_reverse[1][0] = (1/det) * F[0][1];
        F_reverse[1][1] = (1/det) * F[0][0];
        /*F_reverse[0][0] = (1+t)/Math.exp(t);
        F_reverse[0][1] = -t/Math.exp(t);
        F_reverse[1][0] = t/Math.exp(t);
        F_reverse[1][1] = (1-t)/Math.exp(t);*/
        return F_reverse;
    }

    /**
     * Возращает матрицу, обратную фундаментальной,
     * для последующего парсинга
     */
    private String[][] FReverseStr(){
        String[][] F_reverse = new String[2][2];
        String[][] F_str = FStr();
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
    private double[][] E(){
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
    private double[][] A(){
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
    private double[][] K (double t, double tao){
        return matrix_mult(F(t), F_reverse(F(tao)));
    }

    /**
     * Возвращает значение, подсчитанное после
     * синтаксического анализа строки
     * @param t фундаментальная матрица
     * @param expression строка для синтаксического анализа
     */
    private double f(double t, String expression){
        MathParser parser = new MathParser(Ampl, omega, t);
        try {
            return parser.Parse(expression);
        } catch (Exception e) {
            e.printStackTrace();
            return 0.0;
        }
        /*switch (s){
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
        }*/
    }

    /**
     * Возвращает определенный интеграл,
     * подсчитанный методом Симпсона
     * @param s подынтегральное выражение
     * @param a нижний предел интегрирования
     * @param b верхний предел интегрирования
     */
    private double integral(String s, double a, double b){
        int n = 50;
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
    private double[][] integral(String[][] F, double a, double b){
        double[][] temp = new double[2][2];
        for (int i = 0; i < F.length; i++)
            for (int j = 0; j < F[i].length; j++)
                temp[i][j] = integral("("+F[i][j]+")*(A*sin(omega*t))",a,b);
        return temp;
    }

    /**
     * Возвращает сумму матриц
     * @param a первое слагаемое
     * @param b второе слагаемое
     */
    private double[][] matrix_sum(double[][] a, double[][] b){
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
    private static double euc_norm(double[][] sigma){
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

    /**
     * Возвращает произведение матрицы на число
     * @param mA матрица
     * @param mB число
     */
    private double[][] matrix_mult(double[][] mA, double mB){
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

    /**
     * Возвращает сигма - функцию
     * @param t0 левая граница рабочего промежутка
     * @param T правая граница рабочего промежутка
     */
    //метод для получения сигма - функции
    private double[][] form_sigma(double t0, double T){
        double[][] sigma = new double[2][2];
        double[][] F = matrix_mult(F(T),3);
        double[][] I = integral(FReverseStr(),t0,T);
        for (int i = 0; i < sigma.length; i++)
            for (int j = 0; j < sigma[i].length; j++)
                sigma[i][j] = matrix_mult(F, I)[i][j];
        /*double[][] sigma = new double[2][2];
        int j, k, p = 0;
        //сформируем первое слагаемое из формулы
        double[][] K1 = new double[2][2];
        for (j = 0; j < sigma.length; j++)
            for (k = 0; k < sigma[0].length; k++)
            {
                p++;
                K1[j][k] = Ampl *integral("a"+p, t0, T);
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
                sigma[j][k] = matrix_sum(matrix_mult(matrix_mult(F(T),2* Ampl),K2),K1)[j][k];
            }*/
        return sigma;
    }

    /**
     * Формирует оценку погрешности линеаризации
     */
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
