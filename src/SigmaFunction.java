import com.udojava.evalex.Expression;

import java.io.FileWriter;
import java.io.IOException;
import java.math.BigDecimal;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Scanner;

/**
 * Класс сигма - функции
 *
 */

class SigmaFunction {
    /**
     * амплитуда
     */
    private double A;
    /**
     * угловая частота
     */
    private double omega;
    /**
     * шаг выполнения оценки
     */
    private double h;
    /**
     * управляющий параметр
     */
    private double nu;
    /**
     * левая граница рабочего промежутка
     */
    private double t0;
    /**
     * правая граница рабочего промежутка
     */
    private double T;
    /**
     * точность оценки
     */
    private static int SCALE = 10;
    /**
     * переменная для сравнивания чисел типа double
     */
    private double eps = Math.pow(10, -3);
    /**
     * пути создания файлов
     */
    private static String Filename_text = System.getProperty("user.dir") + "\\data\\evaluation.txt";
    private static String Filename_graph_X = System.getProperty("user.dir") + "\\data\\x.txt";
    private static String Filename_graph_Y = System.getProperty("user.dir") + "\\data\\y.txt";
    /**
     * словарь для хранения пар "X - Y координаты"
     */
    private static HashMap<Double, BigDecimal> dictionary = new HashMap<>();

    /**
     * Конструктор класса с уже определенными параметрами
     *
     * @param s1 - амплитуда
     * @param s2 - угловая частота
     * @param s3 - шаг выполнения оценки
     * @param s4 - управляющий параметр
     * @param s5 - левая граница рабочего промежутка
     * @param s6 - правая граница рабочего промежутка
     */
    SigmaFunction(double s1, double s2, double s3, double s4, double s5, double s6) {
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
    SigmaFunction() {
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
    private double[] nuAnalysis() {
        double[] lambda = new double[2];
        if (nu >= 2) {
            lambda[0] = (nu + Math.sqrt(Math.pow(nu, 2) - 4)) / 2;
            lambda[1] = (nu - Math.sqrt(Math.pow(nu, 2) - 4)) / 2;
        } else {
            System.out.print("nu из неверного промежутка");
        }
        return lambda;
    }

    /**
     * Возращает массив с коэффициентами для получения фундаментальной матрицы
     *
     * @param t текущее значение времени
     */
    private BigDecimal[] B(double t) {
        double[] lambda = nuAnalysis();
        BigDecimal b0, b1;
        if (Math.abs(lambda[0] - lambda[1]) >= eps) {
            b1 = (exp(BigDecimal.valueOf(lambda[0] * t)).subtract(exp(BigDecimal.valueOf(lambda[1] * t)))).divide(BigDecimal.valueOf(lambda[0] - lambda[1]), SCALE, BigDecimal.ROUND_HALF_UP);
            b0 = exp(BigDecimal.valueOf(lambda[0] * t)).subtract(b1.multiply(BigDecimal.valueOf(lambda[0])));
        } else {
            b1 = exp(BigDecimal.valueOf(lambda[0] * t)).multiply(BigDecimal.valueOf(t));
            b0 = exp(BigDecimal.valueOf(lambda[0] * t)).subtract(b1.multiply(BigDecimal.valueOf(lambda[0])));
        }
        return new BigDecimal[]{b0, b1};
    }


    /**
     * Возвращает факториал числа
     *
     * @param n число для подсчета
     */
    private static long fact(int n) {
        long ret = 1;
        for (int i = 1; i <= n; ++i) ret *= i;
        return ret;
    }

    /**
     * Возвращает экспоненту
     * числа типа BigDecimal
     *
     * @param t число для подсчета
     */
    private BigDecimal exp(BigDecimal t) {
        BigDecimal result = BigDecimal.valueOf(0);
        //ряд Тейлора для гиперболического синуса и косинуса
        for (int i = 1; i < 2 * SCALE; i += 2) {
            result = result.add((t.pow(i)).divide(BigDecimal.valueOf(fact(i)), SCALE, BigDecimal.ROUND_HALF_UP));
            result = result.add((t.pow(i - 1)).divide(BigDecimal.valueOf(fact(i - 1)), SCALE, BigDecimal.ROUND_HALF_UP));
        }
        return result;
    }

    /**
     * Возращает массив с коэффициентами
     * для получения фундаментальной матрицы
     * для последующего парсинга
     */
    private String[] B() {
        double[] lambda = nuAnalysis();
        String b0, b1;
        if (Math.abs(lambda[0] - lambda[1]) >= eps) {
            b1 = "((e^(" + lambda[0] + "*t)-e^(" + lambda[1] + "*t))/(" + lambda[0] + "-" + lambda[1] + "))";
            b0 = "(e^(" + lambda[0] + "*t)-" + b1 + "*" + lambda[0] + ")";
        } else {
            b1 = "(t*e^(" + lambda[0] + "*t))";
            b0 = "(e^(" + lambda[0] + "*t)-" + b1 + "*" + lambda[0] + ")";
        }
        return new String[]{b0, b1};
    }

    /**
     * Возращает фундаментальную матрицу
     *
     * @param t текущее значение времени
     */
    private BigDecimal[][] F(double t) {
        BigDecimal[] b = B(t);
        return matrixSum(matrixMultiplication(E(), b[0]), matrixMultiplication(A(), b[1]));
    }

    /**
     * Возращает фундаментальную матрицу
     * для последующего парсинга
     */
    private String[][] F() {
        String[][] F_reverse = new String[2][2];
        String[] b = B();
        F_reverse[0][0] = "(" + b[0] + ")";
        F_reverse[0][1] = "(" + b[1] + ")";
        F_reverse[1][0] = "((-1)*" + b[1] + ")";
        F_reverse[1][1] = "(" + b[0] + "+" + nu + "*" + b[1] + ")";
        return F_reverse;
    }

    /**
     * Возращает матрицу, обратную фундаментальной
     *
     * @param F фундаментальная матрица
     */
    private BigDecimal[][] FReverse(BigDecimal[][] F) {
        BigDecimal[][] F_reverse = new BigDecimal[2][2];
        BigDecimal det = (F[0][0].multiply(F[1][1])).subtract(F[0][1].multiply(F[1][0]));
        F_reverse[0][0] = (BigDecimal.valueOf(1).divide(det, SCALE, BigDecimal.ROUND_HALF_UP)).multiply(F[1][1]);
        F_reverse[0][1] = (BigDecimal.valueOf(1).divide(det, SCALE, BigDecimal.ROUND_HALF_UP)).multiply(F[1][0]);
        F_reverse[1][0] = (BigDecimal.valueOf(1).divide(det, SCALE, BigDecimal.ROUND_HALF_UP)).multiply(F[0][1]);
        F_reverse[1][1] = (BigDecimal.valueOf(1).divide(det, SCALE, BigDecimal.ROUND_HALF_UP)).multiply(F[0][0]);
        return F_reverse;
    }

    /**
     * Возращает матрицу, обратную фундаментальной,
     * для последующего парсинга
     */
    private String[][] FReverse() {
        String[][] F_reverse = new String[2][2];
        String[][] F_str = F();
        String det = "(" + F_str[0][0] + "*" + F_str[1][1] + "-" + F_str[0][1] + "*" + F_str[1][0] + ")";
        F_reverse[0][0] = "((1/" + det + ")*" + F_str[1][1] + ")";
        F_reverse[0][1] = "((-1)*(1/" + det + ")*" + F_str[1][0] + ")";
        F_reverse[1][0] = "((-1)*(1/" + det + ")*" + F_str[0][1] + ")";
        F_reverse[1][1] = "((1/" + det + ")*" + F_str[0][0] + ")";
        return F_reverse;
    }

    /**
     * Возвращает единичную матрицу
     */
    private BigDecimal[][] E() {
        BigDecimal[][] E = new BigDecimal[2][2];
        E[0][0] = BigDecimal.valueOf(1);
        E[0][1] = BigDecimal.valueOf(0);
        E[1][0] = BigDecimal.valueOf(0);
        E[1][1] = BigDecimal.valueOf(1);
        return E;
    }

    /**
     * Возвращает матрицу коэффициентов системы
     */
    private BigDecimal[][] A() {
        BigDecimal[][] A = new BigDecimal[2][2];
        A[0][0] = BigDecimal.valueOf(0);
        A[0][1] = BigDecimal.valueOf(1);
        A[1][0] = BigDecimal.valueOf(-1);
        A[1][1] = BigDecimal.valueOf(nu);
        return A;
    }

    /**
     * Возвращает матрицу Коши
     *
     * @param t   первый параметр Матрицы Коши
     * @param tao второй параметр Матрицы Коши
     */
    private BigDecimal[][] K(double t, double tao) {
        return matrixMultiplication(F(t), FReverse(F(tao)));
    }

    /**
     * Возвращает значение, подсчитанное после
     * синтаксического анализа строки
     *
     * @param t          фундаментальная матрица
     * @param expression строка для синтаксического анализа
     */
    private BigDecimal f(double t, String expression) {
        try {
            return new Expression(expression).with("omega", BigDecimal.valueOf(omega)).with("t", BigDecimal.valueOf(t)).eval();
        } catch (Exception e) {
            return BigDecimal.valueOf(0);
        }
    }

    /**
     * Возвращает определенный интеграл,
     * подсчитанный методом Симпсона
     *
     * @param s подынтегральное выражение
     * @param a нижний предел интегрирования
     * @param b верхний предел интегрирования
     */
    private BigDecimal integral(String s, double a, double b) {
        int n = 10;
        BigDecimal sum = BigDecimal.valueOf(0), sum2 = BigDecimal.valueOf(0);
        double[] x = new double[n];
        double h = (b - a) / n;
        for (int i = 0; i < n; i++) {
            x[i] = a + i * h;
        }
        for (int i = 1; i < n; i++) {
            sum = sum.add(f(x[i], s));
            sum2 = sum2.add(f((x[i - 1] + x[i]) / 2, s));
        }
        BigDecimal temp1 = sum.multiply(BigDecimal.valueOf(2));
        BigDecimal temp2 = (sum2.add(BigDecimal.valueOf(b))).multiply(BigDecimal.valueOf(4));
        BigDecimal temp3 = f(a, s).add(f(b, s));
        return BigDecimal.valueOf(h / 6).multiply(temp1.add(temp2.add(temp3)));
    }

    /**
     * Возвращает матрицу определенных интегралов,
     * подсчитанных методом Симпсона
     *
     * @param F матрица подынтегральных выражений
     * @param a нижний предел интегрирования
     * @param b верхний предел интегрирования
     */
    private BigDecimal[][] integral(String[][] F, double a, double b) {
        BigDecimal[][] temp = new BigDecimal[2][2];
        for (int i = 0; i < F.length; i++)
            for (int j = 0; j < F[i].length; j++)
                temp[i][j] = integral("(" + F[i][j] + ")*(sin(omega*t))", a, b);
        return temp;
    }

    /**
     * Возвращает сумму матриц
     *
     * @param a первое слагаемое
     * @param b второе слагаемое
     */
    private double[][] matrixSum(double[][] a, double[][] b) {
        double[][] result = new double[a.length][a[0].length];
        for (int i = 0; i < a.length; i++)
            for (int j = 0; j < a[i].length; j++)
                result[i][j] = a[i][j] + b[i][j];
        return result;
    }

    /**
     * Возвращает сумму матриц
     *
     * @param a первое слагаемое
     * @param b второе слагаемое
     */
    private BigDecimal[][] matrixSum(BigDecimal[][] a, BigDecimal[][] b) {
        BigDecimal[][] result = new BigDecimal[a.length][a[0].length];
        for (int i = 0; i < a.length; i++)
            for (int j = 0; j < a[i].length; j++)
                result[i][j] = a[i][j].add(b[i][j]);
        return result;
    }

    /**
     * Возвращает p - норму матрицы
     *
     * @param sigma входная матрица
     */
    private static BigDecimal pNorm(BigDecimal[][] sigma) {
        BigDecimal max = sigma[0][0];
        for (BigDecimal[] h : sigma)
            for (int j = 0; j < sigma[0].length; j++) {
                if (h[j].abs().compareTo(max) > 0) max = (h[j]).abs();
            }
        return max;
    }

    /**
     * Возвращает произведение матриц
     *
     * @param mA первый множитель
     * @param mB второй множитель
     */
    private double[][] matrixMultiplication(double[][] mA, double[][] mB) {
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
     * Возвращает произведение матриц
     *
     * @param mA первый множитель
     * @param mB второй множитель
     */
    private BigDecimal[][] matrixMultiplication(BigDecimal[][] mA, BigDecimal[][] mB) {
        BigDecimal[][] res = new BigDecimal[mA.length][mB[0].length];
        for (int i = 0; i < mA.length; i++) {
            for (int j = 0; j < mB[0].length; j++) {
                res[i][j] = BigDecimal.valueOf(0);
                for (int k = 0; k < mB.length; k++) {
                    res[i][j] = res[i][j].add(mA[i][k].multiply(mB[k][j]));
                }
            }
        }
        return res;
    }

    /**
     * Возвращает произведение матрицы на число
     *
     * @param mA матрица
     * @param mB число
     */
    private BigDecimal[][] matrixMultiplication(BigDecimal[][] mA, BigDecimal mB) {
        for (int i = 0; i < mA.length; i++) {
            for (int j = 0; j < mA[i].length; j++) {
                mA[i][j] = mA[i][j].multiply(mB);
            }
        }
        return mA;
    }

    /**
     * Возвращает произведение матрицы на число
     *
     * @param mA матрица
     * @param mB число
     */
    private BigDecimal[][] matrixMultiplication(BigDecimal[][] mA, double mB) {
        for (int i = 0; i < mA.length; i++) {
            for (int j = 0; j < mA[i].length; j++) {
                mA[i][j] = mA[i][j].multiply(BigDecimal.valueOf(mB));
            }
        }
        return mA;
    }

    /**
     * Формирует файлы с координатами для графиков
     */
    private void fileGraphic() {
        //сформируем первый файл с X - координатами
        try (FileWriter writer = new FileWriter(Filename_graph_X, false)) {
            System.out.println("Производится запись файла по пути " + Filename_graph_X + "...");
            for (Map.Entry entry : dictionary.entrySet()) {
                writer.write(entry.getKey() + "\r\n");
            }
            writer.flush();
            System.out.println("Запись завершена.");
        } catch (IOException ex) {
            System.out.println(ex.getMessage());
        }
        //сформируем первый файл с Y - координатами
        try (FileWriter writer = new FileWriter(Filename_graph_Y, false)) {
            System.out.println("Производится запись файла по пути " + Filename_graph_Y + "...");
            for (Map.Entry entry : dictionary.entrySet()) {
                writer.write(entry.getValue() + "\r\n");
            }
            writer.flush();
            System.out.println("Запись завершена.");
        } catch (IOException ex) {
            System.out.println(ex.getMessage());
        }
    }

    /**
     * Возвращает массив строк,
     * разбитый по точке
     *
     * @param str строка для обработки
     */
    private String[] split(String str) {
        String[] split = new String[2];
        int index = 0;
        for (int i = 0; i < str.length(); i++) {
            if (str.charAt(i) == '.') {
                index = i;
                break;
            }
        }
        split[0] = str.substring(0, index);
        split[1] = str.substring(index + 1, str.length());
        return split;
    }

    /**
     * Возвращает форматированную строку для вставки в файл
     *
     * @param number число для анализа
     */
    private String formStr(BigDecimal number) {
        number = number.setScale(SCALE, BigDecimal.ROUND_HALF_UP);
        String str = number.toString();
        //нет точки
        if (str.indexOf('.') == -1) {
            str = str.charAt(0) + "." + str.substring(1, SCALE - 1) + "E" + (str.length() - 1);
        }
        //точка есть
        else {
            String[] temp = split(str);
            //число типа "0.001"
            if (temp[0].length() > SCALE - 4) {
                //число типа "10000000.01"
                str = temp[0].charAt(0) + "." + temp[0].substring(1, SCALE - 4) + "E" + (temp[0].length() - 1);
            }
        }
        return str;
    }

    /**
     * Возвращает сигма - функцию
     *
     * @param t0 левая граница рабочего промежутка
     * @param T  правая граница рабочего промежутка
     */
    private BigDecimal[][] formSigma(double t0, double T) {
        return matrixMultiplication(matrixMultiplication(F(T), 3 * A), integral(FReverse(), t0, T));
    }

    /**
     * Формирует оценку погрешности линеаризации
     */
    void evaluation() {
        //подготовка к записи в файл
        try (FileWriter writer = new FileWriter(Filename_text, false)) {
            System.out.println("Производится запись файла по пути " + Filename_text + "...");
            for (double i = t0 + h; i < T; i += h) {
                //сформируем сигма - матрицу
                BigDecimal[][] sigma = formSigma(t0, i);
                //подсчет нормы матрицы
                BigDecimal norma = pNorm(sigma);
                dictionary.put(i, norma);
                //запись данных в буфер
                writer.write("t = " + i + ", ||δ|| = " + formStr(norma) + "\r\n");
                writer.write("δ: " + Arrays.deepToString(sigma) + "\r\n\r\n");
            }
            System.out.println("Завершение расчетов...");
            //непосредственно запись в файл
            writer.flush();
            System.out.println("Запись завершена.");
            fileGraphic();
        } catch (IOException ex) {
            System.out.println(ex.getMessage());
        }
    }
}
