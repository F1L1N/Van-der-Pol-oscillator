import java.io.IOException;
import java.util.Scanner;

public class Program {
    public static void main(String[] args) {
        int flag = -1;
        while (flag != 0) {
            System.out.println("Выберите режим: ");
            System.out.println("---------------------------------------");
            System.out.println("1. Формирование sigma и создание файлов");
            System.out.println("2. Вывод графика");
            System.out.println("0. Выход");
            System.out.println("-------------------------------");

            flag = new Scanner(System.in).nextInt();
            switch (flag) {
                case 1:
                    //выполнение процедуры оценивания c выводом затраченного на это времени
                    SigmaFunction sigma = new SigmaFunction();
                    long startTime = System.currentTimeMillis();
                    sigma.evaluation();
                    long duration = System.currentTimeMillis() - startTime;
                    System.out.println("Затраченное время: " + duration / 1000 + " сек.");
                    break;
                case 2:
                    //вызов метода demonstration из класса Paint
                    try {
                        Paint.demonstration(640, 480);
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                    break;
                case 0:
                    flag = 0;
                    break;
            }
        }
    }
}