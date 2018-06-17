import java.math.BigDecimal;

public class Line {
    private String str;

    Line(String line) {
        str = line;
    }

    String getStr() {
        return str;
    }

    void setStr(String str) {
        this.str = str;
    }

    /**
     * Конструктор класса с уже определенными параметрами
     *
     * @param number - число для обработки
     * @param SCALE - точность
     */
    public Line(BigDecimal number, int SCALE) {
        number = number.setScale(SCALE, BigDecimal.ROUND_HALF_UP);
        str = number.toString();
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
}
