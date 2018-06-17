import java.math.BigDecimal;

public class Matrix {
    private BigDecimal[][] matrix = new BigDecimal[2][2];

    BigDecimal[][] getMatrix() {
        return matrix;
    }

    void setMatrix(BigDecimal[][] Matrix) {
        this.matrix = Matrix;
    }

    /**
     * Возвращает произведение матриц
     *  @param mA первый множитель
     * @param mB второй множитель
     */
    BigDecimal[][] matrixMultiplication(BigDecimal[][] mA, BigDecimal[][] mB) {
        for (int i = 0; i < mA.length; i++) {
            for (int j = 0; j < mB[0].length; j++) {
                matrix[i][j] = BigDecimal.valueOf(0);
                for (int k = 0; k < mB.length; k++) {
                    matrix[i][j] = matrix[i][j].add(mA[i][k].multiply(mB[k][j]));
                }
            }
        }
        return mA;
    }

    /**
     * Возвращает произведение матрицы на число
     * типа BigDecimal
     *
     * @param mA матрица
     * @param mB число
     */
    BigDecimal[][] matrixMultiplication(BigDecimal[][] mA, BigDecimal mB) {
        for (int i = 0; i < mA.length; i++) {
            for (int j = 0; j < mA[i].length; j++) {
                matrix[i][j] = mA[i][j].multiply(mB);
            }
        }
        return matrix;
    }

    /**
     * Возвращает произведение матрицы на число
     * типа double
     *
     * @param mA матрица
     * @param mB число
     */
    BigDecimal[][] matrixMultiplication(BigDecimal[][] mA, double mB) {
        for (int i = 0; i < mA.length; i++) {
            for (int j = 0; j < mA[i].length; j++) {
                matrix[i][j] = mA[i][j].multiply(BigDecimal.valueOf(mB));
            }
        }
        return matrix;
    }

    /**
     * Возвращает сумму матриц
     *
     * @param a первое слагаемое
     * @param b второе слагаемое
     */
    BigDecimal[][] matrixSum(BigDecimal[][] a, BigDecimal[][] b) {
        for (int i = 0; i < a.length; i++)
            for (int j = 0; j < a[i].length; j++)
                matrix[i][j] = a[i][j].add(b[i][j]);
        return matrix;
    }
}
