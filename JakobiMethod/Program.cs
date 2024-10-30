using System;

public class JacobiEigenvalueAlgorithm
{
    // Метод для вывода результатов
    public static void PrintResults(double[] eigenvalues, double[,] eigenvectors, int numOfIterations)
    {
        Console.WriteLine("Собственные значения:");
        foreach (var value in eigenvalues)
            Console.Write($"{value} ");
        Console.WriteLine("\nСобственные векторы:");

        for (int i = 0; i < eigenvectors.GetLength(0); i++)
        {
            for (int j = 0; j < eigenvectors.GetLength(1); j++)
                Console.Write($"{eigenvectors[i, j]} ");
            Console.WriteLine();
        }

        Console.WriteLine($"Общее количество вращений: {numOfIterations}\n");
    }

    // Метод для вычисления суммы вне диагонали и проверки на сходимость
    public static double TresholdConvergence(double[,] a, int n, int itNum)
    {
        double offSum = 0.0;
        for (int j = 0; j < n; j++)
        {
            for (int i = 0; i < j; i++)
            {
                offSum += a[i, j] * a[i, j];
            }
        }

        offSum = Math.Sqrt(offSum) / (4.0 * n);

        if (offSum == 0.0)
            return offSum;

        return offSum;
    }

    // Метод для нахождения наибольшего внедиагонального элемента
    public static (double, int, int) MaxOffDiagonal(double[,] a, int n)
    {
        double maxElem = 0.0;
        int k = 0, l = 0;

        for (int i = 0; i < n - 1; i++)
        {
            for (int j = i + 1; j < n; j++)
            {
                if (Math.Abs(a[i, j]) >= maxElem)
                {
                    maxElem = Math.Abs(a[i, j]);
                    k = i;
                    l = j;
                }
            }
        }

        return (maxElem, k, l);
    }

    // Метод для выполнения вращения и аннулирования элемента a[k, l]
    public static void Rotate(double[,] a, double[,] p, int n, int k, int l)
    {
        double diff = a[l, l] - a[k, k];
        double t;

        if (Math.Abs(a[k, l]) < Math.Abs(diff) * 0.1)
        {
            t = a[k, l] / diff;
        }
        else
        {
            double theta = diff / (2.0 * a[k, l]);
            t = 1.0 / (Math.Abs(theta) + Math.Sqrt(theta * theta + 1.0));

            if (theta < 0.0)
                t = -t;
        }

        double c = 1.0 / Math.Sqrt(t * t + 1.0);
        double s = t * c;
        double tau = s / (1.0 + c);

        double temp = a[k, l];
        a[k, l] = 0.0;
        a[k, k] -= t * temp;
        a[l, l] += t * temp;

        // Обновление элементов матрицы
        for (int i = 0; i < k; i++)
        {
            temp = a[i, k];
            a[i, k] = temp - s * (a[i, l] + tau * temp);
            a[i, l] += s * (temp - tau * a[i, l]);
        }

        for (int i = k + 1; i < l; i++)
        {
            temp = a[k, i];
            a[k, i] = temp - s * (a[i, l] + tau * temp);
            a[i, l] += s * (temp - tau * a[i, l]);
        }

        for (int i = l + 1; i < n; i++)
        {
            temp = a[k, i];
            a[k, i] = temp - s * (a[l, i] + tau * temp);
            a[l, i] += s * (temp - tau * a[l, i]);
        }
        for (int i = 0; i < n; i++)
        {
            for(int j = 0; j < n;j++)
            {
                if(i == j)
                {
                    break;
                }
                a[i,j] = a[j,i];
            }
        }

        // Обновление матрицы преобразования
        for (int i = 0; i < n; i++)
        {
            temp = p[i, k];
            p[i, k] = temp - s * (p[i, l] + tau * temp);
            p[i, l] += s * (temp - tau * p[i, l]);
        }
    }

    // Основной метод для вычисления собственных значений и векторов методом Якоби
    public static (double[], double[,], int) JacobiEigenvalueAlgorithmMethod(double[,] a, int n, int maxIter, double eps = 0.000001)
    {
        double[,] p = new double[n, n];
        for (int i = 0; i < n; i++)
            p[i, i] = 1.0;

        for (int itNum = 0; itNum < maxIter; itNum++)
        {
            var (maxElem, k, l) = MaxOffDiagonal(a, n);
            double offSum = TresholdConvergence(a, n, itNum);

            if (maxElem < eps || offSum == 0.0)
            {
                double[] eigenvalues = new double[n];
                for (int i = 0; i < n; i++)
                    eigenvalues[i] = a[i, i];
                return (eigenvalues, p, itNum);
            }
            for (int i = 0; i < n; i++)
            {
                for(int j = 0; j < n; j++)
                {
                    Console.Write(a[i,j] + " ");
                }
                Console.WriteLine();
            }
            Console.WriteLine();
            Rotate(a, p, n, k, l);
        }

        Console.WriteLine("Метод не сошелся");
        return (null, null, maxIter);
    }

    // Основной метод для тестирования
    public static void Main()
    {
        double[,] matrix1 = {
            { 3.2, 1, 2.2 },
            { 1, 3.7, 2.2 },
            { 2.2, 2.2, 4.2 }
        };

        Console.WriteLine("Матрица 1");
        var result1 = JacobiEigenvalueAlgorithmMethod(matrix1, 3, 100);
        PrintResults(result1.Item1, result1.Item2, result1.Item3);

        double[,] matrix2 = {
            { 1.6, 1.2, -1.1 },
            { 1.2, 1.1, 0.6 },
            { -1.1, 0.6, 0.8}
        };

        Console.WriteLine("Матрица 2");
        var result2 = JacobiEigenvalueAlgorithmMethod(matrix2, 3, 100);
        PrintResults(result2.Item1, result2.Item2, result2.Item3);

        double[,] matrix3 = {
            { 1, 0, 0 },
            { 0, 5, 0 },
            { 0, 0, 3}
        };

        Console.WriteLine("Матрица 3");
        var result3 = JacobiEigenvalueAlgorithmMethod(matrix3, 3, 100);
        PrintResults(result3.Item1, result3.Item2, result3.Item3);
    }
}
