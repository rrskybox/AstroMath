using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace AstroMath
{
    class MatrixMath
    {
        public static double DotProduct(double[] a, double[] b)
        {
            double scalor = 0;
            if (a.Length == b.Length)
            {
                for (int i = 0; i < a.Length; i++)
                { scalor += a[i] * b[i]; }
            }
            return scalor;
        }

        public static double[,] MatrixProduct(double[,] a, double[,] b)
        {
            //figure out sizes of a and b.  Both should be rank = 2
            int ai = a.GetLength(0);  //rows of a
            int aj = a.GetLength(1);  //columns of a
            int bi = b.GetLength(0);  //rows of b
            int bj = b.GetLength(1);  //columns of b
            double[,] c = new double[ai, bj];
            if (aj == bi)
            {
                for (int i = 0; i < ai; i++)
                    for (int j = 0; j < aj; j++)
                        for (int k = 0; k < aj; k++) { c[i, j] += a[i, k] * b[k, i]; }
            }
            return c;
        }
               
        public static double[] CrossProduct(double[] a, double[] b)
        {
            //Must be R3 space
            double[] vector = new double[3];
            vector[0] = a[1] * b[2] - a[2] * b[1];
            vector[1] = a[2] * b[0] - a[0] * b[2];
            vector[2] = a[0] * b[1] - a[1] * b[0];
            return vector;
        }

    }
}
