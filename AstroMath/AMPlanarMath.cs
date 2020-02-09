using System;
using System.Drawing;

namespace AstroMath
{
    public partial class Planar
    {
        public class QuadRoot
        {
            //Structure for returning quadratic root results
            public int nz;
            public double xe;
            public double ye;
            public double zero1;
            public double zero2;

            public QuadRoot()
            {
                nz = 0;  //Number of roots within the interval [-1,1]
                xe = 0;  //Extreme value of X in parabola solution
                ye = 0;  //Extreme value of Y in parabola solution
                zero1 = 0;  //First root within [-1,1] for NZ = 1,2
                zero2 = 0;  //Second root within [-1,1] for NZ = 2
            }
        }

        public static QuadRoot Quad(double yminus, double yzero, double yplus)
        {

            //Finds a parabola through three points 
            //   (-1,yminus), (0,yzero), (1,yplus)
            //   that do not lie on a straight line.

            double a, b, c, dis, dx;
            QuadRoot qr = new QuadRoot();

            qr.nz = 0;
            a = 0.5 * (yminus + yplus) - yzero;
            b = 0.5 * (yplus - yminus);
            c = yzero;
            qr.xe = -b / (2 * a);
            qr.ye = ((a * qr.xe + b) * qr.xe) + c;
            dis = Math.Pow(b, 2) - (4 * a * c);
            if (dis >= 0)
            {
                dx = 0.5 * Math.Sqrt(dis) / Math.Abs(a);
                qr.zero1 = qr.xe - dx;
                qr.zero2 = qr.xe + dx;
                if (Math.Abs(qr.zero1) <= 1)
                {
                    qr.nz = qr.nz + 1;
                }
                if (Math.Abs(qr.zero2) <= 1)
                {
                    qr.nz = qr.nz + 1;
                }
                if (qr.zero1 < -1)
                {
                    qr.zero1 = qr.zero2;
                }
            }
            return qr;
        }

        public static double Frac(double x)
        {
            //   returns fraction Of A less than 1 As positive value
            if (x < 0)
            { x = Math.Abs(x - (int)x); }
            else
            { x = x - (int)x; }
            if (x < 0)
            { return x + 1; }
            else
            { return x; }
        }

        public static Point ThirdPoint(Point C, double circleradius, double Alpha, double ht)
        {
            //Calculates the coordinations (point) for the third point of a isocolese triangle with
            // a height of ht and rotated to an angle (radians)
            double P = Math.Sqrt(Math.Pow(ht, 2) + Math.Pow((circleradius / 2), 2));
            double Beta = Math.Sin(ht / P);
            Point T = new Point((int)(C.X + P * Math.Cos(Alpha + Beta)), (int)(C.Y + P * Math.Sin(Alpha + Beta)));
            return T;
        }

        public static double DotProduct(double[] a, double[] b)
        {
            //Computes the dot product of two vectors
            if (a.Length != b.Length) return 0;
            double dp = 0;
            for (int i = 0; i < a.Length; i++) { dp += a[i] * b[i]; }
            return dp;
        }
    }
}
