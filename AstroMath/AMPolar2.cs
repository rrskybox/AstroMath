using System;
using System.Drawing;

namespace AstroMath
{
    public partial class Polar2D
    {
        public static Point XYRotation(Point xy, double rotation)
        {
            //Result of rotation of point xy through rotation
            // x// = xcos(r) + ysin(r)
            double rotX = ((double)xy.X * Math.Cos(rotation)) + ((double)xy.Y * Math.Sin(rotation));
            // y// = -xsin(r) + ycos(r)
            double rotY = -((double)xy.X * Math.Sin(rotation)) + ((double)xy.Y * Math.Cos(rotation));
            return new Point((int)rotX, (int)rotY);
        }

        public static Point XYTranslation(Point newOrigin, Point xy)
        {
            //Translates an xy point to new origin coordinates
            return new Point(xy.X - newOrigin.X, xy.Y - newOrigin.Y);
        }

        public class RotationMatrix
        {
            double R11, R12, R21, R22;

            public RotationMatrix(double a11, double a12, double a21, double a22)
            {
                R11 = a11;
                R12 = a12;
                R21 = a21;
                R22 = a22;
            }

            public (double, double) Rotate(double a11, double a12)
            {
                //Rotate a11, a12 through matrix
                return (R11 * a11 + R21 * a12, R12 * a11 + R22 * a12);
            }
        }
    }

}

