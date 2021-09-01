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

    }

}

