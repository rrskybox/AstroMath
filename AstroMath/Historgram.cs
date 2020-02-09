using System;
using System.Drawing;

namespace AstroMath
{
    public static class Historgram
    {
        public static Bitmap Stretch16(Bitmap imageIn)
        {
            // g(x,y) = ((f(x,y)-fmin)/(fmax-fmin))*(2^graylevels)
            Bitmap imageOut = new Bitmap(imageIn);
            //step 1: compute fmin,fmax
            //Int16 fVal = 0;
            //Int16 fMin = 0;
            //Int16 fMax = Convert.ToInt16((System.Math.Pow(2, 16)) - 1);
            //for (int iX = 0; iX < imageIn.Width)
            //    for (int iY = 0; iY < imageIn.Height)
            //    {
            //        fVal = imageIn.GetPixel(iX, iY);
            //        if (fVal < fMin) fMin = fVal;
            //        if (fVal > fMax) fMax = fVal;
            //    }



            return imageOut;
        }

        public static Bitmap Equalize(Bitmap imageIn)
        {
            Bitmap imageOut = new Bitmap(imageIn);
            return imageOut;
        }

    }
}
