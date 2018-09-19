using System;

namespace AstroMath
{
    public class Formatters
    {
        public static string HourString(double dvalue)
        //Converts a double value (dvalue) to a string looking like an hour:minutes
        {
            int hr = (int)Math.Truncate(dvalue);
            int min = (int)Math.Truncate((dvalue - hr) * 60);
            return (hr.ToString() + ":" + min.ToString());
        }

     }
}
