using System;

namespace AstroMath
{
    public partial class Transform

    {
        #region Polar Trigometric Functions

        //SIND(Degrees as double) -> Sine as double:
        //   Sine Function in() degrees
        public static double SinD(double x)
        { //-1 to 1
            return Math.Sin(x * Celestial.TWOPI / 360.0);
        }

        //COSD(Degrees as double) -> Cosine as double:
        //   Cosine Function in() degrees
        public static double CosD(double x)
        { //Degrees
            return Math.Cos(x * Celestial.TWOPI / 360.0);
        }

        //RADIANSTODEGREES(Radians as double) -> Degrees as double:
        //   Convert Radians To degrees
        public static double RadiansToDegrees(double rad)
        {
            return (rad * (360.0 / Celestial.TWOPI));
        }

        //DEGREESTORADIANS(Degrees as double) -> Radians as double:
        //   Convert Degrees To radians
        public static double DegreesToRadians(double deg)
        {
            return ((deg % 360.0) * (Celestial.TWOPI / 360.0));
        }

        //HOURSTORADIANS(Hours as double) -> Radians as double:
        //   Convert hours(0 - 24) To radians
        public static double HoursToRadians(double hours)
        {
            //(radians)
            //Convert hours in timespan to radians (15 degrees per hour)
            //Normalize hours 0 - 24,) { convert to radians
            return ((hours % 24.0) / 24.0) * Celestial.TWOPI;
        }

        //RADIANSTOHOURS(Radians as double) -> hours as double:
        //   Convert radians To hours (-24,+24)
        public static double RadiansToHours(double radians)
        {
            //Convert radians to hours (15 degrees).
            return ((radians / Celestial.TWOPI) * 24.0) % 24.0;
        }

        //HOURSTODEGREES(Hours as double) -> Degrees as double:
        //   Convert hours(-24,+24) to degrees (-360, +360)
        public static double HoursToDegrees(double hours)
        {
            //Convert hours in timespan to degrees (15 degrees per hour)
            //Normalize hours 0 - 24,) { convert to degrees 0 - 360.0
            return ((hours % 24.0) / 24.0) * 360.0;
        }

        //DEGREESTOHOURS(Degrees as double) -> Hours as double:
        //   Convert degrees(-360,+360) To hours (-24,+24)
        public static double DegreesToHours(double degrees)
        {
            //Convert degrees to hours (15 degrees).
            return ((degrees / 360.0) * 24.0) % 24.0;
        }

        /// <summary>
        /// HourAngleToPolarAngle translates HourAngle (0 hour at 6 oclock) in hours
        ///     to Polar Coordinate (0 radians at 3 oclock) in radians (-2pi,+2pi)
        /// </summary>
        /// <param name="HA"></param>
        public static double HourAngleToPolarAngle(double haH)
        {
            double haR = HoursToRadians(haH);
            double ha6R = HoursToRadians(6.0);
            double haTrans = -(haR - ha6R);
            haTrans = NormalizeRadianRange(haTrans);
            return (haTrans);
        }
        #endregion

        #region Normalizing Methods

        /// <summary>
        ///  Converts angle in degrees (open) to degrees (0,360)
        /// </summary>
        /// <param name="angleD"></param>
        /// <returns></returns>
        public static double NormalizeDegreeRange(double angleD) => ((angleD % 360.0) + 360.0) % 360.0;

        /// <summary>
        /// Converts angle in radians (open) to radians (-2pi, +2pi)
        /// </summary>
        /// <param name="angleR"></param>
        /// <returns></returns>
        public static double NormalizeRadianRange(double angleR) => ((angleR % Celestial.TWOPI) + Celestial.TWOPI) % Celestial.TWOPI;

        /// <summary>
        /// Converts timespan hours (open) to hours (0,24)
        /// </summary>
        /// <param name="hours"></param>
        /// <returns></returns>
        public static double NormalizeHours(TimeSpan hours) => ((hours.TotalHours % 24.0) + 24.0) % 24.0;

        /// <summary>
        /// Converts hours (open) to hours (-24,+24)
        /// </summary>
        /// <param name="hours"></param>
        /// <returns></returns>
        public static double NormalizeHours(double hours) => ((hours % 24.0) + 24.0) % 24.0;


        #endregion
    }
}
