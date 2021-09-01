using System;
using System.Security.AccessControl;

namespace AstroMath
{
    public partial class Celestial
    {
        ////Module Celestial:  Functions for astronomical calculations
        #region Description
        //Format Conventions::
        //  date:                   datetime instance
        //  time:                   timespan instance
        //  days:                   double
        //  hours:                  double [0 < 24]
        //  centuries:              double 
        //  aradians:               double [0 < 2*PI]
        //  mradians:               double [-PI < +PI]
        //  degrees:                double [0 < 360.0]
        //
        //  RA:                     mradians
        //  Dec:                    mradians
        //  Alt:                    mradians
        //  Azm:                    mradians
        //  Hourangle:              hours
        //  Latitude:               mradians
        //  Longitude:              aradians
        //  Julian Date:            days
        //  Modified Julian Date:   days
        //  J2000 Date:             years    
        //  GMT:                    hours
        //  LST:                    hours
        //  
        //  LatLon                  object of terestrial location.              Properties: Latitude, Longitude
        //  AltAz                   object of horizontal celestial location.    Properties: Alt,      Az
        //  RADec                   object of equitorial celestial location.    Properties  RA,       Dec

        //Math and Trig Conversion Methods: 
        //
        //Planar.Frac(double)            -> double:      returns the fractional remainder of a number
        //SinD(degrees)                  -> double:      returns the sine of an angle in degrees
        //CosD(degrees)                  -> double:      returns the cosine of an angle in degrees
        //RadiansToDegrees(mradians)     -> degrees:     Converts radians to degrees
        //DegreesToRadians(degrees)      -> aradians:    Converts degrees to radians
        //HoursToRadians(hours)          -> aradians:    Converts hours to radians
        //Transform.RadiansToHours(mradians)       -> hours:       Converts radians to hours
        //HoursToDegrees(hours)          -> hours:       Converts hours to degrees
        //Transform.DegreesToHours(degrees)        -> degrees:     Converts degrees to hours

        //Date and Time Conversion Methods
        //
        //DateToJulian(date)             -> days         Converts from Julian Day to UTC date
        //DateToMJD(date)                -> days         Converts from UTC date to Modified Julian Day
        //DateToJ2kD(date)               -> days         Converts from UTC date to J2000 Day
        //DateToJ2kC(date)               -> centuries    Converts from UTC date to J2000 Century
        //JulianToDate(days)             -> date         Converts from Julian Day to UTC date
        //JulianToJ2kD(days)             -> date         Converts from Julian Day to J2000 Day
        //JulianToMJD(days)              -> date         Converts from Julian Day to Modified Julian Day
        //J2kDtoDate(days)               -> date         Converts from J2000 day to UTC date
        //J2kDToJulian(days)             -> days         Converts from J2000 day to Julian Day
        //J2kCToJulian(days)             -> days         Converts from J2000 century to Julian Day
        //MJDtoDate(days)                -> days         Converts from Modified Julian Days to UTC date
        //MJDToJulian(days)              -> days         Converts from Modified Julian Days to Julian Days 
        //TimeToString(hours)            -> string:      Constructs a string formated as "HH:mm" from hours
        //
        //Astronomical Time Conversion Methods
        //
        //DateToGST(date)                -> hours:       returns Greenich Sidereal Time as a function of UTC
        //GSTToLST(hours, LatLon)        -> hours:       returns Local Sidereal Time as a function of Greenich Sidereal Time at longitude 
        //J2KDToLMST(days, lon)          -> hours:       returns Local Sidereal Time as a function of J2000 days at longitude
        //
        //Astronomical Coordinate Conversion Methods
        //
        //HourAngleToRA(hours)           -> aradians:    returns Right Ascension as a function of hour angle at UTC and longitude
        //RAToHourAngle(radians,date,lon)-> hours:       returns Hour Angle as a function of Right Ascension at UTC and longitude
        //ComputeHA(radians,lat,altitude)-> hours:       returns Hour Angle as a function of Declination at altitude and latitude
        //
        //New LatLon(latitude,longitude) -> LatLon:      Creates LatLon object from Latitude (radians), Longitude (radians)
        //New RADec(RA,Dec)              -> RADec:       Creates RADec object from RA (radians), Dec (radians)
        //RADec.MakeAltAz(HA,location)   -> AltAz        Creates new AltAz object from Hour Angle (radians) and location (LatLon)
        //RADec.Altitude(HA,location)    -> aradians     returns altitude (radians) for RADec as a function of Hour Angle (radians) and location (LatLon)
        //RADec.Azimuth                  -> mradians     returns azimuth (mradians) for RADec as a function of Hour Angle (radians) and location (LatLon)
        //RADec.HourAngle                -> mradians     returns Hour Angle (mradians) for RADec as a function of UTC date and location (LatLon)
        //New AltAz(altitude, azimuth)   -> AltAz:       Creates AltAz object from Altitude (radians), Azimuth (radians)
        //AltAz.MakeRADec                -> RADec:       Creates new RADec object as a function of a Apoint object at a specific hourangle
        //AltAz.RightAscension(HA, loc)  -> mradians     returns Right Ascension for AltAz as a function of Hour Angle and location
        //AltAz.Declination(HA, loc)     -> aradians     returns Declination for AltAz as a function of Hour Angle and location

        //    //================================================================
        //    // Manifest constants
        //    //----------------------------------------------------------------

        #endregion

        #region Mathematical Constants

        public const double FIRST_GREGORIAN_YEAR = 1583.0;
        public const double JULIAN_BIAS = 2200000.0;
        public const double SIDEREAL_A = 0.0657098;

        public const double PI = Math.PI;
        public const double TWOPI = Math.PI * 2.0;
        public const double EPOCH2000 = 2451545.0;
        public const double EPOCHMJD = 2400000.5;

        public const double COSEPS = 0.91748;
        public const double SINEPS = 0.39778;
        public const double ARC = 206264.8062;
        #endregion

        #region RA Dec Class

        /// <summary>
        /// Class encapsulating methods and properties of a celestial position structure set by Right Ascension and Declination
        /// </summary>
        public class RADec
        {
            //Represents a sky location in equitorial coordinates of radians
            //
            //    .RA   [ Right Ascension in radians, in [0,+2*pi] ]
            //    .Dec  [ Dec in radians, in [-pi,+pi ]
            ////
            private double dr_RA;
            private double dr_Dec;

            //Empty constructor
            public RADec()
            {
                dr_RA = 0;
                dr_Dec = 0;
                return;
            }

            //RA and Dec parameterized constructor
            public RADec(double RA, double Dec)
            {
                dr_RA = Math.Abs(RA) % TWOPI;
                dr_Dec = Dec % Math.PI;
                return;
            }

            /// <summary>
            /// Right Ascension property
            /// </summary>
            public double RA
            {
                get => (dr_RA);
                set => dr_RA = Math.Abs(value) % TWOPI;
            }

            /// <summary>
            /// Declination property
            /// </summary>
            public double Dec

            {
                get => (dr_Dec);
                set => dr_Dec = value % Math.PI;
            }

            /// <summary>
            /// Create a new AltAz instance from an RADec instance, based on Hour Angle and terrestrial location (Lat/Long)
            /// </summary>
            /// <param name="haR">Hour Angle in radians</param>
            /// <param name="loc">Terretrial location in LatLon</param>
            /// <returns>Celestial location in AltAz instance</returns>
            public AltAz MakeAltAz(double haR, LatLon loc)
            {
                //Create a new AltAz instance from an RADec instance, i.e. convert from equitorial to horizontal.
                //  RADec.MakeAltAz -> AltAz
                //(ha is an object//s hour angle in radians) and
                //(latitude is the observer//s latitude in radians) ->
                //return position in the observer//s sky
                //      in horizon coordinates as an AltAz instance ]
                //   Compute altitude at declination given latitude and hourangle
                //  alt  =  altitude of object as seen from latLon at utc
                double alt = Altitude(haR, loc);
                double czm = Azimuth(haR, loc);
                if (Math.Sin(haR) < 0)
                { return new AltAz(alt, czm); }
                else
                { return new AltAz(alt, TWOPI - Math.Acos(czm)); }
            }

            /// <summary>
            /// Compute celestial altitude fraom RADec instance, based on Hour Ange and terrestrial location 
            /// </summary>
            /// <param name="haR">Hour Angle</param>
            /// <param name="location">LatLon instance</param>
            /// <returns>Altitude in radians</returns>
            public double Altitude(double haR, LatLon location)
            {
                //returns the altitude (radians) of the object at RADec for the given latitude and hour angle, in radians
                //(latitude is the observer//s latitude in radians) ->
                //   Compute altitude at declination given latitude and hourangle
                //  alt  =  altitude of object as seen from latLon at utc
                //
                double alt = Math.Asin((Math.Sin(dr_Dec) * Math.Sin(location.Lat)) + (Math.Cos(dr_Dec) * Math.Cos(location.Lat) * Math.Cos(haR)));
                return alt;
            }

            /// <summary>
            /// Compute celestial azimuth fraom RADec instance, based on Hour Ange and terrestrial location 
            /// </summary>
            /// <param name="haR">Hour Angle in radians</param>
            /// <param name="loc">Terrestrial location in LatLon</param>
            /// <returns>Azimuth in radians</returns>
            public double Azimuth(double haR, LatLon loc)
            {
                //Computes Azimuth (radians) of the RADec object for the given latitude and hour angle
                //  RADec.MakeAltAz -> AltAz
                //(ha is an object//s hour angle in radians (-pi,+pi) and
                //(latitude is the observer//s latitude in radians) ->
                //return position in the observer//s sky
                //  alt  =  altitude of object as seen from latLon at utc
                //double alt = Math.Asin(Math.Sin(dr_Dec) * Math.Sin(loc.Lat) + Math.Cos(dr_Dec) * Math.Cos(loc.Lat) * Math.Cos(haR));
                double alt = Altitude(haR, loc);
                //   az  :=  azimuth of object as seen from latLon at utc 
                double czm = Math.Acos((Math.Sin(dr_Dec) - Math.Sin(alt) * Math.Sin(loc.Lat)) / (Math.Cos(alt) * Math.Cos(loc.Lat)));
                //if (sin(HA) is negative,) { Azm = Azm, otherwise Azm = 2pi - Azm    
                if (Math.Sin(haR) < 0)
                { return czm; }
                else
                { return TWOPI - czm; }
            }

            /// <summary>
            /// Computes Hour Angle of RADec instance at UTC date/time and terrestrial location
            /// </summary>
            /// <param name="utcdate">UTC in datetime</param>
            /// <param name="location">Terrestrial location in LatLon</param>
            /// <returns>Hour Angle in radians</returns>
            public double HourAngle(DateTime utcdate, LatLon location)
            { //Calculate the hour angle(radians -- [-pi,+pi] ) for the current utc time & location (longitude -pi,+pi)
              //int haR = Celestial.HoursToRadians(Celestial.RAToHourAngle(r_RA, ldate, loc.Lon))
              //   Compute hourangle from RA in radians at Universal date & time and east longitude
                double lstH = GSTToLST(DateUTCToGST(utcdate), location.Lon);
                double lstR = Transform.HoursToRadians(lstH);
                double haR = (lstR - dr_RA);
                return haR;
            }

            /// <summary>
            /// Computes HourAngle for a given altitude, celelestial location (RADec) and Celestial Location (lat/lon)
            /// </summary>
            /// <param name="raDec"></param>
            /// <param name="location"></param>
            /// <returns>HourAngle in radians</returns>
            public double HourAngle(double altitude,  LatLon location)
            {
                //Added code 8/13/20 to deal with negative declinations
                double cosHA = (Math.Sin(altitude) - (Math.Sin(dr_Dec) * Math.Sin(location.Lat))) / (Math.Cos(dr_Dec) * Math.Cos(location.Lat));
                if (cosHA > 1) cosHA -= 1;
                if (cosHA < 1) cosHA += 1;
                return Math.Acos(cosHA);
                //return Math.Acos((Math.Sin(altitude) - (Math.Sin(declination) * Math.Sin(location.Lat))) / (Math.Cos(dr_Dec) * Math.Cos(location.Lat)));
            }

            /// <summary>
            /// Computes the local transit time for RADec instance on UTC date
            /// </summary>
            /// <param name="UTCDate">Date in datetime</param>
            /// <param name="location">Terrestrial location in LatLon</param>
            /// <returns>Transit Hour of day in decimal</returns>
            public double TransitTime(DateTime UTCDate, LatLon location)
            {
                //Compute the hour angle for the target at the current location
                //   subract the hour angle to the current time to get the transit time
                double dHA = Transform.RadiansToHours(HourAngle(UTCDate, location));
                DateTime localTransit = UTCDate.ToLocalTime() - TimeSpan.FromHours(dHA);
                double ttH = localTransit.Hour + (localTransit.Minute / 60.0);
                return ttH;
            }
        }
        #endregion

        #region Alt Az Class
        /// <summary>
        /// Class encapsulating methods and properties of a celestial position structure set by Azimuth and Altitude
        /// </summary>
        public class AltAz
        {
            //Represents a sky location in horizon coords. (altitude/azimuth)
            //
            //  Exports/Invariants
            //    .Alt   [ altitude in radians, in [-pi,+pi] ]
            //    .Azm   [ azimuth in radians, in [0,2*pi] ]
            ////
            private double r_alt;
            private double r_azm;

            //Empty constructor
            public AltAz()
            {
                r_alt = 0;
                r_azm = 0;
            }

            //Parameterized constructor
            public AltAz(double alt, double az)
            {
                r_alt = alt % Math.PI;
                r_azm = Math.Abs(az) % (TWOPI);
            }

            public double Alt
            {
                get => (r_alt);
                set => r_alt = value % Math.PI;
            }

            public double Azm
            {
                get => (r_azm);
                set => r_azm = Math.Abs(value) % (TWOPI);
            }

            /// <summary>
            /// Convert AltAz (from location at hourangle) to RaDec at terrestrial location and hour angle
            /// </summary>
            /// <param name="haR">Hour Angle</param>
            /// <param name="loc">Terrestrial location as Lat/Lon instance</param>
            /// <returns>Celstial location as Alt/Az instance</returns>
            public RADec MakeRaDec(double haR, LatLon loc)
            {
                //Convert between horizontal and equatorial coordinates.
                double RA = RightAscension(haR, loc);
                double Dec = Declination(haR, loc);
                return new RADec(RA, Dec);
            }

            /// <summary>
            /// Compute Right Ascension: form Alt/Az at terrestrial location nd hourangle)
            /// </summary>
            /// <param name="ha">Hour Angle in radians</param>
            /// <param name="loc">Terrestrial location as Lat/Long instance</param>
            /// <returns>Right Ascension in radians </returns>
            public double RightAscension(double ha, LatLon loc)
            {
                //Compute RA of object in horizontal coordinates at location LatLon and hour angle haR (radians)
                double RA = Math.Cos(r_alt) * Math.Sin(loc.Lat) * Math.Cos(r_azm) + Math.Sin(r_alt) * Math.Cos(r_azm);
                double Dec = Math.Cos(r_alt) * Math.Sin(r_azm);
                double z = Math.Sin(r_alt) * Math.Sin(loc.Lat) * Math.Sin(loc.Lat) * Math.Cos(r_azm);
                return RA;
            }

            //Declination: (from location at hourangle) to RaDec
            /// <summary>
            /// Compute Declination from Alt/Az instance at given Hour Angle and Terrestrial location
            /// </summary>
            /// <param name="ha">Hour Angle in radians</param>
            /// <param name="loc">Terrestrial location as LatLon instance</param>
            /// <returns>Declination in radians</returns>
            public double Declination(double ha, LatLon loc)
            {
                //Convert between horizontal and equatorial coordinates.
                double RA = Math.Cos(r_alt) * Math.Sin(loc.Lat) * Math.Cos(r_azm) + Math.Sin(r_alt) * Math.Cos(r_azm);
                double Dec = Math.Cos(r_alt) * Math.Sin(r_azm);
                double z = Math.Sin(r_alt) * Math.Sin(loc.Lat) * Math.Sin(loc.Lat) * Math.Cos(r_azm);
                return Dec;
            }
        }
        #endregion

        #region Date and Time Functions

        /// <summary>
        /// Converts civil date to Julian days (double)
        /// </summary>
        /// <param name="civilDate">Civil Date</param>
        /// <returns>Julian Days</returns>
        public static double DateToJulian(DateTime civilDate)
        {
            //   JD = 367K - <(7*(K+<(M+9)/12>))/4> + <(275M)/9> + I + 1721013.5 + UT/24 - 0.5sign(100K+M-190002.5) + 0.5
            double uth = civilDate.Hour + (civilDate.Minute / 60.0) + (civilDate.Second / 3600.0);
            double jd = (367.0 * civilDate.Year) - (int)(7.0 * (civilDate.Year + (int)((civilDate.Month + 9.0) / 12.0)) / 4.0) +
            (int)((275 * civilDate.Month) / 9) +
            civilDate.Day + 1721013.5 + (uth / 24.0) - (0.5 * Math.Sign(100.0 * civilDate.Year + civilDate.Month - 190002.5)) + 0.5;
            return jd;
        }

        /// <summary>
        /// Converts civil date to Modified Julian days (double)
        /// </summary>
        /// <param name="civilDate">Civil Date</param>
        /// <returns>Modified Julian Days</returns>
        public static double DateToMJD(DateTime civilDate)
        {
            //   Compute the Modified Julian Date from current Date/time
            double a;
            int b;

            int tyear = civilDate.Year;
            int tmonth = civilDate.Month;
            int tday = civilDate.Day;
            int thour = (int)((civilDate.Hour + civilDate.Minute / 60.0) + (civilDate.Second / 3600.0));


            a = 10000 * tyear + 100 * tmonth + tday;
            if (tmonth < 2)
            {
                tmonth = tmonth + 12;
                tyear = tyear - 1;
            }
            if (a <= 15821004.1)
            {
                b = -2 + (int)((tyear + 4716) / 4.0) - 1179;
            }
            else
            {
                b = (int)(tyear / 400.0) - (int)(tyear / 100.0) + (int)(tyear / 4.0);
            }
            a = 365 * tyear - 679004;
            double mjd = a + b + (int)(30.6001 * (tmonth + 1)) + tday + thour / 24.0;
            return mjd;
        }

        /// <summary>
        /// Convert civil date to J2000 days (double)
        /// </summary>
        /// <param name="civilDate">Civil Date</param>
        /// <returns>J2000 Days</returns>
        public static double DateToJ2kD(DateTime civilDate)
        {
            //   Convert a Date To J2000 days (Julian days from 1/1/2000, 0:0:0)
            //Compute the Julian days in the current epoch (2000)
            //   Convert date to julian,) { julian to J2000
            return JulianToJ2kD(DateToJulian(civilDate));
        }

        /// <summary>
        /// Converts civil date to J2000 centuries (double)
        /// </summary>
        /// <param name="civilDate">Civil Date</param>
        /// <returns>J2000 Centuries</returns>
        public static double DateToJ2kC(DateTime civilDate)
        {
            //   Convert a Date To J2000 centuries (Julian days from 1/1/2000, 0:0:0)
            //Compute the Julian days in the current epoch (2000)
            //   Convert date to julian,) { julian to J2000
            //int d = DateToJulian(thisdate) - epoch2000)
            //return JulianToJ2KC(DateToJulian(thisdate))
            ////return((DateToJulian(thisdate) - EPOCH2000) / 36525)
            return ((DateToJ2kD(civilDate)) / 36525.0);
        }

        //JULIANTODATE(JD days as datetime) -> Date/Time as datetime:
        /// <summary>
        /// Converts Julian days to a civil date
        /// </summary>
        /// <param name="julianDays">Julian Days</param>
        /// <returns>Civil Date</returns>
        public static DateTime JulianToDate(double julianDays)
        {
            //   Convert Julian Days To a civil date
            int b, d, f;
            double jd0, c, e;
            jd0 = (int)(julianDays + 0.5);
            if (jd0 < 2299161.0)
            {
                c = jd0 + 1524.0;
            }
            else
            {
                b = (int)((jd0 - 1867216.25) / 36524.25);
                c = jd0 + (b - (int)(b / 4.0)) + 1525.0;
            }
            d = (int)((c - 122.1) / 365.25);
            e = 365.0 * d + (int)(d / 4.0);
            f = (int)((c - e) / 30.6001);
            int day = (int)(c - e + 0.5) - (int)(30.6001 * f);
            int month = (int)(f - 1.0 - 12.0 * (int)(f / 14.0));
            int year = d - 4715 - (int)((7.0 + month) / 10.0);
            DateTime dateout = new DateTime(year, month, day);
            dateout = dateout.AddHours(24.0 * (julianDays + 0.5 - jd0));
            return dateout;
        }

        //JULIANTOJ2kD(Julian days as double) -> Time in Julian days since 2000
        /// <summary>
        /// Converts Julian days to J2000 days
        /// </summary>
        /// <param name="julianDays">Julian Days</param>
        /// <returns>J2000 Days</returns>
        public static double JulianToJ2kD(double julianDays)
        {
            //   Convert Julian Days To J2000 days
            return (julianDays - EPOCH2000);
        }

        /// <summary>
        /// Converts Julian days to Julian Centuries
        /// </summary>
        /// <param name="julianDays">Julian Days</param>
        /// <returns>Julian Centuries</returns>
        public static double JulianToJ2KC(double julianDays)
        {
            //   Convert Julian Days To J2000 centuries
            return ((julianDays - EPOCH2000) / 36525.0);
        }

        //JULIANTOMJD(Julian days as double) -> Modified Julian Days as double
        /// <summary>
        /// Converts Julian days to Modified Julian days
        /// </summary>
        /// <param name="julianDays">Julian Days</param>
        /// <returns>Modified Julian Days</returns>
        public static double JulianToMJD(double julianDays)
        {
            //   Convert Modified Julian Days to Julian days
            return (julianDays - EPOCHMJD);
        }

        /// <summary>
        /// Converts J2000 days to civil date
        /// </summary>
        /// <param name="j2k">J2000 Days</param>
        /// <returns>Civil Date</returns>
        public static DateTime J2kDToDate(double j2k)
        {
            //   Convert J2000 Days To a civil date
            int b, d, f;
            double jd, jd0, c, e;


            jd = J2kDToJulian(j2k);
            jd0 = (int)(jd + 0.5);
            if (jd0 < 2299161.0)
            {
                c = jd0 + 1524.0;
            }
            else
            {
                b = (int)((jd0 - 1867216.25) / 36524.25);
                c = jd0 + (b - (int)(b / 4.0)) + 1525.0;
            }
            d = (int)((c - 122.1) / 365.25);
            e = 365.0 * d + (int)(d / 4);
            f = (int)((c - e) / 30.6001);
            int day = (int)(c - e + 0.5) - (int)(30.6001 * f);
            int month = (int)(f - 1.0 - 12.0 * (int)(f / 14.0));
            int year = d - 4715 - (int)((7.0 + month) / 10.0);
            DateTime jdate = new DateTime(year, month, day);
            jdate = jdate.AddHours(24.0 * (jd + 0.5 - jd0));
            return jdate;
        }

        /// <summary>
        /// Converts J2000 days to Julian days
        /// </summary>
        /// <param name="j2k">J2000 Days</param>
        /// <returns>Julian Days</returns>
        public static double J2kDToJulian(double j2k)
        {
            //   Convert J2000 days To Julian Days
            return (EPOCH2000 + j2k);
        }

        /// <summary>
        /// Converts J2000 centuries to Jcivil date
        /// </summary>
        /// <param name="j2kc">J2000 Centuries</param>
        /// <returns>Civil Date</returns>
        public static DateTime J2kCToDate(double j2kc)
        {
            //   Convert J2000 years To Gregorian date
            return (JulianToDate(J2kCToJulian(j2kc)));
        }

        /// <summary>
        /// Converts J2000 centuries to Julian days
        /// </summary>
        /// <param name="j2kc">Julian Centuries</param>
        /// <returns>Julian Days</returns>
        public static double J2kCToJulian(double j2kc)
        {
            //   Convert J2000 years To Julian Days
            return (EPOCH2000 + (j2kc * 36525.0));
        }

        /// <summary>
        /// Converts Modified Julian days to civil date
        /// </summary>
        /// <param name="mjd">MOdified Julian Days</param>
        /// <returns>Civil Date</returns>
        public static DateTime MJDToDate(double mjd)
        {
            //   Convert Modified Julian Days To a civil date
            int b, d, f;
            double jd, jd0, c, e;
            jd = mjd + EPOCHMJD;
            jd0 = (int)(jd + 0.5);
            if (jd0 < 2299161.0)
            {
                c = jd0 + 1524.0;
            }
            else
            {
                b = (int)((jd0 - 1867216.25) / 36524.25);
                c = jd0 + (b - (int)(b / 4)) + 1525.0;
            }
            d = (int)((c - 122.1) / 365.25);
            e = 365.0 * d + (int)(d / 4.0);
            f = (int)((c - e) / 30.6001);
            int day = (int)(c - e + 0.5) - (int)(30.6001 * f);
            int month = f - 1 - 12 * (int)(f / 14.0);
            int year = d - 4715 - (int)((7 + month) / 10.0);
            DateTime dateout = new DateTime(year, month, day);
            dateout = dateout.AddHours(24.0 * (jd + 0.5 - jd0));
            return dateout;
        }

        /// <summary>
        /// Converts Modified Julian days to Julian days
        /// </summary>
        /// <param name="mjd">MOdified Julian Days</param>
        /// <returns>Julian Days</returns>
        public static double MJDToJulian(double mjd)
        {
            //   Convert Modified Julian Days to Julian days
            return (mjd + EPOCHMJD);
        }

        /// <summary>
        /// Converts Modified Julian days to J2000 days
        /// </summary>
        /// <param name="mjd">Modified Julian Days</param>
        /// <returns>J2000 days</returns>
        public static double MJDToJ2kD(double mjd)
        {
            //   Convert MJD Days To J2000 Days
            return JulianToJ2kD(mjd + EPOCHMJD);
        }

        /// <summary>
        /// Converts J2000 days at given longitude to LMST hours
        /// </summary>
        /// <param name="j2k">J2000 Days</param>
        /// <param name="longR">Longitude</param>
        /// <returns>LMST Hours</returns>
        public static double J2kDToLMST(double j2k, double longR)
        {
            //   Convert J2000 Days at longitude to Local Mean Sidereal Time in hours
            return GSTToLST(DateUTCToGST(J2kDToDate(j2k)), longR);
        }

        /// <summary>
        /// Converts Modified Julian days at given longitude to LMST hours
        /// </summary>
        /// <param name="mjd">Modified Julian days</param>
        /// <param name="longR">Longitude</param>
        /// <returns>LMST Hours</returns>
        public static double MJDtoLMST(double mjd, double longR)
        {
            //   Convert Modified Julian Days at longitude to Local Mean Sidereal Time in hours
            double mjd0, t, ut, gmst;
            mjd0 = (int)(mjd);
            ut = (mjd - mjd0) * 24.0;
            t = (mjd0 - 51544.5) / 36525.0;
            gmst = 6.697374558 + 1.0027379093 * ut + (8640184.812866 + (0.093104 - 0.0000062 * t) * t) * t / 360.00;
            double lmst = 24.0 * Planar.Frac((gmst - Transform.RadiansToHours(longR)) / 24.0);
            return lmst;
        }

        /// <summary>
        /// Converts UTC date to GST hours
        /// </summary>
        /// <param name="utcDate">UTC date/time</param>
        /// <returns>Greenwich Sidereal Time in hours</returns>
        public static double DateUTCToGST(DateTime utcDate)
        {
            //   Compute Greenich Sidereal Time in hours from user date & time
            // Greenwich Sidereal Time is translated from UTC from the number of days from J2000.
            double hours = Transform.DegreesToHours((280.6061837 + 360.98564736629 * DateToJ2kD(utcDate)) % 360.0);
            return hours;
        }

        /// <summary>
        /// Converts UTC date/time at given longitude to local time in hours
        /// </summary>
        /// <param name="lst">Local Sidereal Time</param>
        /// <param name="longitudeD">longitude</param>
        /// <returns>Local time in hours</returns>
        public static double LSTToLocalTimeHours(double lst, double longitudeD)
        {
            double JD0 = (int)(DateToJulian(DateTime.UtcNow)) + 0.5;
            double gt = 6.656306 + 0.0657098242 * (JD0 - 2445700.5) + 1.0027379093 * (DateTime.UtcNow.Hour + DateTime.UtcNow.Minute / 60.0);
            double gmst = Transform.NormalizeHours(lst - (longitudeD / 15.0));
            double ut = Planar.Frac(((gmst - 6.697374558) - (0.0657098242 * (JD0 - 2445700.5))) / 1.0027379093) * 24;
            if (ut < 0)
            { ut += 24; }
            TimeSpan ltz = DateTime.Now - DateTime.UtcNow;
            double lt = Transform.NormalizeHours(ut + ltz.TotalHours);
            return lt;
        }

        /// <summary>
        /// Converts GST hours at given longitude to LST hours
        /// </summary>
        /// <param name="gst">Greenwich Sidereal Time</param>
        /// <param name="longitudeR">Longitude</param>
        /// <returns>Local Sidereal Time</returns>
        public static double GSTToLST(double gst, double longitudeR)
        {
            //   Compute Local Sideral Time in hours from Greenich Sidereal Time in hours at a longitude in radians
            //Local Sidereal Time is Greenwich Sidereal Time decremented by 1 hour per 15 degrees of site longitude, independent of date
            //gsth is GST in hours
            //site is longitude in radians -- positive for east longitude (-pi,+pi)
            //double lst = ((gst + Transform.RadiansToHours(longitudeR)) + 24) % 24;
            double lst = gst + Transform.RadiansToHours(longitudeR);
            return lst;
        }

        /// <summary>
        /// Convert LST hours at given longitude to GST hours
        /// </summary>
        /// <param name="lst">Local Sidereal Time</param>
        /// <param name="location">Latitude/Longitude</param>
        /// <returns>Greenwich Sidereal Time in hours</returns>
        public static double LSTToGST(double lst, LatLon location)
        {
            //   Compute Local Sideral Time in hours from Greenich Sidereal Time in hours at a longitude in radians
            //Local Sidereal Time is Greenwich Sidereal Time decremented by 1 hour per 15 degrees of site longitude, independent of date
            //gsth is GST in hours
            //site is longitude in radians -- positive for east longitude
            double gst = (lst + Transform.RadiansToHours(location.Lon) + 24) % 24.0;
            return gst;
        }
        #endregion

        #region Time Comparison Methods

        /// <summary>
        /// Combines a date with hours to make a date/time
        /// </summary>
        /// <param name="thisdate">Civil Date</param>
        /// <param name="someHours">Hours</param>
        /// <returns>Combined date plus hours</returns>
        public static DateTime DayPlusHours(DateTime thisdate, double someHours)
        {
            //  return a datetime for this date at given hour(s)
            DateTime dateday = new DateTime(thisdate.Year, thisdate.Month, thisdate.Day).AddHours(someHours);
            return dateday;
        }

        //TimeInBetween(earliestTime as datetime, latestTime as datetime, thisTime as datetime) as boolean
        /// <summary>
        /// Determines if (this time is later than the earliestTime but sooner than the latestTime, ignoring the date
        /// </summary>
        /// <param name="earliestTime">Starting date/time</param>
        /// <param name="latestTime">Ending date/time</param>
        /// <param name="thisTime">Test date/time</param>
        /// <returns>True if Test is after start and before end, false otherwise</returns>
        public static bool TimeInBetween(DateTime earliestTime, DateTime latestTime, DateTime thisTime)
        {
            if ((earliestTime <= thisTime) && (thisTime <= latestTime))
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        /// <summary>
        /// Computes the hours between rise and set that overlap the hours between dusk and dawn where
        ///    all datetimes are normalized to the same day.
        /// </summary>
        /// <param name="iDusk">Start of night</param>
        /// <param name="iDawn">End of night</param>
        /// <param name="iRise">Start of rise</param>
        /// <param name="iSet">End of rise</param>
        /// <returns>Overlap of intervals in hours</returns>
        public static double IntervalOverlap(DateTime iDusk, DateTime iDawn, DateTime iRise, DateTime iSet)
        {
            //
            //
            //Break the day into two parts - Dusk to midnight and midnight to dawn

            double iHours = 0;

            DateTime Dusk = new DateTime(1, 1, 1, iDusk.Hour, iDusk.Minute, iDusk.Second);
            DateTime Dawn = new DateTime(1, 1, 1, iDawn.Hour, iDawn.Minute, iDawn.Second);
            DateTime Rise = new DateTime(1, 1, 1, iRise.Hour, iRise.Minute, iRise.Second);
            DateTime Set = new DateTime(1, 1, 1, iSet.Hour, iSet.Minute, iSet.Second);
            DateTime BeforeMidnight = new DateTime(1, 1, 1, 23, 59, 59);
            DateTime AfterMidnight = new DateTime(1, 1, 1, 0, 0, 0);
            //TSX returns 12 and 12 for objects that never set.  if so, bump the rise time by a minute just for this algorithm
            //variables used for debug
            double nightDuration = Transform.NormalizeHours(Dawn - Dusk);
            double objectUpDuration = Transform.NormalizeHours(Set - Rise);

            //Dusk period - before midnight
            if (Rise == Set)
            { iHours = nightDuration; }
            else
            {
                //Dusk Period
                if (Rise < Set) //Rise before Set
                {
                    if (Rise < Dusk) //Rise Before Dusk
                    {
                        if (Set < Dusk)
                        //Set before dusk //never up after dusk 
                        { iHours = 0; }
                        else
                        //Set after dusk, must be before or at midnight //dusk to set     
                        { iHours = (Set - Dusk).TotalHours; }
                    }
                    else
                    //Rise after Dusk, set must be before or at midnght//rise to set
                    { iHours = (Set - Rise).TotalHours; }
                }
                else
                // set before rise
                {
                    if (Set < Dusk) //Set before dusk
                    {
                        if (Rise < Dusk) //rise before dusk//dusk to midnight
                        { iHours = (BeforeMidnight - Dusk).TotalHours; }
                        else//rise after dusk up to midnight //rise to midnight 
                        { iHours = (BeforeMidnight - Rise).TotalHours; }
                    }
                    else//set after dusk, rise must be between set and midnight//dusk to set + rise to midnight
                    { iHours = (Set - Dusk).TotalHours + (BeforeMidnight - Rise).TotalHours; }
                }
                //Dawn period
                if (Rise < Set)//Rise before set
                {
                    if (Rise < Dawn)//rise before dawn
                    {
                        if (Set < Dawn)//set before dawn //rise to set
                        { iHours += (Set - Rise).TotalHours; }
                        else //Set after dawn//Rise to dawn
                        { iHours += (Dawn - Rise).TotalHours; }
                    }
                    else//rise after dawn//never up
                    { iHours += 0; }
                }
                else //Set before Rise
                {
                    if (Set < Dawn) //set before dawn
                    {
                        if (Rise < Dawn) //rise before dawn//midnight to set + rise to dawn
                        { iHours += ((Set - AfterMidnight).TotalHours + (Dawn - Rise).TotalHours); }
                        else //rise after dawn//midnight to set
                        { iHours += (Set - AfterMidnight).TotalHours; }
                    }
                    else //set after dawn//midnight to dawn
                    { iHours += (Dawn - AfterMidnight).TotalHours; }
                }
            }
            return iHours;
        }

        /// <summary>
        /// Classifies the type of intersection between a given date/time interval and three other date/time intervals
        ///     (used to determine target, sun and moon positions over three day period
        /// </summary>
        /// <param name="u">Start of test interval</param>
        /// <param name="d">Emd pf test interval</param>
        /// <param name="r1">Start of first period</param>
        /// <param name="s1">End of first period</param>
        /// <param name="r2">Start of second period</param>
        /// <param name="s2">End of second period</param>
        /// <param name="r3">Start of third period</param>
        /// <param name="s3">End of third period</param>
        /// <returns>One of 28 possible intersections, 0 if indeterminate</returns>
        public static int TimeMachine(DateTime u, DateTime d, DateTime r1, DateTime s1, DateTime r2, DateTime s2, DateTime r3, DateTime s3)
        {
            //  return a value for each type of intersection, return 0 if something doesn//t match up
            int stype = 0;
            if (u < r1)
            {
                if (d < r1) { stype = 1; }
                else if (r1 <= d && d <= s1) { stype = 2; }
                else if (s1 <= d && d <= r2) { stype = 3; }
                else if (r2 <= d && d <= s2) { stype = 4; }
                else if (s2 <= d && d <= r3) { stype = 5; }
                else if (r3 <= d && d <= s3) { stype = 6; }
                else if (d > s3) { stype = 7; }
            }
            else if (r1 <= u && u <= s1)
            {
                if (u <= d && d <= s1) { stype = 8; }
                else if (s1 <= d && d <= r2) { stype = 9; }
                else if (r2 <= d && d <= s2) { stype = 10; }
                else if (s2 <= d && d <= r3) { stype = 11; }
                else if (r3 <= d && d <= s3) { stype = 12; }
                else if (d > s3) { stype = 13; }
            }
            else if (s1 <= u && u <= r2)
            {
                if (u <= d && d <= r2) { stype = 14; }
                else if (r2 <= d && d <= s2) { stype = 15; }
                else if (s2 <= d && d <= r3) { stype = 16; }
                else if (r3 <= d && d <= s3) { stype = 17; }
                else if (d > s3) { stype = 18; }
            }
            else if (r2 <= u && u <= s2)
            {
                if (u <= d && d <= s2) { stype = 19; }
                else if (s2 <= d && d <= r3) { stype = 20; }
                else if (r3 <= d && d <= s3) { stype = 21; }
                else if (d > s3) { stype = 22; }
            }
            else if (s2 <= u && u <= r3)
            {
                if (u <= d && d <= r3) { stype = 23; }
                else if (r3 <= d && d <= s3) { stype = 24; }
                else if (d > s3) { stype = 25; }
            }
            else if (r3 <= u && u <= s3)
            {
                if (u <= d && d <= s3) { stype = 26; }
                else if (d > s3) { stype = 27; }
            }
            else if (s3 < u) { stype = 28; }
            return stype;
        }

        /// <summary>
        /// Determines which of three time spans is the longest
        /// </summary>
        /// <param name="a">First time span</param>
        /// <param name="b">Second time span</param>
        /// <param name="c">Third time span</param>
        /// <returns>A value: 1, 2 or 3, whichever is longest </returns>
        public static int LongestPeriod(TimeSpan a, TimeSpan b, TimeSpan c)
        {
            //returns a value, 1, 2, or 3 according to which parameter is the longest timespan
            if (a >= b && a >= c) { return 1; }
            else if (b >= a && b >= c) { return 2; }
            else { return 3; }
        }

        /// <summary>
        /// Determines which of two time spans it the longest
        /// </summary>
        /// <param name="i1">first time span</param>
        /// <param name="i2">Second time span</param>
        /// <returns>Either first or second time span, whichever is longest</returns>
        public static TimeSpan LongestInterval(TimeSpan i1, TimeSpan i2)
        {
            //return the longest of two timespans
            if (i1.TotalHours > i2.TotalHours) { return i1; }
            else { return i2; }
        }
        #endregion

        #region Astronomical Coordinates Conversion Methods

        /// <summary>
        /// Converts an Hour Angle at a given UTC and longitude to a Right Ascension 
        /// </summary>
        /// <param name="ha">Hour Angle</param>
        /// <param name="ut">UTC date/time</param>
        /// <param name="longitude">Longitude</param>
        /// <returns>Right Ascension in radians</returns>
        public static double HourAngleToRA(double ha, DateTime ut, double longitude)
        {
            //  Compute hourangle at Universal Time and longitude
            double lst = GSTToLST(DateUTCToGST(ut), longitude);
            double ra = (lst - ha) % TWOPI;
            return ra;
        }

        #endregion

        #region Astronomical Operations Methods

        /// <summary>
        /// Computes a RA/Dec position which is a given bearing and range (radians) from an initial position (RA,Dec)
        /// </summary>
        /// <param name="initialPosition"></param>
        /// <param name="bearingRadians"></param>
        /// <param name="rangeRadians"></param>
        /// <returns></returns>
        public static Celestial.RADec ComputePositionFromBearingAndRange(Celestial.RADec initialPosition, double bearingRadians, double rangeRadians)
        {
            /* Calulates a new RA/Dec position that is a specific distance and direction from the initial position

            * lat2 = math.asin( math.sin(lat1)*math.cos(d/R) + 
            *                   math.cos(lat1)*math.sin(d/R)*math.cos(brng))
            * lon2 = lon1 + math.atan2(math.sin(brng)*math.sin(d/R)*math.cos(lat1),
            *                          math.cos(d/R)-math.sin(lat1)*math.sin(lat2))
           */
            Celestial.RADec endPosition = new Celestial.RADec();
            endPosition.Dec = Math.Asin(Math.Sin(initialPosition.Dec) * Math.Cos(rangeRadians) +
                                        Math.Cos(initialPosition.Dec) * Math.Sin(rangeRadians) * Math.Cos(bearingRadians));
            endPosition.RA = initialPosition.RA + Math.Atan2(Math.Sin(bearingRadians) * Math.Sin(rangeRadians) * Math.Cos(initialPosition.Dec),
                                                Math.Cos(rangeRadians) - Math.Sin(initialPosition.Dec) * Math.Sin(endPosition.Dec));
            return endPosition;
        }

        #endregion

        #region Lat Long Class

        /// <summary>
        /// Class containing structure, methods and properties for definine terrestrial location in terms of latitude and longitude
        /// </summary>
        public class LatLon
        {
            //Represents a terrestrial location in latitude, longitude
            //
            //  properties
            //    .lat   [ latitude in radians, in [-pi,+pi] ]
            //    .lon    [ longitude in radians, in [0,+2*pi] ]

            private double r_lat;
            private double r_lat_deg;
            private string r_lat_dir;
            private double r_lon;
            private double r_lon_deg;
            private string r_lon_dir;

            /// <summary>
            /// Empty LatLon constructor
            /// </summary>
            public LatLon()
            {
                r_lat = 0;
                r_lat_deg = 0;
                r_lat_dir = "N";
                r_lon = 0;
                r_lon_deg = 0;
                r_lon_dir = "E";
                return;
            }

            /// <summary>
            /// Parameterized constructor -- -pi <= lat <= +pi;  -pi <= lon <= +pi
            /// </summary>
            /// <param name="latitude">Latitude in radians</param>
            /// <param name="longitude">Longitude in radians</param>
            public LatLon(double latitude, double longitude)
            {
                r_lat = latitude % Math.PI;
                r_lon = longitude % TWOPI;

                if (r_lat < 0)
                {
                    r_lat_dir = "S";
                    r_lat_deg = -(360.0 * r_lat / TWOPI);
                }
                else
                {
                    r_lat_dir = "N";
                    r_lat_deg = (360.0 * r_lat / TWOPI);
                }

                if (longitude >= Math.PI)
                {
                    r_lon_dir = "E";
                    r_lon_deg = Transform.RadiansToDegrees(Math.PI - r_lon);
                }
                else
                {
                    r_lon_dir = "W";
                    r_lon_deg = Transform.RadiansToDegrees(r_lon);
                }
            }

            public double Lat => r_lat;
            public double Lon => r_lon;

            public string GetLatitudeString()
            { return (r_lat_deg.ToString() + " " + r_lat_dir); }

            public string GetLongitudeString()
            { return (r_lon_deg.ToString() + " " + r_lon_dir); }
        }
        #endregion

    }
}
