using System;

namespace AstroMath

{
    public class DailyPosition
    {

        //this module contains methods for calculating the daily positions of
        //  astronomical targets based on date, RA and Dec, including sun and moon

        //R. McAlister, V1.0, 12/5/16
        //

        public enum VisibilityState
        {
            UpSome,
            UpAlways,
            Rises,
            Falls,
            DownSome,
            UpNever
        }

        private DateTime t_utcdate;
        private Celestial.RADec t_position;
        private Celestial.LatLon t_location;
        private double t_minalt;
        private VisibilityState t_state;
        private DateTime t_rise;
        private DateTime t_set;
        private DateTime i_rise;//initial start date
        private DateTime i_set;  //initaial end date
        private double t_moonfree;
        private double t_moonphase;

        public DailyPosition()
        {
            //Empty constructor
            t_state = VisibilityState.UpSome;
            //t_utcdate = null;
            //t_rise = null;
            //t_set = null;
            t_moonfree = 0;
            t_moonphase = 0;
            return;
        }

        public DailyPosition(DateTime utcdateStart, DateTime utcdateEnd, Celestial.RADec tradec, Celestial.LatLon tlatlon, double tminAlt)
        {
            //utcdate as utc date (datetime)
            //tradec as position of tar{ get { (RADec)
            //tlatlon as location of observer (LatLon)
            //minAlt as minimum horizon (degrees):  negative is below horizon, positive is above horizon

            bool AboveToStart;
            bool rise;
            bool sett;
            t_state = VisibilityState.UpSome;

            t_moonfree = 0;
            t_moonphase = 0;
            double thour;
            double utriseH = 0;
            double utsetH = 0;
            double yminus, yzero, yplus;

            Planar.QuadRoot solve;

            double sinMinAltD = Transform.SinD(tminAlt);   //Rise at h = 0 arcminutes 
            DateTime udate = utcdateStart;
            rise = false;
            sett = false;
            AboveToStart = false;
            thour = 1;
            double altitude = tradec.Altitude(tradec.HourAngle(udate.AddHours(thour - 1), tlatlon), tlatlon);
            yminus = Math.Sin(tradec.Altitude(tradec.HourAngle(udate.AddHours(thour - 1), tlatlon), tlatlon)) - sinMinAltD;
            if (yminus > 0)
            {
                AboveToStart = true;
            }
            do
            {
                yzero = Math.Sin(tradec.Altitude(tradec.HourAngle(udate.AddHours(thour), tlatlon), tlatlon)) - sinMinAltD;
                yplus = Math.Sin(tradec.Altitude(tradec.HourAngle(udate.AddHours(thour + 1), tlatlon), tlatlon)) - sinMinAltD;

                solve = Planar.Quad(yminus, yzero, yplus);

                switch (solve.nz)
                {
                    case 0:
                        break;
                    case 1:
                        if (yminus < 0)
                        {
                            utriseH = thour + solve.zero1;
                            rise = true;
                        }
                        else
                        {
                            utsetH = thour + solve.zero1;
                            sett = true;
                        }
                        break;
                    case 2:
                        if (solve.ye < 0)
                        {
                            utriseH = thour + solve.zero2;
                            utsetH = thour + solve.zero1;
                        }
                        else
                        {
                            utriseH = thour + solve.zero1;
                            utsetH = thour + solve.zero2;
                        }
                        rise = true;
                        sett = true;
                        break;
                }
                yminus = yplus;
                thour = thour + 2;
            }
            while (!((thour == 25) || (rise && sett) || (utcdateStart.AddHours(thour + 1) > utcdateEnd)));

            t_utcdate = utcdateStart;
            t_position = tradec;
            t_location = tlatlon;
            t_minalt = tminAlt;

            i_rise = utcdateStart;
            i_set = utcdateEnd;

            //Condition abovetostart and rise or not abovetostart and sett cannot occur (at least on paper)
            if (AboveToStart && (!(rise || sett)))
            {
                //Tar{ get { path is always above minimum altitude (e.g. horizon)
                t_state = VisibilityState.UpAlways;
                t_rise = utcdateStart;
                t_set = utcdateEnd;
            }
            else if ((!AboveToStart) && (rise && (!sett)))
            {
                //Tar{ get { path starts below) { ascends above minimum altitude (e.g. horizon)
                t_state = VisibilityState.Rises;
                t_rise = udate.AddHours(utriseH);
                t_set = utcdateEnd;
                //if (t_rise > t_set) {
                //    t_set = t_set.AddDays(1)
                //}
            }
            else if (AboveToStart && ((!rise) && sett))
            {
                //Tar{ get { path starts above) { decends below minimum altitude (e.g. horizon)
                t_state = VisibilityState.Falls;
                t_rise = utcdateStart;
                t_set = udate.AddHours(utsetH);
                if (t_rise > t_set)
                {
                    t_set = t_set.AddDays(1);
                }
            }
            else if (AboveToStart && (rise && sett))
            {
                //Tar{ get { path decends below) { rises above minimum altitude (e.g. horizon)
                //Choose the longer of the two rise/set intervals 
                t_state = VisibilityState.DownSome;
                t_rise = udate.AddHours(utriseH);
                t_set = udate.AddHours(utsetH);
                if ((t_set - i_rise) > (i_set - t_rise))
                {
                    t_rise = i_rise;
                }
                else
                {
                    t_set = i_set;
                }
                //rise should be after set
                //if (t_rise < t_set) {
                //    t_rise += TimeSpan.FromDays(1)
                //}
            }
            else
            {
                if ((!AboveToStart) && (rise && sett))
                {
                    //Tar{ get { path rises above) { decends below minimum altitude (e.g. horizon)
                    //Save the 
                    t_state = VisibilityState.UpSome;
                    t_rise = udate.AddHours(utriseH);
                    t_set = udate.AddHours(utsetH);
                    //set should be after rise
                    if (t_rise > t_set)
                    {
                        t_set += TimeSpan.FromDays(1);
                    }
                }
                else
                { //not above at all
                  //Tar{ get { path is always below minimum altitude (e.g. horizon)
                    t_state = VisibilityState.UpNever;
                    t_rise = utcdateStart;
                    t_set = utcdateStart;
                }

                //The rise time should not precede the beginning of the interval...
                if (t_rise < i_rise)
                {
                    t_rise = i_rise;
                }
                //The set time should not exceed the end of the interval...
                if (t_set > i_set)
                {
                    t_set = i_set;
                }
                return;
            }
        }

        public VisibilityState Visibility => (t_state);

        public DateTime UTCdate
        {
            get => (t_utcdate);
            set => t_utcdate = value;
        }

        public DateTime IntervalStartDate
        {
            get => (i_rise);
            set => i_rise = value;
        }

        public DateTime IntervalEndDate
        {
            get => (i_set);
            set => i_set = value;
        }

        public DateTime Rising
        {
            get => (t_rise);
            set => t_rise = value;
        }

        public DateTime Setting
        {
            get => (t_set);
            set => t_set = value;
        }

        public double MoonFree
        {
            get => (t_moonfree);
            set => t_moonfree = value;
        }

        public Celestial.RADec Position
        {
            get => (t_position);
            set => t_position = value;
        }

        public Celestial.LatLon Location
        {
            get => (t_location);
            set => t_location = value;
        }

        public double MinAlt
        {
            get => (t_minalt);
            set => t_minalt = value;
        }

        public double MoonPhase
        {
            get => (t_moonphase);
            set => t_moonphase = value;
        }

        public DateTime iRise
        {
            get => (i_rise);
            set => i_rise = value;
        }

        public DateTime iSet
        {
            get => (i_set);
            set => i_set = value;
        }

        //MoonPhase:  Compute approximate phase of moon as fraction based on sun ra (radians) and moon ra (radians)
        public void SetMoonPhase(double sunRA, double moonRA)
        {
            //moonazm in radians
            //sunazm in radians
            //The phase of the moon varies from 0 to 100% as sun - moon RA = +/-180 degrees (pi radians)
            //Note:  this is a bit of an approximation, except around the new and full moons
            //
            double dif = 1 - (Math.Abs(Math.PI - (Math.Abs(sunRA - moonRA)))) / Math.PI;
            t_moonphase = dif;
            return;
        }

        //public bool IsUp(DateTime positionTime)
        //{
        //    //Method returns true if (this object daily postion is above the horizon at the 
        //    //  positionTime or false if (it is not
        //    //if (()) {
        //    return true;
        //}

        #region Astronomical Methods

        public static Celestial.RADec SunRADec(double jc2K)
        {
            //jc2k = julian centuries since year 2000
            double t = jc2K;
            double m = Celestial.TWOPI * Planar.Frac(0.993133 + 99.997361 * t);
            double dl = 6893.0 * Math.Sin(m) + 72.0 * Math.Sin(2 * m);
            double l = Celestial.TWOPI * Planar.Frac(0.7859453 + m / Celestial.TWOPI + (6191.2 * t + dl) / 1296000.0);
            double sl = Math.Sin(l);
            double x = Math.Cos(l);
            double y = Celestial.COSEPS * sl;
            double z = Celestial.SINEPS * sl;
            double rho = Math.Sqrt(1.0 - z * z);
            double dec = (360.0 / Celestial.TWOPI) * Math.Atan(z / rho);
            double ra = (48.0 / Celestial.TWOPI) * Math.Atan(y / (x + rho));
            if (ra < 0)
            { ra = ra + 24.0; }
            Celestial.RADec a = new Celestial.RADec(Transform.HoursToRadians(ra), Transform.DegreesToRadians(dec));
            return a;
        }

        public static Celestial.RADec MoonRaDec(double j2kC)
        {
            //Produces approximate RA,Dec position of moon 
            //j2kC in julian centuries since the year 2000
            double lo = Planar.Frac(0.606433 + 1336.855225 * j2kC);     //mean longitude moon (in rev)
            double l = Celestial.TWOPI * Planar.Frac(0.374897 + 1325.55241 * j2kC);   //mean anomaly of the moon
            double ls = Celestial.TWOPI * Planar.Frac(0.993133 + 99.997361 * j2kC);   //mean anomally of the sun
            double d = Celestial.TWOPI * Planar.Frac(0.827361 + 1236.853086 * j2kC);    //diff longitude of Moon-Sun  
            double f = Celestial.TWOPI * Planar.Frac(0.259086 + 1342.227825 * j2kC);    //mean argument of latitude

            double dl = 22640 * Math.Sin(l) - 4586.0 * Math.Sin(l - 2.0 * d) + 2370.0 * Math.Sin(2.0 * d) + 769.0 * Math.Sin(2.0 * l) -
                668.0 * Math.Sin(ls) - 412.0 * Math.Sin(2.0 * f) - 212.0 * Math.Sin(2.0 * l - 2.0 * d) - 206.0 * Math.Sin(l + ls - 2.0 * d) +
                192.0 * Math.Sin(l + 2.0 * d) - 165.0 * Math.Sin(ls - 2.0 * d) - 125.0 * Math.Sin(d) - 110.0 * Math.Sin(l + ls) +
                148.0 * Math.Sin(l - ls) - 55.0 * Math.Sin(2.0 * f - 2.0 * d);

            double s = f + (dl + 412.0 * Math.Sin(2.0 * f) + 541.0 * Math.Sin(ls)) / Celestial.ARC;
            double h = f - 2 * d;
            double n = -526.0 * Math.Sin(h) + 44.0 * Math.Sin(l + h) - 31.0 * Math.Sin(-l + h) - 23.0 * Math.Sin(ls + h) +
                    11.0 * Math.Sin(-ls + h) - 25.0 * Math.Sin(-2.0 * l + f) + 21.0 * Math.Sin(-l + f);
            double l_moon = Celestial.TWOPI * Planar.Frac(lo + dl / 1296000.0); //in rad
            double b_moon = (18520.0 * Math.Sin(s) + n) / Celestial.ARC; //in rad
            double cb = Math.Cos(b_moon);
            double x = cb * Math.Cos(l_moon);
            double v = cb * Math.Sin(l_moon);
            double w = Math.Sin(b_moon);
            double y = Celestial.COSEPS * v - Celestial.SINEPS * w;
            double z = Celestial.SINEPS * v + Celestial.COSEPS * w;
            double rho = Math.Sqrt(1.0 - Math.Pow(z, 2));
            double dec = (360.0 / Celestial.TWOPI) * Math.Atan(z / rho);
            double ra = (48.0 / Celestial.TWOPI) * Math.Atan(y / (x + rho));
            if (ra < 0) { ra = ra + 24.0; }
            return (new Celestial.RADec(Transform.HoursToRadians(ra), Transform.DegreesToRadians(dec)));
        }

        public static double MaxAltitudeOld(DateTime DuskUTC, DateTime DawnUTC, Celestial.RADec position, Celestial.LatLon location)
        {
            //Computes the maximum altitude that a target at position (viewed from location) achieves 
            //   between Dusk and Dawn
            //Returns altitude in radians

            double tMaxAlt = 0;
            double tDuskHArad = position.HourAngle(DuskUTC, location);
            double tDawnHArad = position.HourAngle(DawnUTC, location);
            double tDuskAltrad = position.Altitude(tDuskHArad, location);
            double tDawnAltrad = position.Altitude(tDawnHArad, location);
            double darkHours = ((DawnUTC - DuskUTC).TotalHours); //in hours

            //Any hour angles that are more than nighttime (from meridian) will not peak.
            //or, better said, all hour angles that are <= nighttime will peak during the night.
            //if so, then the maximum altitude is at the hourangle 0.
            //if not, then the maximum altitude is the greater of the HA at dawn and the HA at dusk
            double tDuskHAHrs = Transform.RadiansToHours(tDuskHArad);
            double tDawnHAHrs = Transform.RadiansToHours(tDawnHArad);
            double tDuskAltHrs = Transform.RadiansToHours(tDuskAltrad);
            double tDawnAltHrs = Transform.RadiansToHours(tDawnAltrad);
            if (tDuskHAHrs >= (24 - darkHours))
            { tMaxAlt = position.Altitude(Transform.HoursToRadians(0.0), location); }
            else
            {
                if (tDuskAltHrs > tDawnAltHrs)
                { tMaxAlt = tDuskAltrad; }
                else
                { tMaxAlt = tDawnAltrad; }
            }
            return tMaxAlt;
        }

        public static double MaxAltitude(DateTime DuskUTC, DateTime DawnUTC, Celestial.RADec position, Celestial.LatLon location)
        {
            //Computes the maximum altitude that a target at position (viewed from location) achieves 
            //   between Dusk and Dawn
            //Returns altitude in radians

            //double tMaxAlt = 0;
            double tMaxHArad = Transform.HoursToRadians(0.0);
            double tDuskHArad = position.HourAngle(DuskUTC, location);
            double tDawnHArad = position.HourAngle(DawnUTC, location);
            double tDuskAltrad = position.Altitude(tDuskHArad, location);
            double tDawnAltrad = position.Altitude(tDawnHArad, location);
            double tMaxAltrad = position.Altitude(Transform.HoursToRadians(0.0), location);
            //double darkHours = ((DawnUTC - DuskUTC).TotalHours); //in hours
            //if HA = 0 is between dusk and dawn, then take the greater of the three altitudes
            //  otherwise take the greater of the dusk and dawn altitudes
            //  first pick the greater of the dusk and dawn altitudes
            bool isPeakAtNight = tMaxHArad > tDuskHArad && tMaxHArad < tDawnAltrad;
            if (tDuskAltrad > tDawnAltrad)  //dusk high
                if (tMaxAltrad > tDuskAltrad)  //dusk high, but HA0 higher
                    if (isPeakAtNight)
                        return tMaxAltrad; //dusk high, but HA0 higher and at night
                    else
                        return tDuskAltrad; //dusk high, but HA0 higher but not at night
                else
                    return tDuskAltrad; //dusk highest
            else //dawn high
                if (tMaxAltrad > tDawnAltrad)  //dawn high, but HA0 higher
                if (isPeakAtNight)
                    return tMaxAltrad; //dawn high, but HA0 higher and at night
                else
                    return tDawnAltrad; //dawn high, but HA0 higher but not at night
            else
                return tDawnAltrad; //dawn highest
        }

        public static double MinAltitude(DateTime DuskUTC, DateTime DawnUTC, Celestial.RADec position, Celestial.LatLon location)
        {
            //Computes the maximum altitude that a target at position (viewed from location) achieves 
            //   between Dusk and Dawn
            //Returns altitude in radians

            //double tMaxAlt = 0;
            double tMinHArad = Transform.HoursToRadians(12.0);
            double tDuskHArad = position.HourAngle(DuskUTC, location);
            double tDawnHArad = position.HourAngle(DawnUTC, location);
            double tDuskAltrad = position.Altitude(tDuskHArad, location);
            double tDawnAltrad = position.Altitude(tDawnHArad, location);
            double tLowestAltrad = position.Altitude(Transform.HoursToRadians(12.0), location);
            //if HA = 0 is between dusk and dawn, then take the greater of the three altitudes
            //  otherwise take the greater of the dusk and dawn altitudes
            //  first pick the greater of the dusk and dawn altitudes
            bool isLowAtNight = tMinHArad > tDuskHArad && tMinHArad < tDawnAltrad;
            if (tDuskAltrad < tDawnAltrad)  //dusk lower than dawn
                if (tLowestAltrad < tDuskAltrad)  //dusk low, but HA12 lower
                    if (isLowAtNight)
                        return tLowestAltrad; //dusk low, but HA12 lower and at night
                    else
                        return tDuskAltrad; //dusk low, but HA0 lower but not at night
                else
                    return tDuskAltrad; //dusk lowest
            else //dawn lower than dusk
                if (tLowestAltrad < tDawnAltrad)  //dawn low, but HA12 lower
                if (isLowAtNight)
                    return tLowestAltrad; //dawn low, but HA12 lower and at night
                else
                    return tDawnAltrad; //dawn low, but HA12 lower but not at night
            else
                return tDawnAltrad; //dawn lowest
        }

        #endregion
    }
}
