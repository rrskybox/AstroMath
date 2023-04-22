using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace AstroMath
{
    public static class Spherical
    {
        //Class to manipulate spherical coordinates
        public struct PointSph
        {
            public double Radius { get; set; }
            public double Phi { get; set; } //Radians
            public double Theta { get; set; } //Radians
        }

        public struct PointXYZ
        {
            public double X { get; set; }
            public double Y { get; set; }
            public double Z { get; set; }
        }

        public struct PointRangeRADec  //Radians
        {
            public double Range { get; set; }
            public double RA { get; set; } //Radians
            public double Dec { get; set; } //Radians
        }

        public static PointSph TranslateSpherical(PointSph sphPnt, PointSph sphNewOrigin)
        {
            //Translates sphOld spherical origin to sphNex spherical origin
            //convert points to cartesian
            PointXYZ oldLocXYZ = SphericalToCartesian(sphPnt);
            PointXYZ newOriginLocXYZ = SphericalToCartesian(sphNewOrigin);
            //Translate old point to new point
            PointXYZ newTransXYZ = new PointXYZ
            {
                X = oldLocXYZ.X - newOriginLocXYZ.X,
                Y = oldLocXYZ.Y - newOriginLocXYZ.Y,
                Z = oldLocXYZ.Z - newOriginLocXYZ.Z
            };
            //convert back to spherical
            PointSph topoLocSph = CartesianToSpherical(newTransXYZ);
            return topoLocSph;
        }

        public static PointRangeRADec TranslateRADecRate(PointRangeRADec preRangeRADec, double dRA, double dDec, PointSph sphNewOrigin)
        {
            //Computes translated (to new origin) difference of two spherical points to new difference
            //Create endpoint for pre and post with rate
            //convert points to cartesian
            PointXYZ preLocXYZ = RADecToCartesian(preRangeRADec);
            PointRangeRADec postRangeRADec = new PointRangeRADec
            {
                RA = preRangeRADec.RA + preRangeRADec.RA * dRA,
                Dec = preRangeRADec.Dec + preRangeRADec.Dec * dDec,
                Range = preRangeRADec.Range  //Assume small dRA,dDec
            };
            PointXYZ postLocXYZ = RADecToCartesian(postRangeRADec);
            PointSph preSph = CartesianToSpherical(preLocXYZ);
            PointSph postSph = CartesianToSpherical(postLocXYZ);
            //Translate origin for pre and post spherical coordinates
            PointSph transPreSph = TranslateSpherical(preSph, sphNewOrigin);
            PointSph transPostSph = TranslateSpherical(postSph, sphNewOrigin);
            //Convert back to RA/Dec/Range
            PointRangeRADec transPreRangeRADec = SphericalToRADec(transPreSph);
            PointRangeRADec transPostRangeRADec = SphericalToRADec(transPostSph);
            //Create new difference Range/RA/Dec
            PointRangeRADec transRateRangeRADec = new PointRangeRADec
            {
                RA = transPostRangeRADec.RA - transPreRangeRADec.RA,
                Dec = transPostRangeRADec.Dec - transPreRangeRADec.Dec,
                Range = transPostRangeRADec.Range - transPreRangeRADec.Range
            };
 
            return transRateRangeRADec;
        }

        public static PointSph CartesianToSpherical(PointXYZ xyz)
        {
            //Converts Cartesian XYZ to Spherical coordinates
            PointSph conv = new PointSph();
            conv.Radius = Math.Sqrt(Math.Pow(xyz.X, 2) + Math.Pow(xyz.Y, 2) + Math.Pow(xyz.Z, 2));
            if (conv.Radius != 0)
                conv.Phi = Math.Acos(xyz.Z / conv.Radius);
            else
                conv.Phi = 0;
            conv.Theta = Math.Atan2(xyz.Y, xyz.X);
            return conv;
        }

        public static PointXYZ SphericalToCartesian(PointSph ria)
        {
            //Converts Spherical coordinates to Cartesian XYZ
            PointXYZ conv = new PointXYZ();
            conv.X = ria.Radius * Math.Sin(ria.Phi) * Math.Cos(ria.Theta);
            conv.Y = ria.Radius * Math.Sin(ria.Phi) * Math.Sin(ria.Theta);
            conv.Z = ria.Radius * Math.Cos(ria.Phi);
            return conv;
        }

        public static PointSph RADecToSpherical(PointRangeRADec rad)
        {
            //converts RA/Dec coordinates to Spherical coordinates
            PointSph pntSph = new PointSph();
            pntSph.Radius = rad.Range;
            pntSph.Phi = (Celestial.PI / 2) - rad.Dec;
            pntSph.Theta = Celestial.PI - rad.RA;
            return pntSph;
        }

        public static PointRangeRADec SphericalToRADec(PointSph pntSph)
        {
            //converts RA/Dec coordinates to Spherical coordinates
            PointRangeRADec pntRD = new PointRangeRADec();
            pntRD.Range = pntSph.Radius;
            pntRD.Dec = (Celestial.PI / 2) - pntSph.Phi;
            pntRD.RA = Celestial.PI - pntSph.Theta;
            return pntRD;
        }

        public static PointXYZ RADecToCartesian(PointRangeRADec rad)
        {
            //converts RA/Dec coordinates to Spherical coordinates
            return SphericalToCartesian(RADecToSpherical(rad));
        }

        public static PointRangeRADec CartesianToRangeRADec(PointXYZ pnt)
        {
            //converts RA/Dec coordinates to Spherical coordinates
            return SphericalToRADec(CartesianToSpherical(pnt));
        }

    }
}
