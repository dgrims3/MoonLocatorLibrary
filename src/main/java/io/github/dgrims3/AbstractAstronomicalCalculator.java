package io.github.dgrims3;

import java.time.LocalDate;
import java.time.LocalDateTime;
import java.time.ZoneId;
import java.time.ZonedDateTime;

public abstract class AbstractAstronomicalCalculator {
    public static final double PI = 3.1415926535897932384626433832795028841971;

    /**
     * Converts radians to degrees.
     *
     * @param radians The angle in radians
     * @return The angle in degrees
     */
    public double radToDeg(double radians) {
        return (180.0 / PI) * radians;
    }

    /**
     * Converts degrees to radians.
     *
     * @param degrees The angle in degrees
     * @return The angle in radians
     */
    public double degToRad(double degrees) {
        return (PI / 180.0) * degrees;
    }

    /**
     * Limits an angle to the range [0, 360) degrees.
     * Normalizes any angle value to be within a single rotation.
     *
     * @param degrees The angle in degrees (can be any value)
     * @return The normalized angle in the range [0, 360) degrees
     */
    public double limit_degrees(double degrees) {
        double limited;
        degrees /= 360.0;
        limited = 360.0 * (degrees - Math.floor(degrees));
        if (limited < 0) limited += 360.0;

        return limited;
    }

    /**
     * Calculates the obliquity of the ecliptic (tilt of Earth's axis) including nutation correction.
     * This is the angle between Earth's equatorial plane and the ecliptic plane.
     *
     * @param jce Julian centuries from J2000.0
     * @return The obliquity of the ecliptic in degrees
     */
    public double obliquityOfEcliptic(double jce) {
        double e0 = calcMeanObliquityOfEcliptic(jce);
        double omega = 125.04452 - 1934.136261 * jce;
        return e0 + 0.00256 * Math.cos(degToRad(omega));
    }

    /**
     * Calculates the mean obliquity of the ecliptic without nutation correction.
     *
     * @param jce Julian centuries from J2000.0
     * @return The mean obliquity in degrees
     */
    public double calcMeanObliquityOfEcliptic(double jce) {
        double seconds = 21.448 - jce * (46.8150 + jce * (0.00059 - jce * (0.001813)));
        return 23.0 + (26.0 + (seconds / 60.0)) / 60.0;
    }

    /**
     * Calculates the longitude of the Moon's ascending node.
     * The ascending node is where the Moon crosses the ecliptic plane moving northward.
     *
     * @param jce Julian centuries from J2000.0
     * @return The longitude of the ascending node in degrees
     */
    public double moon_ascending_node(double jce) {
        return limit_degrees((125.04452 - (1934.136261 * jce) + (.0020708 * Math.pow(jce, 2)) + (Math.pow(jce, 3) / 450000)));
    }

    /**
     * Calculates the local hour angle for a celestial object.
     * The hour angle measures the time since the object crossed the observer's meridian.
     *
     * @param siderealInDegrees Greenwich sidereal time in degrees
     * @param observerLongitude Observer's longitude in degrees (positive east)
     * @param ascension Right ascension of the celestial object in degrees
     * @return The local hour angle in degrees
     */
    public double localHourAngle(double siderealInDegrees, double observerLongitude, double ascension) {
        double hourAngle = siderealInDegrees - observerLongitude - ascension;
        if (hourAngle < 0) {
            return hourAngle + 360;
        }
        return hourAngle;
    }

    /**
     * Calculates the altitude (elevation) of a celestial object above the horizon.
     *
     * @param localHourAngle The local hour angle in degrees
     * @param observersLatitude Observer's latitude in degrees
     * @param declination Declination of the celestial object in degrees
     * @return The altitude in degrees (positive above horizon, negative below)
     */
    public double localAltitude(double localHourAngle, double observersLatitude, double declination) {
        double altitude = Math.asin((Math.sin(degToRad(observersLatitude)) * Math.sin(degToRad(declination))) + (Math.cos(degToRad(observersLatitude)) * Math.cos(degToRad(declination)) * Math.cos(degToRad(localHourAngle))));
        return radToDeg(altitude);
    }

    /**
     * Calculates the azimuth (compass direction) of a celestial object.
     *
     * @param localHourAngle The local hour angle in degrees
     * @param observersLatitude Observer's latitude in degrees
     * @param declination Declination of the celestial object in degrees
     * @return The azimuth in degrees (measured from north, positive eastward)
     */
    public double azimuth(double localHourAngle, double observersLatitude, double declination) {
        double numerator = Math.sin(degToRad(localHourAngle));
        double denominator = ((Math.cos(degToRad(localHourAngle)) * Math.sin(degToRad(observersLatitude))) - (Math.tan(degToRad(declination)) * Math.cos(degToRad(observersLatitude))));
        double azimuth = Math.atan2(numerator, denominator);
        return radToDeg(azimuth);
    }

    /**
     * Calculates the declination from ecliptic coordinates.
     * Converts ecliptic latitude and longitude to equatorial declination.
     *
     * @param lat Ecliptic latitude in degrees
     * @param lng Ecliptic longitude in degrees
     * @param jce Julian centuries from J2000.0
     * @return The declination in degrees
     */
    public double getDeclination(double lat, double lng, double jce) {
        double obOfEclip = degToRad(obliquityOfEcliptic(jce));
        lat = degToRad(lat);
        lng = degToRad(lng);
        return radToDeg(Math.asin((Math.sin(lat) * Math.cos(obOfEclip)) + (Math.cos(lat) * Math.sin(obOfEclip) * Math.sin(lng))));
    }

    /**
     * Calculates the right ascension from ecliptic coordinates.
     * Converts ecliptic latitude and longitude to equatorial right ascension.
     *
     * @param lat Ecliptic latitude in degrees
     * @param lng Ecliptic longitude in degrees
     * @param jce Julian centuries from J2000.0
     * @return The right ascension in degrees
     */
    public double getRightAscension(double lat, double lng, double jce) {
        double obOfEclip = degToRad(obliquityOfEcliptic(jce));
        lat = degToRad(lat);
        lng = degToRad(lng);
        //  α=atan2(cosβsinλcosε−sinβsinε,cosβcosλ)
        double inpt1 = (Math.cos(lat) * Math.sin(lng) * Math.cos(obOfEclip)) - (Math.sin(lat) * Math.sin(obOfEclip));
        double inpt2 = Math.cos(lat) * Math.cos(lng);
        return radToDeg(Math.atan2(inpt1, inpt2));
    }

    /**
     * Normalizes a decimal value representing a fraction of a day to the range [0, 1).
     *
     * @param decimal The decimal value (typically a fraction of a day)
     * @return The normalized value in the range [0, 1)
     */
    public double getInRange(double decimal) {
        if (decimal < 0) {
            return decimal + 1;
        } else if (decimal > 1) {
            return decimal - 1;
        } else {
            return decimal;
        }
    }

    /**
     * Calculates the nutation in longitude using the IAU 1980 nutation model.
     * Nutation is a periodic oscillation in Earth's axis of rotation.
     *
     * @param jce Julian centuries from J2000.0
     * @return The nutation in longitude in degrees
     */
    public double nutation_in_longitude(double jce) {
        double T2 = Math.pow(jce, 2);
        double T3 = Math.pow(jce, 3);

        double D = degToRad(297.85036 + 445267.11148 * jce - 0.0019142 * T2 + (T3 / 189474));
        double M = degToRad(357.52772 + 35999.05034 * jce - 0.0001603 * T2 - (T3 / 300000));
        double MPR = degToRad(134.96298 + 477198.867398 * jce + 0.0086972 * T2 + (T3 / 56250));
        double F = degToRad(93.27191 + 483202.017538 * jce - 0.0036825 * T2 + (T3 / 327270));
        double omega = degToRad(moon_ascending_node(jce));

        double w = Math.sin(omega) * (-174.2 * jce - 171996);
        w = w + Math.sin(2 * (F + omega - D)) * (-1.6 * jce - 13187);
        w = w + Math.sin(2 * (F + omega)) * (-2274 - 0.2 * jce);
        w = w + Math.sin(2 * omega) * (0.2 * jce + 2062);
        w = w + Math.sin(M) * (1426 - 3.4 * jce);
        w = w + Math.sin(MPR) * (0.1 * jce + 712);
        w = w + Math.sin(2 * (F + omega - D) + M) * (1.2 * jce - 517);
        w = w + Math.sin(2 * F + omega) * (-0.4 * jce - 386);
        w = w + Math.sin(2 * (F + omega - D) - M) * (217 - 0.5 * jce);
        w = w + Math.sin(2 * (F - D) + omega) * (129 + 0.1 * jce);
        w = w + Math.sin(MPR + omega) * (0.1 * jce + 63);
        w = w + Math.sin(omega - MPR) * (-0.1 * jce - 58);
        w = w + Math.sin(2 * M) * (17 - 0.1 * jce);
        w = w + Math.sin(2 * (M + F + omega - D)) * (0.1 * jce - 16);
        w = w - 301 * Math.sin(2 * (F + omega) + MPR);
        w = w - 158 * Math.sin(MPR - 2 * D);
        w = w + 123 * Math.sin(2 * (F + omega) - MPR);
        w = w + 63 * Math.sin(2 * D);
        w = w - 59 * Math.sin(2 * (D + F + omega) - MPR);
        w = w - 51 * Math.sin(2 * F + MPR + omega);
        w = w + 48 * Math.sin(2 * (MPR - D));
        w = w + 46 * Math.sin(2 * (F - MPR) + omega);
        w = w - 38 * Math.sin(2 * (D + F + omega));
        w = w - 31 * Math.sin(2 * (MPR + F + omega));
        w = w + 29 * Math.sin(2 * MPR);
        w = w + 29 * Math.sin(2 * (F + omega - D) + MPR);
        w = w + 26 * Math.sin(2 * F);
        w = w - 22 * Math.sin(2 * (F - D));
        w = w + 21 * Math.sin(2 * F + omega - MPR);
        w = w + 16 * Math.sin(2 * D - MPR + omega);
        w = w - 15 * Math.sin(M + omega);
        w = w - 13 * Math.sin(MPR + omega - 2 * D);
        w = w - 12 * Math.sin(omega - M);
        w = w + 11 * Math.sin(2 * (MPR - F));
        w = w - 10 * Math.sin(2 * (F + D) + omega - MPR);
        w = w - 8 * Math.sin(2 * (F + D + omega) + MPR);
        w = w + 7 * Math.sin(2 * (F + omega) + M);
        w = w - 7 * Math.sin(MPR - 2 * D + M);
        w = w - 7 * Math.sin(2 * (F + omega) - M);
        w = w - 7 * Math.sin(2 * D + 2 * F + omega);
        w = w + 6 * Math.sin(2 * D + MPR);
        w = w + 6 * Math.sin(2 * (MPR + F + omega - D));
        w = w + 6 * Math.sin(2 * (F - D) + MPR + omega);
        w = w - 6 * Math.sin(2 * (D - MPR) + omega);
        w = w - 6 * Math.sin(2 * D + omega);
        w = w + 5 * Math.sin(MPR - M);
        w = w - 5 * Math.sin(2 * (F - D) + omega - M);
        w = w - 5 * Math.sin(omega - 2 * D);
        w = w - 5 * Math.sin(2 * (MPR + F) + omega);
        w = w + 4 * Math.sin(2 * (MPR - D) + omega);
        w = w + 4 * Math.sin(2 * (F - D) + M + omega);
        w = w + 4 * Math.sin(MPR - 2 * F);
        w = w - 4 * Math.sin(MPR - D);
        w = w - 4 * Math.sin(M - 2 * D);
        w = w - 4 * Math.sin(D);
        w = w + 3 * Math.sin(2 * F + MPR);
        w = w - 3 * Math.sin(2 * (F + omega - MPR));
        w = w - 3 * Math.sin(MPR - D - M);
        w = w - 3 * Math.sin(M + MPR);
        w = w - 3 * Math.sin(2 * (F + omega) + MPR - M);
        w = w - 3 * Math.sin(2 * (D + F + omega) - M - MPR);
        w = w - 3 * Math.sin(2 * (F + omega) + 3 * MPR);
        w = w - 3 * Math.sin(2 * (D + F + omega) - M);
        w = radToDeg(w);
        return w / 36000000.0;
    }

    /**
     * Calculates Julian centuries from the J2000.0 epoch.
     *
     * @param jd Julian day number
     * @return Julian centuries since J2000.0
     */
    public double calcTimeJulianCent(double jd) {
        return (jd - 2451545.0) / 36525.0;
    }

    /**
     * Converts a calendar date to Julian day number.
     *
     * @param year The year
     * @param month The month (1-12)
     * @param day The day of the month (can include fractional day)
     * @return The Julian day number
     */
    public double getJDFromCalenderDate(int year, int month, double day) {
        if (month <= 2) {
            year -= 1;
            month += 12;
        }
        double A = Math.floor(year / 100.0);
        double B = 2 - A + Math.floor(A / 4);
        return Math.floor(365.25 * (year + 4716)) + Math.floor(30.6001 * (month + 1)) + day + B - 1524.5;
    }

    /**
     * Converts a Julian day number to a calendar date.
     *
     * @param jd The Julian day number
     * @return The calendar date as a LocalDate
     */
    public LocalDate getCalendarDateFromJD(double jd) {
        jd += .5;
        double A;
        double B;
        double C;
        double D;
        double E;
        double F = jd % 1;
        double Z = Math.floor(jd);
        if (Z >= 2299161) {
            double val = Math.floor(((Z - 1867216.25) / 36524.25));
            A = Math.floor(Z + 1 + val - Math.floor((val / 4)));
        } else {
            A = Z;
        }
        B = A + 1524;
        C = Math.floor(((B - 122.1) / 365.25));
        D = Math.floor(365.25 * C);
        E = Math.floor((B - D) / 30.6001);

        int Day = (int) Math.floor(B - D - (30.6001 * E) + F);
        int Month = (int) Math.floor(E < 14 ? E - 1 : E - 13);
        int Year = (int) Math.floor(Month > 2 ? C - 4716 : C - 4715);
        return LocalDate.of(Year, Month, Day);
    }

    /**
     * Calculates Greenwich Mean Sidereal Time at 0h UT.
     * Sidereal time measures Earth's rotation relative to distant stars.
     *
     * @param jce Julian centuries from J2000.0
     * @return Greenwich mean sidereal time in degrees
     */
    public double greenwichMeanSiderealTime(double jce) {
        return limit_degrees(((((float) 1 / 38710000) * jce + .000387933) * jce + 36000.770053608) * jce + 100.46061837);
    }

    /**
     * Calculates Greenwich Apparent Sidereal Time including nutation correction.
     *
     * @param jce Julian centuries from J2000.0
     * @return Greenwich apparent sidereal time in degrees
     */
    public double greenwichApparentSiderealTime(double jce) {
        double nutation = degToRad(nutation_in_longitude(jce));
        double obOfEcliptic = degToRad(obliquityOfEcliptic(jce));
        double meanSidereal = greenwichMeanSiderealTime(jce);
        double equationOfEquinox = radToDeg(nutation * (Math.cos(obOfEcliptic)));
        return meanSidereal + equationOfEquinox;
    }

    /**
     * Calculates sidereal time at a specific instant during the day.
     *
     * @param siderealTime Sidereal time at 0h UT in degrees
     * @param fractionOfDay Fraction of the day (0.0 = midnight, 0.5 = noon)
     * @return Sidereal time at the instant in degrees
     */
    public double siderealTimeAtInstantAtGreenwichInDegrees(double siderealTime, double fractionOfDay) {
        return limit_degrees(siderealTime + (360.985647 * fractionOfDay));
    }

    /**
     * Converts degrees to hours, minutes, and seconds format.
     * Commonly used for converting right ascension from degrees to time units.
     *
     * @param degrees The angle in degrees
     * @return Array containing [hours, minutes, seconds]
     */
    public double[] degreesToHoursMinutesSeconds(double degrees) {
        double hours = Math.floor(degrees / 15);
        double minutes = Math.floor(((degrees / 15) - hours) * 60);
        double seconds = ((((degrees / 15) - hours) * 60) - minutes) * 60;
        // return String.format(Locale.getDefault(), "%f, %f, %f", hours, minutes, seconds);
        return new double[]{hours, minutes, seconds};
    }

    /**
     * Converts a LocalDateTime to a fraction of a day.
     *
     * @param localDateTime The local date and time
     * @return Fraction of the day (0.0 = midnight, 0.5 = noon)
     */
    public double localDateTimeToFractionOfDay(LocalDateTime localDateTime) {
        return (localDateTime.getHour() + ((localDateTime.getMinute() + (localDateTime.getSecond() / 60.0)) / 60.0)) / 24;
    }

    /**
     * Converts a UTC fraction of day to local ZonedDateTime.
     *
     * @param fractionOfDay Fraction of the day in UTC
     * @param date The date
     * @return The corresponding ZonedDateTime in the system's default time zone
     */
    public ZonedDateTime localFractionOfDayFromUTCToLocal(double fractionOfDay, LocalDate date) {
        ZoneId zoneId = ZoneId.systemDefault();

        double doubleHour = fractionOfDay * 24;
        int hour = (int) Math.floor(doubleHour);
        double doubleMinute = 60 * (doubleHour - hour);
        int minute = (int) Math.floor(doubleMinute);
        ZonedDateTime dateTimeInUTC = LocalDateTime.of(date.getYear(), date.getMonthValue(), date.getDayOfMonth(), hour, minute).atZone(ZoneId.of("UTC"));

        return dateTimeInUTC.withZoneSameInstant(zoneId);
    }

    /**
     * Converts a UTC fraction of day to ZonedDateTime in a specified time zone.
     *
     * @param fractionOfDay Fraction of the day in UTC
     * @param date The date
     * @param zoneId The time zone to convert to
     * @return The corresponding ZonedDateTime in the specified time zone
     */
    public ZonedDateTime localFractionOfDayFromUTCToLocal(double fractionOfDay, LocalDate date, ZoneId zoneId) {
        double doubleHour = fractionOfDay * 24;
        int hour = (int) Math.floor(doubleHour);
        double doubleMinute = 60 * (doubleHour - hour);
        int minute = (int) Math.floor(doubleMinute);
        ZonedDateTime dateTimeInUTC = LocalDateTime.of(date.getYear(), date.getMonthValue(), date.getDayOfMonth(), hour, minute).atZone(ZoneId.of("UTC"));

        return dateTimeInUTC.withZoneSameInstant(zoneId);
    }
}
