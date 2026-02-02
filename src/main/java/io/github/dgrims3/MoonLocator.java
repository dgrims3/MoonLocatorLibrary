package io.github.dgrims3;

import java.time.LocalDateTime;
import java.time.ZoneId;
import java.time.ZonedDateTime;
import java.util.List;

/**
 * References:
 *     - Reda, I. (2010). Solar Eclipse Monitoring for Solar Energy Applications
 *       Using the Solar and Moon Position Algorithms. NREL/TP-3B0-47681.
 *     - Meeus, J. (1998). Astronomical Algorithms. 2nd Ed.
 */
public class MoonLocator extends MoonCalculator {
    private final int offSet;
    private final ZonedDateTime dateTime;
    private final double latitude;
    private final double longitude;

    /**
     * Constructs a MoonLocator for calculating moon positions and properties at a specific location and time.
     *
     * @param dateTime The zoned date and time for moon calculations
     * @param latitude The observer's latitude in degrees (positive north, negative south)
     * @param longitude The observer's longitude in degrees (positive east, negative west)
     */
    public MoonLocator(ZonedDateTime dateTime, double latitude, double longitude) {
        this.dateTime = dateTime;
        this.latitude = latitude;
        this.longitude = longitude;
        this.offSet = dateTime.getOffset().getTotalSeconds() / 3600;
    }

    /**
     * Constructs a MoonLocator using the system's default time zone.
     *
     * @param dateTime The local date and time for moon calculations
     * @param latitude The observer's latitude in degrees (positive north, negative south)
     * @param longitude The observer's longitude in degrees (positive east, negative west)
     */
    public MoonLocator(LocalDateTime dateTime, double latitude, double longitude) {
        this(dateTime.atZone(ZoneId.systemDefault()), latitude, longitude);
    }

    /**
     * Constructs a MoonLocator with a specific time zone.
     *
     * @param dateTime The local date and time for moon calculations
     * @param latitude The observer's latitude in degrees (positive north, negative south)
     * @param longitude The observer's longitude in degrees (positive east, negative west)
     * @param zoneId The time zone for the given dateTime
     */
    public MoonLocator(LocalDateTime dateTime, double latitude, double longitude, ZoneId zoneId) {
        this(dateTime.atZone(zoneId), latitude, longitude);
    }

    /**
     * Gets the time zone offset from UTC in hours.
     *
     * @return The UTC offset in hours
     */
    public int getOffSet() {
        return offSet;
    }

    /**
     * Gets the zoned date and time configured for this MoonLocator.
     *
     * @return The zoned date and time
     */
    public ZonedDateTime getDateTime() {
        return dateTime;
    }

    /**
     * Gets the time zone configured for this MoonLocator.
     *
     * @return The time zone as a ZoneId
     */
    public ZoneId getZoneId() {
        return dateTime.getZone();
    }

    /**
     * Gets the observer's latitude.
     *
     * @return The latitude in degrees (positive north, negative south)
     */
    public double getLatitude() {
        return latitude;
    }

    /**
     * Gets the observer's longitude.
     *
     * @return The longitude in degrees (positive east, negative west)
     */
    public double getLongitude() {
        return longitude;
    }

    /**
     * Calculates the moon rise time for the configured date and location.
     * Moon rise is when the moon crosses the horizon ascending.
     *
     * @return The ZonedDateTime when the moon rises, in the configured time zone
     */
    public ZonedDateTime getMoonRise() {
        List<ZonedDateTime> riseTransitSet = moonRisingSettingTransitPrecise(latitude, longitude, dateTime.toLocalDate());
        return riseTransitSet.get(time.Rise.ordinal()).withZoneSameInstant(dateTime.getZone());
    }

    /**
     * Calculates the moon set time for the configured date and location.
     * Moon set is when the moon crosses the horizon descending.
     *
     * @return The ZonedDateTime when the moon sets, in the configured time zone
     */
    public ZonedDateTime getMoonSet() {
        List<ZonedDateTime> riseTransitSet = moonRisingSettingTransitPrecise(latitude, longitude, dateTime.toLocalDate());
        return riseTransitSet.get(time.Set.ordinal()).withZoneSameInstant(dateTime.getZone());
    }

    /**
     * Calculates the moon transit time for the configured date and location.
     * Moon transit is when the moon reaches its highest point in the sky (culmination).
     *
     * @return The ZonedDateTime when the moon transits, in the configured time zone
     */
    public ZonedDateTime getMoonTransit() {
        List<ZonedDateTime> riseTransitSet = moonRisingSettingTransitPrecise(latitude, longitude, dateTime.toLocalDate());
        return riseTransitSet.get(time.Transit.ordinal()).withZoneSameInstant(dateTime.getZone());
    }

    /**
     * Calculates the moon's azimuth at the configured date, time, and location.
     * Azimuth is the compass direction from north (0°) measured clockwise.
     *
     * @return The azimuth in degrees (0° = north, 90° = east, 180° = south, 270° = west)
     */
    public double getAzimuth() {
        LocalDateTime localDateTime = dateTime.toLocalDateTime();
        double[] azimuthAltitude = getAzimuthAndAltitudeForMoonAtInstant(localDateTime, offSet, latitude, longitude);
        return azimuthAltitude[0];
    }

    /**
     * Calculates the moon's altitude at the configured date, time, and location.
     * Altitude is the angular height above the horizon.
     *
     * @return The altitude in degrees (0° = horizon, 90° = zenith, negative = below horizon)
     */
    public double getAltitude() {
        LocalDateTime localDateTime = dateTime.toLocalDateTime();
        double[] azimuthAltitude = getAzimuthAndAltitudeForMoonAtInstant(localDateTime, offSet, latitude, longitude);
        return azimuthAltitude[1];
    }

    /**
     * Calculates the Earth-Moon distance at the configured date and time.
     *
     * @return The distance from Earth to Moon in kilometers
     */
    public double getDistance() {
        LocalDateTime utcDateTime = dateTime.toLocalDateTime().minusHours(offSet);
        double jd = getJDFromCalenderDate(utcDateTime.getYear(), utcDateTime.getMonthValue(), utcDateTime.getDayOfMonth());
        double fractionOfDay = localDateTimeToFractionOfDay(utcDateTime);
        double[] ascDecDistance = lunarAscensionDeclinationDistance(jd + fractionOfDay);
        return ascDecDistance[coords.DISTANCE.ordinal()];
    }

    /**
     * Calculates the moon's illuminated fraction at the configured date and time.
     * This represents how much of the moon's visible disk is illuminated by the sun.
     *
     * @return The illumination percentage (0 = new moon, 100 = full moon)
     */
    public double getIllumination() {
        LocalDateTime utcDateTime = dateTime.toLocalDateTime().minusHours(offSet);
        double jd = getJDFromCalenderDate(utcDateTime.getYear(), utcDateTime.getMonthValue(), utcDateTime.getDayOfMonth());
        double fractionOfDay = localDateTimeToFractionOfDay(utcDateTime);
        return moonIlluminatedFractionOfDisk(jd + fractionOfDay);
    }

    /**
     * Calculates the moon's right ascension at the configured date and time.
     * Right ascension is the celestial equivalent of longitude in the equatorial coordinate system.
     *
     * @return The right ascension in degrees
     */
    public double getRightAscension() {
        LocalDateTime utcDateTime = dateTime.toLocalDateTime().minusHours(offSet);
        double jd = getJDFromCalenderDate(utcDateTime.getYear(), utcDateTime.getMonthValue(), utcDateTime.getDayOfMonth());
        double fractionOfDay = localDateTimeToFractionOfDay(utcDateTime);
        double[] ascDecDistance = lunarAscensionDeclinationDistance(jd + fractionOfDay);
        return ascDecDistance[position.ASCENSION.ordinal()];
    }

    /**
     * Calculates the moon's declination at the configured date and time.
     * Declination is the celestial equivalent of latitude in the equatorial coordinate system.
     *
     * @return The declination in degrees (positive north of celestial equator, negative south)
     */
    public double getDeclination() {
        LocalDateTime utcDateTime = dateTime.toLocalDateTime().minusHours(offSet);
        double jd = getJDFromCalenderDate(utcDateTime.getYear(), utcDateTime.getMonthValue(), utcDateTime.getDayOfMonth());
        double fractionOfDay = localDateTimeToFractionOfDay(utcDateTime);
        double[] ascDecDistance = lunarAscensionDeclinationDistance(jd + fractionOfDay);
        return ascDecDistance[position.DECLINATION.ordinal()];
    }

    /**
     * Calculates the moon's local hour angle at the configured date, time, and location.
     * The hour angle indicates the moon's position relative to the observer's meridian.
     *
     * @return The local hour angle in degrees (0° = on meridian, increases westward)
     */
    public double getLocalHourAngle() {
        LocalDateTime utcDateTime = dateTime.toLocalDateTime().minusHours(offSet);
        double fractionOfDay = localDateTimeToFractionOfDay(utcDateTime);
        double jd = getJDFromCalenderDate(utcDateTime.getYear(), utcDateTime.getMonthValue(), utcDateTime.getDayOfMonth());
        double jce = calcTimeJulianCent(jd);
        double sidereal = greenwichApparentSiderealTime(jce);
        double siderealAtInstant = siderealTimeAtInstantAtGreenwichInDegrees(sidereal, fractionOfDay);
        double[] ascDec = lunarAscensionDeclinationDistance(jd + fractionOfDay);
        return localHourAngle(siderealAtInstant, -1 * longitude, ascDec[position.ASCENSION.ordinal()]);
    }

    /**
     * Calculates the moon's geographic latitude at the configured date and time.
     * This is the latitude of the point on Earth where the moon appears directly overhead
     * (the moon's sub-point). The geographic latitude equals the moon's declination.
     *
     * @return The moon's geographic latitude in degrees (positive north, negative south)
     */
    public double getMoonGeographicLatitude() {
        LocalDateTime utcDateTime = dateTime.toLocalDateTime().minusHours(offSet);
        double fractionOfDay = localDateTimeToFractionOfDay(utcDateTime);
        double jd = getJDFromCalenderDate(utcDateTime.getYear(), utcDateTime.getMonthValue(), utcDateTime.getDayOfMonth());
        return getMoonLatAtInstant(jd, fractionOfDay);
    }

    /**
     * Calculates the moon's geographic longitude at the configured date and time.
     * This is the longitude of the point on Earth where the moon appears directly overhead
     * (the moon's sub-point). The longitude is determined by the moon's right ascension
     * and Greenwich Mean Sidereal Time.
     *
     * @return The moon's geographic longitude in degrees (positive east, negative west)
     */
    public double getMoonGeographicLongitude() {
        LocalDateTime utcDateTime = dateTime.toLocalDateTime().minusHours(offSet);
        double fractionOfDay = localDateTimeToFractionOfDay(utcDateTime);
        double jd = getJDFromCalenderDate(utcDateTime.getYear(), utcDateTime.getMonthValue(), utcDateTime.getDayOfMonth());
        return getMoonLongAtInstant(jd, fractionOfDay);
    }
}
